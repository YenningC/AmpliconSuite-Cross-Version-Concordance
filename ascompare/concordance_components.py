"""Separate complete AA reconstruction from AC classification and cache each.

This never enables upstream AA's reuse_intermediate_files: no segment, edge or
coverage estimate from one AA commit is substituted into another AA commit.
"""
import copy
import json
from pathlib import Path
import re

from .common import atomic_json, digest, now
from . import concordance_cache as cache
from .concordance_component_identity import component_identity, effective_environment


def _files(attempt, component):
    rows = [p.relative_to(attempt).as_posix() for p in (attempt / component).rglob('*')
            if p.is_file() and not p.is_symlink()]
    names = [component + '.log']
    if component == 'aa':
        names += ['aa_background.json', 'selected_seeds.bed', 'reference/coverage.stats']
        names += [p.name for p in attempt.glob('bam_access.*.json')]
    return rows + [n for n in names if (attempt / n).is_file() and not (attempt / n).is_symlink()]


def _aa_input_files(attempt):
    return {p.name: digest(p) for p in sorted((attempt / 'aa').iterdir())
            if p.is_file() and p.name.endswith(('_summary.txt', '_graph.txt', '_cycles.txt'))}


def _validate_aa(attempt, lengths):
    """Validate native AA structure before any expensive AC call or publication."""
    from .parsing import aa_amplicon_count, _amplicon_files, _read_graph_pair
    graphs = _amplicon_files(sorted((attempt / 'aa').glob('*_amplicon*_graph.txt')), 'graph')
    cycles = _amplicon_files(sorted((attempt / 'aa').glob('*_amplicon*_cycles.txt')), 'cycles')
    count = aa_amplicon_count(attempt)
    if count is None or count != len(graphs) or count != len(cycles) or graphs.keys() != cycles.keys():
        raise ValueError('Incomplete AA output: summary must match all graph/cycles pairs')
    for amp in graphs:
        _read_graph_pair(graphs[amp], cycles[amp], lengths, ())
    return count


def _staged(entry, metadata):
    stage = copy.deepcopy(entry['stage'])
    stage.update(cache_hit=True, cache_key=metadata['key'],
                 producer_run_dir=metadata['producer_run_dir'],
                 execution='reused_verified_artifacts')
    return stage


def _new_attempt(api, out):
    # Atomic mkdir closes the standalone-caller race; workflow calls additionally
    # hold the existing whole-output lock.
    while True:
        attempt = api._new_attempt(out, create=False)
        try:
            attempt.mkdir()
            return attempt
        except FileExistsError:
            continue


def run_component_downstream(sample, seeds, runtime, entry, verified, config, out, measurement=None):
    from . import concordance as api
    from .parsing import parse_run, aa_amplicon_count
    from .concordance_timing import measure
    from .versions import _check_source

    out = Path(out)
    supported = entry.get('supported_labels', entry.get('supported_classes', []))
    sid, vid = sample['sample_id'], runtime['id']
    stage, stages, component_cache = 'input_check', {}, {}
    attempt = None
    record = None
    cache_root = Path(config.get('_component_cache_dir', out / '.component_cache'))
    parameter_policy = config.get('parameter_policy', 'custom')
    official = parameter_policy == 'suite_defaults'
    aa_env = effective_environment('aa', runtime, config['reference'], pin_hash_seed=not official)
    ac_env = effective_environment('ac', runtime, config['reference'])
    hash_seed = aa_env.get('PYTHONHASHSEED')
    background = measurement or dict(mode='native_measurement', reason='No eligible same-version seed-stage measurement')
    if (not official and background['mode'] == 'same_version_seed_stage'
            and any(c.isspace() for c in str(Path(sample['alignment']).resolve()))
            and any(c.isspace() for c in str(out.absolute()))):
        background = dict(mode='native_measurement',
            reason='BAM and working paths both contain whitespace; native AA cannot match a private cache alias')
    reuse_stats = background['mode'] == 'same_version_seed_stage'
    try:
        fixed_hash = bool(hash_seed and re.fullmatch(r'[0-9]{1,10}', hash_seed)
                          and int(hash_seed) <= 4294967295)
        if not official and not fixed_hash:
            raise ValueError('Bounded AA requires a numeric PYTHONHASHSEED between 0 and 4294967295 for reproducible downsampling')
        # Launch/source verification is mandatory even when every output is cached.
        if not verified.get('ready'):
            raise ValueError('Runtime source/version checks failed: ' + '; '.join(verified.get('diagnostics', [])))
        for name in ('aa', 'ac'):
            if verified.get(name, {}).get('source_root'):
                _check_source(verified[name]['source_root'], entry[name])
        policy = json.loads(Path(config['_io_policy_path']).read_text())
        # Refresh only the bounded BAM identity (three small blocks) and index,
        # so a long multi-version run cannot reuse outputs after its input moves.
        current_bam = api.fast_fingerprint(sample['alignment'])
        approved = {str(Path(p).resolve()): str(Path(i).resolve())
                    for p, i in policy.get('alignments', {}).items()}
        index = approved.get(current_bam['path'])
        if not index:
            raise ValueError('AA BAM is absent from the verified access policy')
        current_index = dict(path=index, sha256=digest(index), method='index_file_full_SHA256')
        for observed in (current_bam, current_index):
            previous = next((p for p in sample['_input_identities'] if p.get('path') == observed['path']), None)
            if previous and any(previous.get(k) != observed.get(k) for k in
                                ('bytes', 'mtime_ns', 'sampled_sha256', 'sha256') if k in previous):
                raise ValueError('BAM/index identity changed since bounded input validation')
        if reuse_stats:
            from .concordance_aa_stats import validated_values
            if (digest(background['path']) != background.get('sha256') or
                    validated_values(Path(background['path']), Path(sample['alignment']).resolve()) != background.get('fields')):
                raise ValueError('Native AA measurement changed after the seed stage')
        params = config.get('aa_parameters', {})
        runmode = 'FULL' if official else params.get('runmode', 'CYCLES')
        downsample = 10 if official else params.get('downsample', -1)
        if runmode not in {'CYCLES', 'FULL'}:
            raise ValueError('AA runmode must be CYCLES or FULL')
        aa_identity = dict(component=component_identity('aa', runtime, entry, verified, aa_env),
            seeds=[{k: s[k] for k in ('chrom', 'start', 'end', 'seed_id')} for s in seeds],
            inputs=[current_bam, current_index], reference=config['reference'],
            reference_identity=[r for r in config.get('_reference_identity', [])
                                if r.get('kind') != 'local_cnv_background_bed'],
            build=sample['reference_build'], policy=policy,
            background=dict(mode=background['mode'], raw_stats_sha256=background.get('sha256')),
            parameters=dict(policy=parameter_policy, downsample=downsample, runmode=runmode, extendmode='EXPLORE'),
            executor_sha256=digest(__file__),
            policy_executor_sha256=digest(Path(__file__).with_name('concordance_policy.py')),
            stats_handoff_sha256=digest(Path(__file__).with_name('concordance_aa_stats.py')))
        aa_identity['output_validator_sha256'] = digest(Path(__file__).with_name('parsing.py'))
        if official and not fixed_hash:
            # A fresh interpreter may choose another read subset. Resume the
            # same completed Suite job but avoid substitution across jobs.
            aa_identity['unseeded_job_scope'] = str(out.resolve())
        aa_key = cache.signature(aa_identity)
        ac_component = component_identity('ac', runtime, entry, verified, ac_env)
        from .concordance_import import ac_cli_flags
        ac_flags = ac_cli_flags(config, vid)
        complete_signature = cache.signature(dict(aa=aa_key, ac=ac_component,
            ac_flags=ac_flags,
            ac_input_writer_sha256=digest(Path(__file__).with_name('runner.py')),
            parser_sha256=digest(Path(__file__).with_name('parsing.py'))))
        # Preserve the existing attempt path on an unchanged, fully valid resume.
        with measure('downstream_cache', measurement_kind='wall', sample_id=sid, version=vid,
                     cache_status='check') as timing:
            try:
                cached = api._cached_result(out, complete_signature)
            except (OSError, ValueError, TypeError, KeyError, AttributeError):
                cached = None
            timing.mark(cache_status='hit' if cached else 'miss')
        if cached:
            fresh = parse_run(cached['run_dir'], sid, '__selected_seeds__', vid, supported,
                              reference_lengths=config.get('_reference_lengths'))
            if fresh['status'] in ('ok', 'not_detected'):
                component_cache = copy.deepcopy(cached.get('component_cache', {}))
                for info in component_cache.values():
                    info['hit'] = True
                fresh.update(component_cache=component_cache, downstream_cache_hit=True)
                record, attempt = fresh, Path(cached['run_dir'])
                print('  verified AA and AC outputs reused', flush=True)
                return _finish(record, attempt, sample, entry, component_cache)
        attempt = _new_attempt(api, out)
        (attempt / 'aa').mkdir()
        (attempt / 'ac').mkdir()
        bed = attempt / 'selected_seeds.bed'
        # AA artifacts are published before attempting AC.  An AC failure never
        # turns the successful AA stage into a failed cache or forces AA rerun.
        stage = 'aa'
        with cache.stage_lock(cache_root, 'aa', aa_key):
            with measure('aa_cache', measurement_kind='wall', sample_id=sid, version=vid,
                         cache_status='check') as timing:
                aa_cached = cache.read_entry(cache_root, 'aa', aa_key)
                if not aa_cached and config.get('reuse_aa_from'):
                    aa_cached = cache.read_entry(Path(config['reuse_aa_from']) / '.component_cache', 'aa', aa_key)
                if aa_cached:
                    cache.materialize(aa_cached, attempt)
                    timing.status = 'complete'
                timing.mark(cache_status='hit' if aa_cached else 'miss',
                            producer_id=aa_cached['producer_run_dir'] if aa_cached else None)
            if aa_cached:
                component_cache['aa'] = cache.metadata(aa_cached, True)
                stages['aa'] = _staged(aa_cached, component_cache['aa'])
                print('  verified exact-commit AA output reused', flush=True)
            else:
                bed.write_text(''.join(f'{s["chrom"]}\t{s["start"]}\t{s["end"]}\t{s["seed_id"]}\n' for s in seeds))
                env = dict(aa_env, AA_DATA_REPO=str(api.local_reference(
                    config['reference']['aa_data_repo'], attempt / 'reference')),
                    MPLCONFIGDIR=str(attempt / 'mplcache'))
                aaargs = ['--bam', sample['alignment'], '--bed', str(bed),
                    '--out', str(attempt / 'aa' / 'reconstruction'), '--ref', sample['reference_build'],
                    '--runmode', runmode, '--extendmode', 'EXPLORE',
                    '--downsample', str(downsample)]
                help_text = verified.get('aa', {}).get('help', '')
                evidence = {k: v for k, v in background.items() if k != 'fields'}
                evidence.update(original_bam=sample['alignment'], reuse_confirmed=False,
                                python_hash_seed=hash_seed)
                if reuse_stats:
                    aaargs[1] = api.stage_measurement(background, attempt, sample['alignment'],
                                                      allow_alias=not official)
                    evidence.update(reader_bam_path=aaargs[1],
                        staged_sha256=digest(attempt / 'reference' / 'coverage.stats'))
                elif not official and '--no_cstats' in help_text:
                    aaargs.append('--no_cstats')
                solver = aa_identity['component']['solver']
                if '--solver' in help_text:
                    aaargs.extend(['--solver', solver])
                elif solver != 'mosek':
                    raise ValueError('Historical AA solver cannot be silently replaced')
                atomic_json(attempt / 'aa_background.json', evidence)
                command = api.guarded_command(runtime['aa_python'], runtime['aa_script'], aaargs, config['_io_policy_path'])
                stages['aa'] = dict(status='running', argv=command, started=now(),
                                    parameter_policy=parameter_policy,
                                    solver_selection=verified.get('official_solver'),
                                    cache_hit=False, cache_key=aa_key, execution='subprocess')
                atomic_json(attempt / 'stages.json', stages)
                with measure('aa', measurement_kind='subprocess', sample_id=sid, version=vid,
                             log_path=attempt / 'aa.log', cache_status='miss') as timing:
                    rc = api.execute(command, attempt, env, attempt / 'aa.log',
                                     timeout=config.get('job_timeout_seconds', 3600))
                    timing.status = 'complete' if rc == 0 else 'failed'
                stages['aa'].update(status='complete' if rc == 0 else 'failed', returncode=rc, ended=now())
                atomic_json(attempt / 'stages.json', stages)
                aa_log = (attempt / 'aa.log').read_text(errors='replace')
                if reuse_stats:
                    evidence['reuse_confirmed'] = ('Reusing cstats from ' + str(attempt / 'reference' / 'coverage.stats') in aa_log
                        and digest(attempt / 'reference' / 'coverage.stats') == evidence['staged_sha256'])
                    atomic_json(attempt / 'aa_background.json', evidence)
                for audit in sorted(attempt.glob('bam_access.*.json')):
                    access = json.loads(audit.read_text())
                    if access.get('blocked'):
                        raise RuntimeError(access['blocked'][0] + '; see ' + audit.name)
                if rc:
                    raise RuntimeError('Bounded AA failed; see aa.log')
                if reuse_stats and not evidence['reuse_confirmed'] and not official:
                    raise RuntimeError('Native AA did not confirm use of its staged coverage.stats')
                if 'Traceback (most recent call last)' in aa_log:
                    raise RuntimeError('AA logged an exception despite exit zero')
                # Strict validation gates publication, while ordinary result
                # parsing retains unaffected amplicons from a partial native run.
                # A process exit of zero is not by itself a reusable AA success.
                validation_error = None
                try:
                    _validate_aa(attempt, config.get('_reference_lengths'))
                except (OSError, ValueError) as exc:
                    validation_error = str(exc)
                    stages['aa']['output_validation'] = dict(status='unevaluable', reason=str(exc))
                cycles = sorted((attempt / 'aa').glob('*_amplicon*_cycles.txt'))
                if aa_amplicon_count(attempt) != len(cycles):
                    raise ValueError('Incomplete AA output')
                if cycles:
                    api.prepare_ac_input(attempt)  # also checks complete, nonempty graph/cycles pairs
                if validation_error is None:
                    aa_cached = cache.publish(cache_root, 'aa', aa_identity, attempt,
                                              _files(attempt, 'aa'), stages['aa'])
                    component_cache['aa'] = cache.metadata(aa_cached, False)
                else:
                    component_cache['aa'] = dict(hit=False, key=aa_key, published=False,
                        producer_run_dir=str(attempt), validation_issue=validation_error)
            atomic_json(attempt / 'stages.json', stages)
        cycles = sorted((attempt / 'aa').glob('*_amplicon*_cycles.txt'))
        if cycles:
            stage = 'ac'
            inputfile = api.prepare_ac_input(attempt)
            ac_identity = dict(component=ac_component, aa_key=aa_key,
                aa_outputs=_aa_input_files(attempt), build=sample['reference_build'],
                reference=config['reference'],
                reference_identity=[r for r in config.get('_reference_identity', [])
                                    if r.get('kind') != 'local_cnv_background_bed'],
                ac_flags=ac_flags,
                input_writer_sha256=digest(Path(__file__).with_name('runner.py')),
                parser_sha256=digest(Path(__file__).with_name('parsing.py')))
            ac_key = cache.signature(ac_identity)
            with cache.stage_lock(cache_root, 'ac', ac_key):
                with measure('ac_cache', measurement_kind='wall', sample_id=sid, version=vid,
                             cache_status='check') as timing:
                    ac_cached = cache.read_entry(cache_root, 'ac', ac_key)
                    if ac_cached:
                        cache.materialize(ac_cached, attempt)
                        timing.status = 'complete'
                    timing.mark(cache_status='hit' if ac_cached else 'miss',
                                producer_id=ac_cached['producer_run_dir'] if ac_cached else None)
                if ac_cached:
                    component_cache['ac'] = cache.metadata(ac_cached, True)
                    stages['ac'] = _staged(ac_cached, component_cache['ac'])
                else:
                    command = [runtime['ac_python'], runtime['ac_script'], '-i', str(inputfile),
                        '--ref', sample['reference_build'], '-o', str(attempt / 'ac' / 'classification'),
                        *ac_flags]
                    stages['ac'] = dict(status='running', argv=command, started=now(),
                                        parameter_policy=parameter_policy,
                                        cache_hit=False, cache_key=ac_key, execution='subprocess')
                    atomic_json(attempt / 'stages.json', stages)
                    env = dict(ac_env, AA_DATA_REPO=str(api.local_reference(
                        config['reference']['aa_data_repo'], attempt / 'reference')),
                        MPLCONFIGDIR=str(attempt / 'mplcache'))
                    with measure('ac', measurement_kind='subprocess', sample_id=sid, version=vid,
                                 log_path=attempt / 'ac.log', cache_status='miss') as timing:
                        rc = api.execute(command, attempt, env, attempt / 'ac.log',
                                         timeout=config.get('job_timeout_seconds', 3600))
                        timing.status = 'complete' if rc == 0 else 'failed'
                    stages['ac'].update(status='complete' if rc == 0 else 'failed', returncode=rc, ended=now())
                    atomic_json(attempt / 'stages.json', stages)
                    if rc:
                        raise RuntimeError('AC failed; see ac.log')
                    log = (attempt / 'ac.log').read_text(errors='replace')
                    native_log = attempt / 'ac' / 'classification.log'
                    if native_log.is_file():
                        log += '\n' + native_log.read_text(errors='replace')
                    if 'Traceback (most recent call last)' in log:
                        raise RuntimeError('AC traceback despite process exit status')
                    api.check_ac_runtime(log, entry['ac']['version'])
                    if entry['ac']['version'].startswith('2.') and (not native_log.is_file() or not re.search(
                            r' - INFO - Complete\s*$', native_log.read_text(errors='replace'))):
                        raise RuntimeError('AC2 final completion marker missing')
                atomic_json(attempt / 'stages.json', stages)
                stage = 'result_reading'
                record = parse_run(attempt, sid, '__selected_seeds__', vid, supported,
                                   reference_lengths=config.get('_reference_lengths'))
                if not ac_cached and record['status'] in ('ok', 'not_detected'):
                    ac_cached = cache.publish(cache_root, 'ac', ac_identity, attempt,
                                              _files(attempt, 'ac'), stages['ac'])
                    component_cache['ac'] = cache.metadata(ac_cached, False)
        else:
            stages['ac'] = dict(status='skipped_zero_amplicons', execution='not_required')
            atomic_json(attempt / 'stages.json', stages)
            stage = 'result_reading'
            record = parse_run(attempt, sid, '__selected_seeds__', vid, supported,
                               reference_lengths=config.get('_reference_lengths'))
    except KeyboardInterrupt:
        attempt = attempt or _new_attempt(api, out)
        stages.setdefault(stage, {}).update(status='interrupted', reason='User interruption; no result inferred')
        atomic_json(attempt / 'stages.json', stages)
        raise
    except Exception as exc:
        attempt = attempt or _new_attempt(api, out)
        stages.setdefault(stage, {}).update(status='failed', reason=str(exc))
        atomic_json(attempt / 'stages.json', stages)
        record = api.empty_record(sid, '__selected_seeds__', vid, 'failed', stage, str(exc), supported_labels=supported)
    finally:
        if attempt:
            alias = attempt / 'coverage_input.bam'
            if alias.is_symlink():
                alias.unlink()
    record = _finish(record, attempt, sample, entry, component_cache)
    if record['status'] in ('ok', 'not_detected'):
        api._save_cache(out, complete_signature, record)
        # Preserve the old cache.json schema while also protecting the stage's
        # command inputs, root logs and access audits on whole-job resume.
        path = out / 'cache.json'
        complete = json.loads(path.read_text())
        for rel in _files(attempt, 'aa') + _files(attempt, 'ac'):
            complete['files'][str(attempt / rel)] = digest(attempt / rel)
        atomic_json(path, complete)
    return record


def _finish(record, attempt, sample, entry, component_cache):
    record.update(run_dir=str(attempt), analysis_kind='aa_ac', source_version=sample['source_version'],
        source_aa_version=entry['aa']['version'], provenance_kind='bounded_reconstruction',
        comparison_eligible=record['status'] in ('ok', 'not_detected'),
        reference_scope='bounded_indexed_BAM_connected_reconstruction',
        component_versions={k: entry[k]['version'] for k in ('aa', 'ac')},
        component_cache=component_cache)
    return record
