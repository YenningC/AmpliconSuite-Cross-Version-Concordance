"""Coordinate genome-wide seed discovery and paired, version-native reconstruction."""
import json
import os
import shlex
import random
from pathlib import Path

from .common import atomic_json
from .concordance_cnv import load_builtin_map, LocalFailure
from .concordance_timing import timed_call, record_native_aa
from .runner import empty_record


def selected_native(seeds, loci):
    """Select original BED rows, preserving their order, fields and boundaries."""
    spans = {}
    for locus in loci:
        spans.setdefault(locus['chrom'], []).append((locus['start'], locus['end']))
    return [s for s in seeds if any(s['start'] < end and start < s['end']
                                    for start, end in spans.get(s['chrom'], ()))]


def _coordinates(seeds):
    return {(s['chrom'], s['start'], s['end']) for s in seeds}


def _historical_seed_command_matches(sample, config):
    """A user-supplied BED alone cannot establish reproduction of the original run."""
    evidence = sample.get('source_cnvkit_evidence', {})
    if evidence.get('origin') != 'source_result_run_metadata' or not sample.get('seeds_bed'):
        return False
    try:
        metadata = json.loads(Path(evidence['path']).read_text())
        from .concordance_global_cnv import _historical
        spec = sample['source_result']
        archive = _historical(spec.get('result_dir') or spec['run_dir'], sample, config,
                              sample['source_cnvkit_version'], spec)
        wrapper = shlex.split(metadata['AmpliconSuite-pipeline_command'])
        if set(wrapper) & {'--purity', '--ploidy', '--no_filter', '--cnv_bed'}:
            return False
        args = shlex.split(metadata['amplified_intervals_cmd'])
        def arg(name):
            return args[args.index(name) + 1]
        cnv_dir = Path(archive['file']).parent
        base = Path(arg('--bam')).stem
        expected = {cnv_dir / (base + '_CNV_CALLS' + ending)
                    for ending in ('_pre_filtered.bed', '_unfiltered_gains.bed')}
        thresholds = config.get('full_suite', {})
        return (Path(arg('--bed')).resolve() in expected and Path(arg('--bed')).is_file()
                and float(arg('--gain')) == thresholds.get('cngain', 4.5)
                and float(arg('--cnsize_min')) == thresholds.get('cnsize_min', 50000)
                and arg('--ref') == sample['reference_build']
                and os.path.samefile(arg('--bam'), sample['alignment'])
                and os.path.samefile(arg('--out') + '.bed', sample['seeds_bed']))
    except (KeyError, ValueError, IndexError, OSError):
        return False


def _source_evidence(sample, baseline, sets, joints, selected, out, config):
    from .concordance import _bridge
    source = sample['source_version']
    result = sets[source]
    original = _coordinates(sample['seeds'])
    observed = _coordinates(result.get('seeds', []))
    comparable = bool(result['status'] == 'ok' and sample.get('seeds_bed')
                      and sample.get('source_cnvkit_version') == result.get('cnvkit_version')
                      and _historical_seed_command_matches(sample, config))
    calibration = ('passed' if original == observed else 'failed') if comparable else 'unverified'
    sample['global_source_seed_calibration'] = calibration
    bridge = []
    if source in joints:
        bridge = _bridge(baseline, joints[source], selected)
    sample['source_bridge_status'] = 'unverified'
    if calibration == 'failed':
        sample['source_bridge_status'] = 'failed'
    elif (calibration == 'passed' and bridge and baseline.get('comparison_eligible')
          and joints.get(source, {}).get('comparison_eligible')):
        evaluable = [r for r in bridge if r.get('evaluable')]
        if any(not r['passed'] for r in evaluable):
            sample['source_bridge_status'] = 'failed'
        elif len(evaluable) == len(bridge):
            sample['source_bridge_status'] = 'passed'
    evidence = dict(source_suite=source, historical_cnvkit=sample.get('source_cnvkit_version'),
                    panel_cnvkit=result.get('cnvkit_version'), comparable=comparable,
                    original_seed_bed=sample.get('seeds_bed'), exact_coordinate_agreement=calibration,
                    source_bridge_status=sample['source_bridge_status'],
                    observed_seed_coordinates_match=original == observed if sample.get('seeds_bed') and result['status'] == 'ok' else None,
                    original_seed_count=len(original) if sample.get('seeds_bed') else None,
                    replay_seed_count=len(observed) if result['status'] == 'ok' else None,
                    historical_only=sorted(original - observed) if comparable else None,
                    replay_only=sorted(observed - original) if comparable else None)
    atomic_json(out / 'source_bridge' / (sample['sample_id'] + '.json'), bridge)
    atomic_json(out / 'source_bridge' / (sample['sample_id'] + '.seeds.json'), evidence)


def _run_native_aa(sample, native, call, runtime, entry, verified, settings, seed_settings, out):
    from . import concordance as api
    from .concordance_global_seeds import global_seed_measurement

    sid, vid = sample['sample_id'], runtime['id']
    try:
        measurement = global_seed_measurement(Path(call['seed_attempt']), call['seeds'],
                                               sample, runtime, verified, seed_settings)
        joint = timed_call('selected_seed_downstream', api.run_bounded_downstream, sample, native,
                           runtime, entry, verified, settings, out / 'downstream_jobs' / sid / vid,
                           measurement, timing=dict(sample_id=sid, version=vid))
        if joint.get('run_dir'):
            cache = joint.get('component_cache', {}).get('aa', {})
            record_native_aa(joint['run_dir'], sid, vid, role='producer_reference' if cache.get('hit')
                             else 'current_execution', producer_id=cache.get('producer_run_dir', ''))
        return joint
    except (OSError, ValueError, KeyError, RuntimeError, LocalFailure) as exc:
        return empty_record(sid, '__selected_seeds__', vid, 'failed', getattr(exc, 'stage', 'aa'), str(exc),
                            analysis_kind='aa_ac', source_version=sample['source_version'],
                            comparison_eligible=False, detected=None)


def _recover_global_sample(sample, baseline, versions, entries, config, out, verified, verified_runtime, recovery):
    from .concordance_global_report import global_report_units

    sid = sample['sample_id']
    prior, sampling = recovery['previous_sample'], recovery['sampling']
    sets = sample['global_seed_sets'] = prior['global_seed_sets']
    for name in ('upstream_seed_ids', 'selected_seed_ids', 'global_executed_seed_ids'):
        sample[name] = prior[name]
    selected = sampling['selected']
    settings = dict(config, _component_cache_dir=str(out / '.component_cache'))
    seed_settings = dict(settings, _io_policy_path=prior['_policy_path'])
    prior_policy = json.loads(Path(prior['_policy_path']).read_text())
    rescue_limit = config['aa_rescue_unique_bp']
    rescue_policy = dict(prior_policy, max_unique_bp=rescue_limit)
    joints = dict(recovery['completed'])
    final_states = dict(prior['global_version_statuses'])
    retries = recovery['retries']
    audit = dict(source=str(Path(config['recover_aa_from']).resolve()),
                 prior_report_sha256=recovery['checksums'], prior_seed_sha256=recovery['seeds'],
                 completed_versions_reused=list(recovery['completed']), unique_rescue_limit_bp=rescue_limit,
                 cnvkit_and_seed_filter_executed=False, retries=retries)
    audit_path = out / 'aa_recovery' / (sid + '.json')
    atomic_json(out / 'sampling' / (sid + '.json'), sampling)
    atomic_json(out / 'global_upstream' / sid / 'discovery.json', sets)
    atomic_json(audit_path, audit)
    print(f'  {sid}: verified prior CNVkit and seed files; AA recovery will retry {len(retries)} Suite paths', flush=True)
    for vid, evidence in retries.items():
        call = sets[vid]
        native = selected_native(call['seeds'], selected)
        if [s['seed_id'] for s in native] != sample['global_executed_seed_ids'][vid]:
            raise ValueError(f'AA recovery changed the original AA seed rows for {sid} / {vid}')
        runtime = verified_runtime(vid, sid)
        original_solver = evidence['solver']
        policy_path = out / 'policies' / f'{sid}.{vid}.aa_recovery.json'
        atomic_json(policy_path, rescue_policy)
        evidence.update(attempt_policy=str(policy_path), effective_requested_limit_bp=rescue_limit)
        print(f'  {sid} / {vid}: retrying AA with max_unique_bp={rescue_limit:,}', flush=True)
        if original_solver and verified[vid].get('official_solver', {}).get('requested') != original_solver['requested']:
            joint = empty_record(sid, '__selected_seeds__', vid, 'failed', 'runtime_check',
                                 'AA recovery could not verify the original solver',
                                 analysis_kind='aa_ac', source_version=sample['source_version'],
                                 comparison_eligible=False, detected=None)
        else:
            joint = _run_native_aa(sample, native, call, runtime, entries[vid], verified[vid],
                                   dict(settings, _io_policy_path=str(policy_path)), seed_settings, out)
        joints[vid] = joint
        stage = joint.get('stage', 'result_reading')
        final_states[vid] = dict(status=joint['status'], stage=stage)
        evidence.update(status=joint['status'], stage=stage, reason=joint.get('reason'),
                        run_dir=joint.get('run_dir'))
        atomic_json(audit_path, audit)
    sample['global_version_statuses'] = final_states
    _source_evidence(sample, baseline, sets, joints, selected, out, config)
    upstream, newly_projected = global_report_units(sample, joints, versions)
    new_rows = {(r['version'], r['seed_id']): r for r in newly_projected}
    downstream = [new_rows[vid, row['seed_id']] if vid in retries else row
                  for vid in versions for row in recovery['previous_downstream'][vid]]
    return upstream, downstream


def run_global_sample(sample, baseline, versions, entries, config, out, runtimes, verified, verified_runtime,
                      *, aa_recovery=None):
    if config.get('recover_aa_from'):
        if aa_recovery is None:
            from .concordance_aa_recovery import verify_recovery_source
            aa_recovery = verify_recovery_source(config, out)
        return _recover_global_sample(sample, baseline, versions, entries, config, out, verified,
                                      verified_runtime, aa_recovery[sample['sample_id']])
    from . import concordance as api
    from .concordance_global_cnv import global_cns
    from .concordance_global_report import global_loci, global_report_units
    from .concordance_global_seeds import global_suite_seeds

    sid = sample['sample_id']
    settings = dict(config, _io_policy_path=sample['_policy_path'],
                    _component_cache_dir=str(out / '.component_cache'))
    sets = sample['global_seed_sets'] = {}
    jobs = out / 'global_upstream' / sid
    mapped = load_builtin_map()['versions']
    for vid in versions:
        entry = entries[vid]
        expected = mapped[vid]['cnvkit_version']
        record = sets[vid] = dict(status='missing', stage='runtime_check', reason='No configured Suite runtime',
                                  seeds=[], cnvkit_version=expected, native_recipe_match=False,
                                  reference_mode='matched_normal' if sample.get('normal_alignment') and sample['reference_build'] != 'GRCh38_viral' else 'suite_reference',
                                  supported_labels=entry.get('supported_labels', entry.get('supported_classes', [])))
        if vid not in runtimes:
            continue
        runtime = runtimes[vid]
        timing = dict(sample_id=sid, version=vid)
        stage = 'cnvkit'
        try:
            print(f'  {sid} / {vid}: global CNVkit {expected} and native seed selection', flush=True)
            cnv = timed_call('global_cnvkit', global_cns, sample, runtime, entry, settings,
                             out / 'global_cnv_cache', timing=timing)
            if cnv['cnvkit_version'] != expected:
                raise ValueError('Global output CNVkit version differs from the pinned Suite recipe')
            record.update(cns=str(cnv['cns']), cnvkit_provenance=str(cnv['provenance']), reuse=cnv['reuse'],
                          evidence_level=cnv['evidence_level'],
                          historical_provenance_limited=cnv['evidence_level'] != 'native_verified')
            stage = 'seed_processing'
            attempt = api._new_attempt(jobs / vid)
            record['seed_attempt'] = str(attempt)
            runtime = verified_runtime(vid, sid)
            seed_config = dict(settings, _seed_verified_runtime=verified[vid])
            seeds = timed_call('global_suite_seed_filter', global_suite_seeds, sample, runtime, entry,
                               seed_config, cnv['cns'], attempt, timing=timing)
            record.update(status='ok', stage=stage, seeds=seeds, native_recipe_match=True,
                          reason='Global segmentation and version-specific native seed filter completed')
        except (OSError, ValueError, KeyError, RuntimeError, LocalFailure) as exc:
            record.update(status=getattr(exc, 'status', 'failed'), stage=getattr(exc, 'stage', stage),
                          reason=str(exc))
        atomic_json(jobs / 'discovery.json', sets)

    # Freeze the common comparison units before looking at any AA/AC output.
    loci = global_loci(sets)
    sampling = config.get('sampling', {})
    cap = sampling.get('max_seeds', 10)
    rng = random.Random(str(sampling.get('random_seed', 20261117)) + ':' + sid)
    selected = loci if cap == 'all' or len(loci) <= cap else rng.sample(loci, cap)
    selected = sorted(selected, key=lambda s: (s['chrom'], s['start'], s['end']))
    sample['upstream_seed_ids'] = [s['seed_id'] for s in loci]
    sample['selected_seed_ids'] = [s['seed_id'] for s in selected]
    executed = sample['global_executed_seed_ids'] = {}
    joints = {}
    for vid in versions:
        result = sets[vid]
        if result['status'] != 'ok':
            continue
        native = selected_native(result['seeds'], selected)
        executed[vid] = [s['seed_id'] for s in native]
        if not native:
            continue
        print(f'  {sid} / {vid}: AA/AC from {len(native)} native seed intervals', flush=True)
        runtime = verified_runtime(vid, sid)
        entry = entries[vid]
        joint = _run_native_aa(sample, native, result, runtime, entry, verified[vid], settings, settings, out)
        joints[vid] = joint
    atomic_json(out / 'sampling' / (sid + '.json'), dict(all_seeds=loci, selected=selected, parameters=sampling,
                native_selected_ids=executed, rule='paired_sampling_from_union_of_global_overlapping_seed_loci'))
    atomic_json(jobs / 'discovery.json', sets)
    _source_evidence(sample, baseline, sets, joints, selected, out, config)
    final_states = {v:dict(status=r['status'], stage=r['stage']) for v,r in sets.items()}
    for vid, joint in joints.items():
        final_states[vid] = dict(status=joint['status'], stage=joint.get('stage', 'result_reading'))
    sample['global_version_statuses'] = final_states
    return global_report_units(sample, joints, versions)
