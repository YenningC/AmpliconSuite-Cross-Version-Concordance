"""AmpliconSuite Cross-Version Concordance: result-first, bounded, local version comparison."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import random
import re
import sys
from . import __version__
from .common import atomic_json, digest, now, table
from .inputs import safe_id
from .parsing import merge_intervals
from .runner import empty_record, execute, local_reference, check_ac_runtime
from .versions import load_lock, verify_runtime
from .benchmark_policy import guarded_command

PROGRAM_VERSION = __version__
FIXED_REFERENCE = 'suite-1.6.2'
PATH_KEYS = {'alignment', 'index', 'normal_alignment', 'normal_index', 'seeds_bed',
             'source_seed_bed', 'baseline_cns', 'aa_dir', 'ac_dir', 'run_dir',
             'result_dir', 'stages_json', 'baseline_provenance', 'fasta', 'aa_data_repo', 'background_bed',
             'aa_script', 'ac_script', 'suite_script', 'cnvkit_script', 'version_lock'}


def fast_fingerprint(path):
    """Three bounded 64 KiB reads, never a full BAM digest."""
    p = Path(path).resolve(strict=True)
    stat = p.stat()
    h = hashlib.sha256()
    offsets = sorted(set((0, max(0, stat.st_size // 2 - 32768), max(0, stat.st_size - 65536))))
    with p.open('rb') as handle:
        for offset in offsets:
            handle.seek(offset)
            h.update(str(offset).encode())
            h.update(handle.read(65536))
    return dict(path=str(p), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                sampled_sha256=h.hexdigest(), method='three_64KiB_blocks_not_full_integrity')


def _resolve(value, parent):
    if isinstance(value, list):
        return [_resolve(v, parent) for v in value]
    if not isinstance(value, dict):
        return value
    result = {}
    for k, v in value.items():
        if v and (k in PATH_KEYS or k.endswith('_python')):
            p = Path(v).expanduser()
            p = parent / p if not p.is_absolute() else p
            result[k] = str(p.absolute() if k.endswith('_python') else p.resolve())
        elif k == 'extra_path':
            result[k] = [str((parent / Path(p).expanduser()).resolve()) for p in v]
        else:
            result[k] = _resolve(v, parent)
    return result


def discover_result(spec):
    """Only inspect known Suite layout levels; no recursive patient-tree crawl."""
    spec = copy.deepcopy(spec)
    if spec.get('run_dir'):
        root = Path(spec['run_dir'])
        spec.setdefault('aa_dir', str(root / 'aa'))
        spec.setdefault('ac_dir', str(root / 'ac'))
        if (root / 'stages.json').is_file():
            spec.setdefault('stages_json', str(root / 'stages.json'))
    if spec.get('result_dir'):
        root = Path(spec['result_dir'])
        for key, pattern in (('aa_dir', '*_AA_results'), ('ac_dir', '*_classification')):
            if not spec.get(key):
                found = [p for p in root.glob(pattern) if p.is_dir()]
                direct = root / ('aa' if key == 'aa_dir' else 'ac')
                if direct.is_dir():
                    found.append(direct)
                if len(found) == 1:
                    spec[key] = str(found[0])
        if not spec.get('seeds_bed'):
            found = list(root.glob('*_AA_CNV_SEEDS.bed'))
            if len(found) == 1:
                spec['seeds_bed'] = str(found[0])
    return spec


def read_seed_bed(path, lengths=None):
    rows, seen = [], set()
    with Path(path).open() as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f'Invalid seed BED line {n}')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if not 0 <= start < end or lengths is not None and (chrom not in lengths or end > lengths[chrom]):
                raise ValueError(f'Seed outside declared reference at line {n}')
            key = (chrom, start, end)
            if key in seen:
                raise ValueError('Duplicate seed coordinates')
            seen.add(key)
            ident = 'seed_' + hashlib.sha256(f'{chrom}:{start}-{end}'.encode()).hexdigest()[:12]
            rows.append(dict(chrom=chrom, start=start, end=end, seed_id=ident,
                             source_fields=fields[3:], source_line=n))
    return sorted(rows, key=lambda r: (r['chrom'], r['start'], r['end']))


def _signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _reference_lengths(config):
    fasta = config.get('reference', {}).get('fasta')
    if not fasta:
        return None
    fai = Path(fasta + '.fai')
    if not fai.is_file():
        raise ValueError('Existing FASTA index required; automatic indexing is disabled')
    return {f[0]: int(f[1]) for line in fai.read_text().splitlines() if (f := line.split())}


def load_manifest(path):
    path = Path(path).resolve(strict=True)
    config = _resolve(json.loads(path.read_text()), path.parent)
    if config.get('mode') not in ('with_bam', 'results_only'):
        raise ValueError('mode must be with_bam or results_only')
    if config.get('reference_version', FIXED_REFERENCE) != FIXED_REFERENCE:
        raise ValueError('This release fixes its reference at ' + FIXED_REFERENCE)
    config['reference_version'] = FIXED_REFERENCE
    if config.get('analysis_scope') == 'full_suite' or config.get('alignment_fingerprint') == 'full':
        raise ValueError('Whole-genome operations and full BAM hashing are prohibited')
    lock = load_lock(config.get('version_lock'))
    if lock['reference_version'] != FIXED_REFERENCE:
        raise ValueError('Version lock reference differs from this program release')
    entries = {e['id']: e for e in lock['versions'] if e.get('suite')}
    targets = config.get('target_versions', list(entries))
    if not isinstance(targets, list) or not targets or len(targets) != len(set(targets)):
        raise ValueError('target_versions must be a nonempty unique list')
    if any(v not in entries for v in targets):
        raise ValueError('Every target must have a built-in verified source lock')
    if FIXED_REFERENCE not in targets:
        targets = [FIXED_REFERENCE, *targets]
    config['target_versions'] = targets
    samples = config.get('samples')
    if not isinstance(samples, list) or not samples:
        raise ValueError('At least one sample with existing Suite output is required')
    ids = [safe_id(s['sample_id']) for s in samples]
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate sample_id')
    if any('patient_id' in sample for sample in samples):
        for sample in samples:
            patient = sample.get('patient_id')
            if not isinstance(patient, str) or not patient:
                raise ValueError('When patient_id is used, every sample must have a nonempty patient_id')
            safe_id(patient)
    for s in samples:
        if not s.get('source_version'):
            raise ValueError('source_version is required in both input modes')
        safe_id(s['source_version'])
        if s['source_version'] not in entries:
            raise ValueError('Source Suite is not in the source lock; add an audited lock entry, never guess components')
        if not s.get('reference_build'):
            raise ValueError('Every sample must declare its reference_build')
        s['source_result'] = discover_result(s.get('source_result', {}))
        if s['source_result'].get('version', s['source_version']) != s['source_version']:
            raise ValueError('source_result.version contradicts source_version')
        if not (s['source_result'].get('aa_dir') and s['source_result'].get('ac_dir')):
            raise ValueError('source_result must identify AA and AC directories')
        s['comparison_results'] = [discover_result(r) for r in s.get('comparison_results', [])]
        versions = [s['source_version']] + [r['version'] for r in s['comparison_results']]
        if len(versions) != len(set(versions)):
            raise ValueError('Duplicate imported result version for a sample')
        if any(v not in entries for v in versions):
            raise ValueError('Imported comparison version absent from lock')
        if config['mode'] == 'with_bam':
            if not s.get('alignment') or not s.get('index'):
                raise ValueError('with_bam requires alignment and existing index')
            if not s['alignment'].lower().endswith('.bam'):
                raise ValueError('This bounded release accepts indexed BAM only; CRAM conversion is prohibited')
            seed = s.get('seeds_bed') or s['source_result'].get('seeds_bed')
            if not seed:
                raise ValueError('Existing final Suite seed BED required; no automatic genome-wide discovery')
            s['seeds_bed'] = seed
    if config['mode'] == 'with_bam' and not all(config.get('reference', {}).get(k) for k in ('build', 'fasta', 'aa_data_repo')):
        raise ValueError('with_bam requires a fixed reference build, FASTA and AA resource directory')
    for runtime in config.get('runtimes', []):
        if runtime['id'] not in entries:
            raise ValueError('Runtime absent from source lock')
        if runtime.get('threads', 1) != 1:
            raise ValueError('Bounded BAM callers use one process so budgets cannot be bypassed by workers')
    runtime_ids = [r['id'] for r in config.get('runtimes', [])]
    if len(runtime_ids) != len(set(runtime_ids)):
        raise ValueError('Duplicate runtime IDs')
    if config.get('threads', 1) != 1:
        raise ValueError('The bounded benchmark currently executes serially')
    for key, allowed in (('aa_parameters', {'downsample'}), ('ac_parameters', {'min_size', 'min_flow'})):
        if not isinstance(config.get(key, {}), dict) or set(config.get(key, {})) - allowed:
            raise ValueError('Unsupported parameters in ' + key + '; no parameters are silently ignored')
    import math
    for key, default in (('min_size', 50000), ('min_flow', 1.0)):
        value = float(config.get('ac_parameters', {}).get(key, default))
        if not math.isfinite(value) or value <= 0:
            raise ValueError('AC parameters must be positive finite numbers')
    ds = float(config.get('aa_parameters', {}).get('downsample', -1))
    if not math.isfinite(ds) or ds != -1 and ds <= 0:
        raise ValueError('AA downsample must be -1 or positive finite coverage')
    from .benchmark_report import _configuration
    _configuration(config)  # Reject invalid analysis plans before launching component jobs.
    return config, lock, entries


def _validate_bam(sample, config, lengths):
    import pysam
    identities = []
    for key, index_key in (('alignment', 'index'), ('normal_alignment', 'normal_index')):
        if not sample.get(key):
            continue
        p, index = sample[key], sample.get(index_key)
        if not index or not Path(index).is_file() or not p.lower().endswith('.bam'):
            raise ValueError('Existing indexed BAM required for ' + key)
        # Header and index queries only. Never flagstat, sort, quickcheck or index.
        with pysam.AlignmentFile(p, 'rb', index_filename=index) as bam:
            if not bam.has_index():
                raise ValueError('Unreadable BAM index')
            if bam.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
                raise ValueError('BAM header does not declare coordinate sorting')
            found = dict(zip(bam.references, bam.lengths))
            for chrom, size in lengths.items():
                if chrom in found and found[chrom] != size:
                    raise ValueError('BAM/reference contig length mismatch')
            for seed in sample['seeds']:
                if found.get(seed['chrom']) != lengths.get(seed['chrom']):
                    raise ValueError('Seed chromosome absent or inconsistent in BAM')
        identities.append(fast_fingerprint(p))
        identities.append(dict(path=index, sha256=digest(index), method='index_file_full_SHA256'))
    return identities


def _overlap(a, b):
    return a[0] == b[0] and min(a[2], b[2]) > max(a[1], b[1])


def project_record(record, seed):
    """Select complete reconstructed amplicons anchored at a seed; never clip structures."""
    r = copy.deepcopy(record)
    r['seed_id'], r['seed'] = seed['seed_id'], seed
    if r['status'] not in ('ok', 'not_detected'):
        return r
    interval = [seed['chrom'], seed['start'], seed['end']]
    amps = [a for a in r.get('amplicons', []) if any(_overlap(interval, i) for i in a.get('reconstruction_intervals', []))]
    ids = {a['amplicon_id'] for a in amps}
    r['amplicons'] = amps
    r['features'] = [f for f in r.get('features', []) if f['amplicon_id'] in ids]
    r['intervals'] = merge_intervals([i for f in r['features'] for i in f['intervals']])
    r['reconstruction_intervals'] = merge_intervals([i for a in amps for i in a['reconstruction_intervals']])
    r['labels'] = sorted({label for a in amps for label in a['labels']})
    r['status'] = 'ok' if r['labels'] else 'not_detected'
    r['detected'] = bool(r['labels'])
    r['projection'] = 'overlapping_AA_reconstruction_complete_amplicons_retained'
    return r


def _source_spec(spec, entry):
    spec = copy.deepcopy(spec)
    supplied = spec.get('component_versions', {})
    inferred = not bool(supplied or spec.get('aa_version') or spec.get('ac_version'))
    declared, matching = {}, True
    for name in ('aa', 'ac'):
        nested = supplied.get(name, supplied.get(name.upper()))
        if isinstance(nested, dict):
            nested = nested.get('version')
        direct = spec.get(name + '_version')
        if direct and nested and str(direct) != str(nested):
            matching = False
        declared[name] = str(direct or nested or entry[name]['version'])
        matching = matching and declared[name] == entry[name]['version']
    spec['component_versions'] = declared
    spec['validated_declared_components'] = matching
    spec['component_identity_level'] = 'inferred_from_declared_suite' if inferred else 'user_declared_matching_lock'
    return spec


def _cached_result(job, signature):
    path = Path(job) / 'cache.json'
    if not path.is_file():
        return None
    old = json.loads(path.read_text())
    if old.get('signature') != signature:
        return None
    files = old.get('files', {})
    if not files or not all(Path(p).is_file() and digest(p) == sha for p, sha in files.items()):
        return None
    if old['record'].get('status') not in ('ok', 'not_detected'):
        return None
    return old['record']


def _save_cache(job, signature, record):
    files = {}
    run = Path(record.get('run_dir', job))
    for directory in (run / 'aa', run / 'ac'):
        if directory.is_dir():
            for p in directory.rglob('*'):
                if p.is_file() and not p.is_symlink():
                    files[str(p)] = digest(p)
    if (run / 'stages.json').is_file():
        files[str(run / 'stages.json')] = digest(run / 'stages.json')
    atomic_json(Path(job) / 'cache.json', dict(signature=signature, record=record, files=files))


def _new_attempt(job, create=True):
    job = Path(job)
    job.mkdir(parents=True, exist_ok=True)
    existing = [int(p.name[8:]) for p in job.glob('attempt_*') if p.name[8:].isdigit()]
    p = job / f'attempt_{max(existing, default=0)+1:04d}'
    if create:
        p.mkdir()
    return p


def run_bounded_downstream(sample, seeds, runtime, entry, verified, config, out):
    """One joint selected-seed job per version. Connected expansion is not clipped."""
    from .parsing import parse_run, aa_amplicon_count
    supported = entry.get('supported_labels', entry.get('supported_classes', []))
    vid, sid = runtime['id'], sample['sample_id']
    identity = {name: {k: verified.get(name, {}).get(k) for k in ('version', 'manifest_sha256', 'sha256')}
                for name in ('aa', 'ac')}
    identity['packages'] = {name: verified.get(name, {}).get('packages', {}).get('stdout') for name in ('aa', 'ac')}
    identity['bfb'] = verified.get('bfb_probe', {}).get('stdout')
    identity['runtime_variants'] = verified.get('runtime_variants', [])
    effective_env = dict(os.environ, **runtime.get('env', {}))
    identity['solver_environment'] = {k: effective_env.get(k) for k in
        ('MOSEKLM_LICENSE_FILE', 'GRB_LICENSE_FILE', 'DYLD_LIBRARY_PATH', 'LD_LIBRARY_PATH', 'PATH')}
    identity['python_hash_seed'] = effective_env.get('PYTHONHASHSEED')
    license_paths = [Path.home() / 'mosek' / 'mosek.lic']
    for k in ('MOSEKLM_LICENSE_FILE', 'GRB_LICENSE_FILE'):
        value = effective_env.get(k)
        if value and '@' not in value:
            license_paths.append(Path(value).expanduser())
    identity['license_file_metadata_only'] = [{
        'path': str(p), 'exists': p.is_file(),
        'bytes': p.stat().st_size if p.is_file() else None,
        'mtime_ns': p.stat().st_mtime_ns if p.is_file() else None} for p in license_paths]
    signature = _signature(dict(seeds=seeds, inputs=sample['_input_identities'], runtime=runtime,
        verified=identity, reference=config['reference'], reference_identity=config.get('_reference_identity'), policy=config.get('io_policy', {}),
        runner_sha256=digest(__file__), policy_executor_sha256=digest(Path(__file__).with_name('benchmark_policy.py')),
        aa_parameters=config.get('aa_parameters', {}), ac_parameters=config.get('ac_parameters', {})))
    cached = _cached_result(out, signature)
    if cached:
        print('  verified AA+AC cache reused', flush=True)
        fresh = parse_run(cached['run_dir'], sid, '__selected_seeds__', vid, supported,
                          reference_lengths=config.get('_reference_lengths'))
        fresh.update(run_dir=cached['run_dir'], analysis_kind='aa_ac', source_version=sample['source_version'],
            source_aa_version=entry['aa']['version'], provenance_kind='bounded_reconstruction',
            reference_scope='bounded_indexed_BAM_connected_reconstruction',
            component_versions={k: entry[k]['version'] for k in ('aa', 'ac')})
        fresh['comparison_eligible'] = fresh['status'] in ('ok', 'not_detected')
        return fresh
    attempt = _new_attempt(out)
    for name in ('aa', 'ac'):
        (attempt / name).mkdir()
    stages, stage = {}, 'input_check'
    record = None
    try:
        if not verified.get('ready'):
            raise ValueError('Runtime source/version checks failed: ' + '; '.join(verified.get('diagnostics', [])))
        bed = attempt / 'selected_seeds.bed'
        bed.write_text(''.join(f'{s["chrom"]}\t{s["start"]}\t{s["end"]}\t{s["seed_id"]}\n' for s in seeds))
        env = dict(os.environ)
        env.update(runtime.get('env', {}))
        env.update(AA_DATA_REPO=str(local_reference(config['reference']['aa_data_repo'], attempt / 'reference')),
                   PYTHONDONTWRITEBYTECODE='1', MPLCONFIGDIR=str(attempt / 'mplcache'),
                   OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
        bins = [str(Path(runtime[k]).absolute().parent) for k in ('aa_python', 'ac_python')]
        env['PATH'] = os.pathsep.join(bins + runtime.get('extra_path', []) + [env.get('PATH', '')])
        params = config.get('aa_parameters', {})
        aaargs = ['--bam', sample['alignment'], '--bed', str(bed), '--out', str(attempt / 'aa' / 'reconstruction'),
                  '--ref', sample['reference_build'], '--runmode', 'FULL', '--extendmode', 'EXPLORE',
                  '--downsample', str(params.get('downsample', -1))]
        help_text = verified.get('aa', {}).get('help', '')
        if '--no_cstats' in help_text:
            aaargs.append('--no_cstats')
        solver = runtime.get('solver', 'clarabel' if entry['aa']['version'].startswith('1.6.') else 'mosek')
        if '--solver' in help_text:
            aaargs.extend(['--solver', solver])
        elif solver != 'mosek':
            raise ValueError('Historical AA solver cannot be silently replaced')
        stage = 'aa'
        command = guarded_command(runtime['aa_python'], runtime['aa_script'], aaargs, config['_io_policy_path'])
        stages[stage] = dict(status='running', argv=command, started=now())
        atomic_json(attempt / 'stages.json', stages)
        rc = execute(command, attempt, env, attempt / 'aa.log', timeout=config.get('job_timeout_seconds', 3600))
        stages[stage].update(status='complete' if rc == 0 else 'failed', returncode=rc, ended=now())
        atomic_json(attempt / 'stages.json', stages)
        if rc:
            raise RuntimeError('Bounded AA failed or exceeded its I/O budget; see aa.log')
        if 'Traceback (most recent call last)' in (attempt / 'aa.log').read_text(errors='replace'):
            raise RuntimeError('AA logged an exception despite exit zero')
        cycles = sorted((attempt / 'aa').glob('*_amplicon*_cycles.txt'))
        if aa_amplicon_count(attempt) != len(cycles):
            raise ValueError('Incomplete AA output')
        if cycles:
            inputfile = attempt / 'ac.input'
            inputfile.write_text(''.join('reconstruction\t' + str(c) + '\t' + str(c.with_name(c.name.replace('_cycles.txt', '_graph.txt'))) + '\n' for c in cycles))
            p = config.get('ac_parameters', {})
            command = [runtime['ac_python'], runtime['ac_script'], '-i', str(inputfile),
                       '--ref', sample['reference_build'], '-o', str(attempt / 'ac' / 'classification'),
                       '--min_size', str(p.get('min_size', 50000)), '--min_flow', str(p.get('min_flow', 1.0))]
            stage = 'ac'
            stages[stage] = dict(status='running', argv=command, started=now())
            atomic_json(attempt / 'stages.json', stages)
            acenv = dict(env)
            acenv['PATH'] = os.pathsep.join([bins[1], bins[0]] + runtime.get('extra_path', []) + [os.environ.get('PATH', '')])
            rc = execute(command, attempt, acenv, attempt / 'ac.log', timeout=config.get('job_timeout_seconds', 3600))
            stages[stage].update(status='complete' if rc == 0 else 'failed', returncode=rc, ended=now())
            atomic_json(attempt / 'stages.json', stages)
            if rc:
                raise RuntimeError('AC failed; see ac.log')
            log = (attempt / 'ac.log').read_text(errors='replace')
            native_log = attempt / 'ac' / 'classification.log'
            if native_log.is_file():
                log += '\n' + native_log.read_text(errors='replace')
            if 'Traceback (most recent call last)' in log:
                raise RuntimeError('AC traceback despite process exit status')
            check_ac_runtime(log, entry['ac']['version'])
            if entry['ac']['version'].startswith('2.') and (not native_log.is_file() or not re.search(r' - INFO - Complete\s*$', native_log.read_text(errors='replace'))):
                raise RuntimeError('AC2 final completion marker missing')
        else:
            stages['ac'] = dict(status='skipped_zero_amplicons')
            atomic_json(attempt / 'stages.json', stages)
        stage = 'result_reading'
        record = parse_run(attempt, sid, '__selected_seeds__', vid, supported,
                           reference_lengths=config.get('_reference_lengths'))
    except KeyboardInterrupt:
        stages[stage] = dict(status='interrupted', reason='User interruption; no result inferred')
        atomic_json(attempt / 'stages.json', stages)
        raise
    except Exception as exc:
        stages[stage] = dict(status='failed', reason=str(exc))
        atomic_json(attempt / 'stages.json', stages)
        record = empty_record(sid, '__selected_seeds__', vid, 'failed', stage, str(exc), supported_labels=supported)
    record.update(run_dir=str(attempt), analysis_kind='aa_ac', source_version=sample['source_version'],
                  source_aa_version=entry['aa']['version'], provenance_kind='bounded_reconstruction',
                  comparison_eligible=record['status'] in ('ok', 'not_detected'),
                  reference_scope='bounded_indexed_BAM_connected_reconstruction',
                  component_versions={k: entry[k]['version'] for k in ('aa', 'ac')})
    _save_cache(out, signature, record)
    return record


def _bridge(original, replay, seeds):
    from .statistics import _pair, AGREEMENT_METRICS
    rows = []
    for seed in seeds:
        a, b = project_record(replay, seed), project_record(original, seed)
        row = _pair(a, b, FIXED_REFERENCE, 'source_reconstruction_bridge')
        row['passed'] = bool(a.get('comparison_eligible') and b.get('comparison_eligible')
            and row['evaluable'] and all(row.get(m) == 0 for m in AGREEMENT_METRICS))
        if row['evaluable'] and not row.get('signal_positive'):
            row['passed'] = bool(a.get('comparison_eligible') and b.get('comparison_eligible'))
        row['scope'] = 'source-Suite reduced-seed replay agrees on reported metrics only; not structural truth'
        rows.append(row)
    return rows


def run_benchmark(config, lock, entries, out, check_only=False):
    from .benchmark_import import import_result, reclassify_result
    from .benchmark_cnv import run_local_upstream
    from .benchmark_report import write_benchmark_reports
    out = Path(out).resolve()
    # Reject overlap before any write, including beneath an original result tree.
    originals = []
    def collect(value):
        if isinstance(value, dict):
            for k, v in value.items():
                if k in PATH_KEYS and isinstance(v, str):
                    originals.append(Path(v).resolve())
                else:
                    collect(v)
        elif isinstance(value, list):
            for v in value:
                collect(v)
    collect(config)
    for p in originals:
        if p == out or out in p.parents or p.is_dir() and p in out.parents:
            raise ValueError('Output directory overlaps an input/resource directory')
    out.mkdir(parents=True, exist_ok=True)
    import fcntl
    handle = (out / '.benchmark.lock').open('a+')
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return _run_locked(config, lock, entries, out, check_only)
    finally:
        handle.close()


def _run_locked(config, lock, entries, out, check_only):
    from .benchmark_import import import_result, reclassify_cached, _artifacts, MAX_FILE_BYTES, MAX_TOTAL_BYTES
    from .benchmark_cnv import run_local_upstream
    from .benchmark_report import write_benchmark_reports
    config = copy.deepcopy(config)
    if config.get('runtimes') and re.search(r'\s', str(out)):
        raise ValueError('Native AC input tables do not support whitespace in working/output paths; choose a separate path without whitespace')
    lengths = _reference_lengths(config)
    config['_reference_lengths'] = lengths
    run_identity = _signature(dict(config=config, lock=lock, program_version=PROGRAM_VERSION))
    provenance = out / 'provenance.json'
    if provenance.exists() and json.loads(provenance.read_text())['manifest_hash'] != run_identity:
        raise ValueError('Manifest changed; choose a new output directory')
    atomic_json(provenance, dict(program='AmpliconSuite Cross-Version Concordance', program_version=PROGRAM_VERSION,
        created_utc=now(), manifest_hash=run_identity, reference_version=FIXED_REFERENCE,
        mode=config['mode'], whole_genome_bam_operations='prohibited', synthetic=config.get('_synthetic', False)))
    atomic_json(out / 'versions.lock.json', lock)
    runtimes = {r['id']: r for r in config.get('runtimes', [])}
    checks, upstream, downstream, source_files = [], [], [], []
    import_checks = []
    ref = config.get('reference', {})
    reference_identity = []
    if ref.get('fasta'):
        reference_identity.append(fast_fingerprint(ref['fasta']))
        reference_identity.append(dict(path=ref['fasta']+'.fai', sha256=digest(ref['fasta']+'.fai')))
    if ref.get('aa_data_repo') and ref.get('build'):
        resource_root = Path(ref['aa_data_repo']) / ref['build']
        for p in sorted(resource_root.glob('*')):
            if p.is_file() and p.suffix not in ('.fa', '.fasta'):
                reference_identity.append(dict(path=str(p), sha256=digest(p)) if p.stat().st_size <= 16*1024**2 else fast_fingerprint(p))
    background = config.get('local_cnv', {}).get('background_bed')
    if config['mode'] == 'with_bam' and background and Path(background).is_file():
        reference_identity.append(dict(path=background, sha256=digest(background), kind='local_cnv_background_bed'))
    config['_reference_identity'] = reference_identity
    source_files.extend(reference_identity)
    for sample in config['samples']:
        sid, source = sample['sample_id'], sample['source_version']
        print(f'Input {sid}: source={source}, mode={config["mode"]}', flush=True)
        try:
            build = config.get('reference', {}).get('build', sample['reference_build'])
            if build != sample['reference_build']:
                raise ValueError('Reference build differs between sample and configuration')
            seeds_path = sample.get('seeds_bed') or sample['source_result'].get('seeds_bed')
            seeds = read_seed_bed(seeds_path, lengths) if seeds_path else []
            sample['seeds'] = seeds
            sample['upstream_seed_ids'] = [s['seed_id'] for s in seeds]
            if seeds_path:
                sample['seeds_bed'] = seeds_path
                source_files.append(dict(path=seeds_path, sha256=digest(seeds_path)))
            specs = [(source, sample['source_result'])] + [(r['version'], r) for r in sample['comparison_results']]
            for vid, spec in specs:
                if check_only:
                    from .benchmark_import import inspect_import_format
                    entry = entries[vid]
                    categories = entry.get('supported_labels', entry.get('supported_classes', []))
                    result = inspect_import_format(spec, categories)
                    import_checks.append(dict(sample_id=sid, version=vid, **result))
                try:
                    if spec.get('reference_build', sample['reference_build']) != sample['reference_build']:
                        raise ValueError('Imported result reference_build differs from the sample declaration')
                    if not Path(spec.get('aa_dir', '')).is_dir() or not Path(spec.get('ac_dir', '')).is_dir():
                        raise ValueError('Imported AA/AC directory is missing')
                    if len(list(Path(spec['aa_dir']).glob('*_summary.txt'))) != 1:
                        raise ValueError('Expected one AA summary for each imported result')
                    artifact_total = 0
                    for path, _ in _artifacts(spec.get('aa_dir', ''), spec.get('ac_dir')):
                        if path.is_symlink() or not path.is_file():
                            raise ValueError('Only regular small result files can be imported')
                        artifact_total += path.stat().st_size
                        if path.stat().st_size > MAX_FILE_BYTES or artifact_total > MAX_TOTAL_BYTES:
                            raise ValueError('Result artifact size limit exceeded')
                        source_files.append(dict(path=str(path.resolve()), sha256=digest(path)))
                    for key in ('seeds_bed', 'stages_json'):
                        if spec.get(key) and Path(spec[key]).is_file():
                            source_files.append(dict(path=spec[key], sha256=digest(spec[key])))
                except (OSError, ValueError, KeyError) as exc:
                    if spec is sample['source_result']:
                        raise
                    spec['_input_error'] = str(exc)
            if config['mode'] == 'with_bam':
                for key in ('baseline_cns', 'source_seed_bed', 'baseline_provenance'):
                    if sample.get(key) and Path(sample[key]).is_file():
                        source_files.append(dict(path=sample[key], sha256=digest(sample[key]), kind=key))
                sample['_input_identities'] = _validate_bam(sample, config, lengths)
                source_files.extend(sample['_input_identities'])
                policy = dict(config.get('io_policy', {}))
                policy['alignments'] = {sample[k]: sample[i] for k, i in
                    (('alignment', 'index'), ('normal_alignment', 'normal_index')) if sample.get(k)}
                policy_path = out / 'policies' / (sid + '.json')
                atomic_json(policy_path, policy)
                sample['_policy_path'] = str(policy_path)
            checks.append(dict(sample_id=sid, status='ok', seed_count=len(seeds),
                reference_build=build, reference_identity='declared_build_and_contig_lengths_only',
                reason='No whole-genome BAM validation performed'))
        except Exception as exc:
            checks.append(dict(sample_id=sid, status='failed', reason=str(exc), stage='input_check'))
            sample['_input_failed'] = str(exc)
    input_identity = _signature(source_files)
    identity_file = out / 'input_identity.json'
    if identity_file.exists() and json.loads(identity_file.read_text())['hash'] != input_identity:
        raise ValueError('Input BAM/seed identity changed; choose a new output directory')
    atomic_json(identity_file, dict(hash=input_identity, sources=source_files))
    atomic_json(out / 'input_check.json', checks)
    table(out / 'input_check.tsv', checks)
    if check_only:
        atomic_json(out / 'import_format_check.json', import_checks)
        table(out / 'import_format_check.tsv', import_checks)
        atomic_json(out / 'config.resolved.json', config)
        return dict(check_only=True, input_check=checks, import_format_check=import_checks)
    verified = {}
    for sample in config['samples']:
        sid, source = sample['sample_id'], sample['source_version']
        versions = list(dict.fromkeys([source, *config['target_versions']]))
        if sample.get('_input_failed'):
            for vid in versions:
                r = empty_record(sid, '__input_check__', vid, 'failed', 'input_check', sample['_input_failed'])
                r.update(analysis_kind='aa_ac', source_version=source, comparison_eligible=False)
                downstream.append(r)
            continue
        seeds = sample['seeds']
        sampling = config.get('sampling', {})
        max_seeds = int(sampling.get('max_seeds', 10))
        if max_seeds < 1:
            raise ValueError('sampling.max_seeds must be positive')
        rng = random.Random(str(sampling.get('random_seed', 20261117)) + ':' + sid)
        chosen = sorted(rng.sample(seeds, min(max_seeds, len(seeds))), key=lambda x: x['seed_id'])
        sample['selected_seed_ids'] = [s['seed_id'] for s in chosen]
        atomic_json(out / 'sampling' / (sid + '.json'), dict(all_seeds=seeds, selected=chosen,
                    rule='uniform_without_replacement_before_results', parameters=sampling))
        imported = {}
        specs = [(source, sample['source_result'])] + [(r['version'], r) for r in sample['comparison_results']]
        for vid, spec in specs:
            entry = entries[vid]
            if spec.get('_input_error'):
                r = empty_record(sid, '__whole_result__', vid, 'missing', 'input_check', spec['_input_error'])
                r.update(source_version=source, analysis_kind='aa_ac', comparison_eligible=False)
                imported[vid] = r
                continue
            spec = _source_spec(spec, entry)
            spec.setdefault('reference_build', sample['reference_build'])
            job = out / 'imports' / sid / vid
            # Imports are small and read afresh on resume to validate provenance.
            dest = _new_attempt(job, create=False)
            r = import_result(spec, sid, vid, dest, entry.get('supported_labels', entry.get('supported_classes', [])), lengths)
            r.update(source_version=source, analysis_kind='aa_ac',
                     component_identity_level=spec['component_identity_level'],
                     reference_scope='declared_reference_' + sample['reference_build'])
            imported[vid] = r
        baseline = imported[source]
        if config['mode'] == 'results_only':
            for vid in versions:
                if vid in imported:
                    r = imported[vid]
                else:
                    r = empty_record(sid, '__whole_result__', vid, 'missing', 'result_reading', 'No existing AA result for this version; AA cannot be rerun without BAM')
                    r.update(analysis_kind='aa_ac', comparison_eligible=False, source_version=source)
                r['seed_id'] = '__whole_result__'
                downstream.append(r)
            # Fixed AA graph identity: every AC-only comparison includes a source anchor.
            anchor = copy.deepcopy(baseline)
            anchor.update(analysis_kind='ac_only', seed_id='__whole_result__', source_aa_version=baseline.get('source_aa_version'))
            downstream.append(anchor)
            for vid in versions:
                if vid == source:
                    continue
                if vid not in runtimes:
                    r = empty_record(sid, '__whole_result__', vid, 'missing', 'ac', 'No configured AC runtime')
                    r.update(analysis_kind='ac_only', comparison_eligible=False)
                else:
                    r = reclassify_cached(baseline, runtimes[vid], entries[vid], config, out / 'reclassifications' / sid / vid)
                r.update(source_version=source, seed_id='__whole_result__', source_aa_version=baseline.get('source_aa_version'))
                r['reference_scope'] = baseline.get('reference_scope', 'declared_reference_' + sample['reference_build'])
                downstream.append(r)
        else:
            # The per-sample fields are top-level overrides. Copying every
            # manifest sample here would make cohort preparation quadratic.
            local_config = dict(config)
            local_config['_io_policy_path'] = sample['_policy_path']
            local_config['_cnv_cache_dir'] = str(out / 'cnv_cache')
            for vid in versions:
                print(f'  {sid} / {vid}: local upstream and selected-seed downstream', flush=True)
                runtime, entry = runtimes.get(vid), entries[vid]
                if runtime:
                    if vid not in verified:
                        env = dict(os.environ, AA_DATA_REPO=config['reference']['aa_data_repo'], PYTHONDONTWRITEBYTECODE='1')
                        try:
                            verified[vid] = verify_runtime(runtime, entry, env)
                        except Exception as exc:
                            verified[vid] = dict(ready=False, diagnostics=[str(exc)], id=vid)
                        atomic_json(out / 'runtime_check.json', verified)
                    aa_verified = verified[vid].get('aa', {})
                    if aa_verified.get('version') == entry['aa']['version'] and aa_verified.get('sha256'):
                        upstream_attempt = _new_attempt(out / 'upstream_jobs' / sid / vid)
                        upstream.extend(run_local_upstream(sample, runtime, entry, local_config, upstream_attempt))
                    else:
                        for seed in seeds or [dict(seed_id='__no_seeds__')]:
                            upstream.append(dict(sample_id=sid, seed_id=seed['seed_id'], version=vid, source_version=source,
                                status='unsupported', reason='AA seed-filter source/runtime could not be verified',
                                comparison_eligible=False, calibration_status='unverified', cnvkit_version=runtime.get('cnvkit_version')))
                    if chosen:
                        joint = run_bounded_downstream(sample, chosen, runtime, entry, verified[vid], local_config, out / 'downstream_jobs' / sid / vid)
                        if vid == source:
                            bridge = _bridge(baseline, joint, chosen)
                            sample['source_bridge_status'] = 'passed' if bridge and all(r['passed'] for r in bridge) else 'failed'
                            atomic_json(out / 'source_bridge' / (sid + '.json'), bridge)
                        for seed in chosen:
                            downstream.append(project_record(joint, seed))
                    else:
                        r = empty_record(sid, '__no_seeds__', vid, 'no_seeds', 'seed_processing', 'No supplied candidates; no BAM discovery permitted')
                        r.update(analysis_kind='aa_ac', source_version=source, comparison_eligible=False)
                        downstream.append(r)
                else:
                    if vid == source:
                        sample['source_bridge_status'] = 'unverified'
                    for seed in seeds or [dict(seed_id='__no_seeds__')]:
                        upstream.append(dict(sample_id=sid, seed_id=seed['seed_id'], version=vid, source_version=source,
                            status='missing' if seeds else 'no_seeds', reason='No configured runtime', comparison_eligible=False,
                            calibration_status='unverified', cnvkit_version=None))
                    for seed in chosen or [dict(seed_id='__no_seeds__')]:
                        r = empty_record(sid, seed['seed_id'], vid, 'missing', 'aa', 'No configured AA+AC runtime')
                        r.update(analysis_kind='aa_ac', source_version=source, comparison_eligible=False)
                        downstream.append(r)
        atomic_json(out / 'downstream_records.json', downstream)
        atomic_json(out / 'upstream_records.json', upstream)
    atomic_json(out / 'config.resolved.json', config)
    atomic_json(out / 'downstream_records.json', downstream)
    atomic_json(out / 'upstream_records.json', upstream)
    table(out / 'run_metadata.tsv', [dict(program_version=PROGRAM_VERSION, created_utc=now(),
        reference_version=FIXED_REFERENCE, mode=config['mode'], manifest_hash=run_identity,
        whole_genome_bam_operations='prohibited', **s) for s in source_files])
    summary = write_benchmark_reports(out, downstream, upstream, config)
    atomic_json(out / 'report_inputs.sha256.json', {name: digest(out / name) for name in
        ('downstream_records.json', 'upstream_records.json', 'config.resolved.json')})
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description='AmpliconSuite Cross-Version Concordance: existing results, bounded local BAM operations')
    parser.add_argument('--version', action='version', version=PROGRAM_VERSION)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('check', 'run'):
        p = sub.add_parser(name)
        p.add_argument('--manifest', required=True)
        p.add_argument('--out', required=True)
    p = sub.add_parser('report'); p.add_argument('--out', required=True)
    p = sub.add_parser('demo'); p.add_argument('--out', required=True)
    p = sub.add_parser('prepare-panel', help='Fetch verified upstream sources without reading samples')
    p.add_argument('--out', required=True)
    p.add_argument('--versions', nargs='+')
    p.add_argument('--lock')
    sub.add_parser('versions')
    args = parser.parse_args(argv)
    try:
        if args.command == 'prepare-panel':
            from .versions import prepare_panel
            print('Saved ' + str(prepare_panel(args.out, args.versions, args.lock)))
            return 0
        if args.command == 'versions':
            from .benchmark_cnv import load_builtin_map
            lock = load_lock()
            recipe_map = load_builtin_map()['versions']
            print(json.dumps(dict(reference_version=lock['reference_version'], versions=[
                {**{k: (v['version'] if isinstance(v, dict) and 'version' in v else v)
                    for k, v in e.items() if k in ('id', 'suite', 'aa', 'ac', 'bfb', 'source_status', 'runtime_status')},
                 'cnvkit_recipe': recipe_map.get(e['id'], {}).get('cnvkit_version')}
                for e in lock['versions'] if e.get('suite')]), indent=2))
            return 0
        if args.command == 'demo':
            from .benchmark_demo import demo
            demo(args.out)
            return 0
        if args.command == 'report':
            from .benchmark_report import write_benchmark_reports
            out = Path(args.out).resolve()
            run_version = json.loads((out / 'provenance.json').read_text()).get('program_version')
            if run_version != PROGRAM_VERSION:
                raise ValueError('Use the original program version to redraw this report, or run the manifest into a new output directory')
            expected = json.loads((out / 'report_inputs.sha256.json').read_text())
            if not all((out / name).is_file() and digest(out / name) == value for name, value in expected.items()):
                raise ValueError('Report inputs changed after execution; rerun the same manifest to revalidate, or start a new analysis')
            write_benchmark_reports(out, json.loads((out / 'downstream_records.json').read_text()),
                json.loads((out / 'upstream_records.json').read_text()), json.loads((out / 'config.resolved.json').read_text()))
        else:
            config, lock, entries = load_manifest(args.manifest)
            result = run_benchmark(config, lock, entries, args.out, args.command == 'check')
            if args.command == 'check' and any(r['status'] != 'ok' for r in
                    result['input_check'] + result['import_format_check']):
                return 2
        print('Saved concordance results: ' + str(Path(args.out).resolve()))
        return 0
    except KeyboardInterrupt:
        print('Interrupted. Completed results retained; repeat command to resume.', file=sys.stderr)
        return 130
    except Exception as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
