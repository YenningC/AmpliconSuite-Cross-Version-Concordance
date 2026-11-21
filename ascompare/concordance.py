"""AmpliconSuite Cross-Version Concordance: local, paired version comparison."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import random
import re
import sys
from collections import Counter
from . import __version__
from .common import atomic_json, digest, now, table
from .inputs import safe_id
from .parsing import merge_intervals
from .runner import empty_record, execute, local_reference, check_ac_runtime, prepare_ac_input
from .versions import load_lock, verify_runtime
from .concordance_policy import DEFAULT_MAX_UNIQUE_BP, guarded_command
from .concordance_aa_stats import seed_stage_measurement, stage_measurement
from .concordance_timing import TimingSession, measure, timed_call, record_native_aa

PROGRAM_VERSION = __version__
FIXED_REFERENCE = 'suite-1.6.2'
PATH_KEYS = {'alignment', 'index', 'normal_alignment', 'normal_index', 'seeds_bed',
             'source_seed_bed', 'baseline_cns', 'aa_dir', 'ac_dir', 'run_dir',
             'result_dir', 'stages_json', 'baseline_provenance', 'fasta', 'aa_data_repo', 'background_bed',
             'aa_script', 'ac_script', 'suite_script', 'cnvkit_script', 'version_lock', 'reuse_aa_from',
             'reuse_global_cnv_from', 'recover_aa_from'}


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
        if k == 'cnvkit_runtimes':
            from .local_cnv_comparator import resolve_registry
            result[k] = resolve_registry(v, parent)
        elif v and (k in PATH_KEYS or k.endswith('_python')):
            p = Path(v).expanduser()
            p = parent / p if not p.is_absolute() else p
            result[k] = str(p.absolute() if k.endswith('_python') else p.resolve())
        elif k in {'extra_path', 'cnvkit_extra_path'}:
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


def cnv_scope(config):
    # An old manifest which explicitly configured bounded windows keeps its
    # scope. This prevents an upgrade from silently scanning those BAMs.
    return config.get('cnv_scope', 'local' if 'local_cnv' in config else 'global')


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
    global_run = config['mode'] == 'with_bam' and cnv_scope(config) == 'global'
    if config['mode'] == 'with_bam':
        config['cnv_scope'] = cnv_scope(config)
        if config['cnv_scope'] not in ('local', 'global'):
            raise ValueError('cnv_scope must be local or global')
        if global_run:
            if 'local_cnv' in config:
                raise ValueError('Remove local_cnv window settings when choosing global discovery')
            config.setdefault('parameter_policy', 'suite_defaults')
            gsettings = config.get('global_cnv', {})
            if not isinstance(gsettings, dict) or set(gsettings) - {'processes', 'timeout_seconds'}:
                raise ValueError('global_cnv accepts processes and timeout_seconds')
            if any(type(v) is not int or v < 1 for v in gsettings.values()):
                raise ValueError('global_cnv processes and timeout_seconds must be positive integers')
            fraction = config.get('seed_match_fraction', 1.0)
            if type(fraction) not in (int, float) or not 0 < fraction <= 1:
                raise ValueError('seed_match_fraction must be in (0, 1]')
            cap = config.get('sampling', {}).get('max_seeds', 10)
            if cap != 'all' and (type(cap) is not int or cap < 1):
                raise ValueError('sampling.max_seeds must be a positive integer or "all"')
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
        if global_run:
            archives = s.get('global_cnv_results', [])
            if not isinstance(archives, list) or not all(
                    isinstance(a, dict) and set(a) <= {'result_dir','version','cns_sha256'} and a.get('result_dir')
                    for a in archives):
                raise ValueError('global_cnv_results must list historical Suite result_dir entries')
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
        if not (s['source_result'].get('aa_dir') and s['source_result'].get('ac_dir')
                or global_run and s['source_result'].get('result_dir')):
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
            if not seed and not global_run:
                raise ValueError('Existing final Suite seed BED required; no automatic genome-wide discovery')
            if seed:
                s['seeds_bed'] = seed
            if bool(s.get('normal_alignment')) != bool(s.get('normal_index')):
                raise ValueError('normal_alignment and normal_index must be supplied together')
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
    from .local_cnv_comparator import is_exact_version
    extra_cnv = config.get('local_cnv', {}).get('comparison_versions', [])
    if not isinstance(extra_cnv, list) or not all(is_exact_version(v) for v in extra_cnv):
        raise ValueError('local_cnv.comparison_versions requires a list of exact CNVkit versions')
    if config.get('threads', 1) != 1:
        raise ValueError('The bounded concordance currently executes serially')
    parameter_policy = config.get('parameter_policy', 'custom')
    if parameter_policy not in ('suite_defaults', 'custom'):
        raise ValueError('parameter_policy must be suite_defaults or custom')
    for key, allowed in (('aa_parameters', {'downsample', 'runmode'}), ('ac_parameters', {'min_size', 'min_flow'})):
        if not isinstance(config.get(key, {}), dict) or set(config.get(key, {})) - allowed:
            raise ValueError('Unsupported parameters in ' + key + '; no parameters are silently ignored')
    if parameter_policy == 'suite_defaults':
        if config.get('aa_parameters') or config.get('ac_parameters') or any(
                key in config for key in ('ac_min_size', 'ac_min_flow', 'full_suite')):
            raise ValueError('suite_defaults uses official AA/AC parameters; remove manifest overrides')
        for key, default in (('cngain', 4.5), ('cnsize_min', 50000)):
            if key in config.get('local_cnv', {}) and config['local_cnv'][key] != default:
                raise ValueError(f'suite_defaults requires the official seed-filter default {key}={default}')
        for runtime in config.get('runtimes', []):
            if 'PYTHONHASHSEED' in runtime.get('env', {}):
                raise ValueError('suite_defaults inherits the Python hash environment; remove runtime PYTHONHASHSEED')
            if config['mode'] == 'with_bam' and 'solver' in runtime:
                raise ValueError('suite_defaults selects the official AA solver; remove runtime solver overrides')
    if config.get('aa_parameters', {}).get('runmode', 'CYCLES') not in ('CYCLES', 'FULL'):
        raise ValueError('aa_parameters.runmode must be CYCLES or FULL')
    if 'reuse_aa_from' in config:
        previous = config['reuse_aa_from']
        if (config['mode'] != 'with_bam' or not isinstance(previous, str) or not previous
                or not (Path(previous) / '.component_cache' / 'aa').is_dir()):
            raise ValueError('reuse_aa_from requires an existing with_bam output with a completed AA cache')
    if 'reuse_global_cnv_from' in config:
        origin = config['reuse_global_cnv_from']
        if not global_run or not isinstance(origin, str) or not (Path(origin) / 'global_cnv_cache').is_dir():
            raise ValueError('reuse_global_cnv_from needs a prior global with_bam output')
    if 'recover_aa_from' in config or 'aa_rescue_unique_bp' in config:
        origin = config.get('recover_aa_from')
        if not global_run or not isinstance(origin, str) or not (Path(origin) / 'report_inputs.sha256.json').is_file():
            raise ValueError('recover_aa_from needs a completed global with_bam report directory')
        policy = config.get('io_policy', {})
        limit = config.get('aa_rescue_unique_bp')
        if (type(limit) is not int or not int(policy.get('max_unique_bp', DEFAULT_MAX_UNIQUE_BP))
                < limit <= int(policy.get('max_total_fetch_bp', 500_000_000))):
            raise ValueError('aa_rescue_unique_bp must exceed the original unique limit and fit within max_total_fetch_bp')
    import math
    for key, default in (('min_size', 50000), ('min_flow', 1.0)):
        value = float(config.get('ac_parameters', {}).get(key, default))
        if not math.isfinite(value) or value <= 0:
            raise ValueError('AC parameters must be positive finite numbers')
    ds = float(config.get('aa_parameters', {}).get('downsample', -1))
    if not math.isfinite(ds) or ds != -1 and ds <= 0:
        raise ValueError('AA downsample must be -1 or positive finite coverage')
    from .concordance_report import _configuration
    _configuration(config)  # Reject invalid analysis plans before launching component jobs.
    return config, lock, entries


def _validate_bam(sample, config, lengths):
    import pysam
    identities = []
    discovery_lengths = {}
    if config.get('cnv_scope') == 'global':
        ref = config['reference']
        build = ref['build']
        genome = Path(ref['aa_data_repo']) / build
        sizes = genome / (build + '_noAlt.fa.fai')
        for line in sizes.read_text().splitlines():
            if line.strip():
                chrom, size, *_ = line.split()
                if chrom in discovery_lengths or lengths.get(chrom) != int(size):
                    raise ValueError(f'Global discovery reference contig differs from FASTA: {chrom}')
                discovery_lengths[chrom] = int(size)
        if not discovery_lengths:
            raise ValueError('Global discovery contig list is empty')
        if not sample.get('normal_alignment') or build == 'GRCh38_viral':
            if not (genome / (build + '_cnvkit_filtered_ref.cnn')).is_file():
                raise ValueError('Missing Suite tumor-only CNVkit reference for ' + build)
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
            for chrom, size in discovery_lengths.items():
                if found.get(chrom) != size:
                    raise ValueError(f'{key} lacks the expected global discovery contig: {chrom}')
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
    for name in ('stages.json', 'ac.input', 'ac_summary_map.txt', 'aa_background.json', 'reference/coverage.stats'):
        if (run / name).is_file():
            files[str(run / name)] = digest(run / name)
    atomic_json(Path(job) / 'cache.json', dict(signature=signature, record=record, files=files))


def _new_attempt(job, create=True):
    job = Path(job)
    job.mkdir(parents=True, exist_ok=True)
    existing = [int(p.name[8:]) for p in job.glob('attempt_*') if p.name[8:].isdigit()]
    p = job / f'attempt_{max(existing, default=0)+1:04d}'
    if create:
        p.mkdir()
    return p


def run_bounded_downstream(sample, seeds, runtime, entry, verified, config, out, measurement=None):
    """Run or integrity-check each exact AA and AC stage independently."""
    from .concordance_components import run_component_downstream
    return run_component_downstream(sample, seeds, runtime, entry, verified, config, out, measurement)


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


def run_concordance(config, lock, entries, out, check_only=False):
    from .concordance_import import import_result, reclassify_result
    from .concordance_cnv import run_local_upstream
    from .concordance_report import write_concordance_reports
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
        progress_path = out / 'sample_progress.json'
        previous_progress = progress_path.read_bytes() if progress_path.is_file() else None
        try:
            with TimingSession(out):
                return _run_locked(config, lock, entries, out, check_only)
        except (Exception, KeyboardInterrupt) as exc:
            try:
                if progress_path.is_file() and progress_path.read_bytes() != previous_progress:
                    progress = json.loads(progress_path.read_text())
                    for row in progress:
                        if row['state'] in {'checking', 'running'}:
                            row.update(state='interrupted' if isinstance(exc, KeyboardInterrupt) else 'aborted', updated_utc=now())
                    atomic_json(progress_path, progress)
                    table(out / 'sample_progress.tsv', progress)
            except (OSError, ValueError, TypeError, KeyError):
                pass  # Preserve the original failure even if progress cannot be saved.
            raise
    finally:
        handle.close()


def _run_locked(config, lock, entries, out, check_only):
    from .concordance_import import import_result, reclassify_cached, _artifacts, MAX_FILE_BYTES, MAX_TOTAL_BYTES
    from .local_cnv_comparator import run_local_cnv_comparators
    from .concordance_cnv import run_local_upstream
    from .concordance_report import write_concordance_reports
    config = copy.deepcopy(config)
    global_run = config['mode'] == 'with_bam' and config.get('cnv_scope') == 'global'
    bam_scope = 'global_CNVkit_only; bounded_seed_filter_and_AA' if global_run else 'prohibited'
    if config['mode'] == 'with_bam':
        from .concordance_cnv_source import identify_source_cnvkit
        for sample in config['samples']:
            identify_source_cnvkit(sample)
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
        mode=config['mode'], parameter_policy=config.get('parameter_policy', 'custom'),
        whole_genome_bam_operations=bam_scope, synthetic=bool(config.get('synthetic', config.get('_synthetic', False)))))
    atomic_json(out / 'versions.lock.json', lock)
    runtimes = {r['id']: r for r in config.get('runtimes', [])}
    checks, upstream, downstream, source_files = [], [], [], []
    local_cnv_records = []
    progress = [dict(order=i, sample_id=s['sample_id'], source_version=s['source_version'],
                     state='queued', updated_utc=now()) for i, s in enumerate(config['samples'], 1)]

    def record_progress(sample, state):
        item = next(r for r in progress if r['sample_id'] == sample['sample_id'])
        local = [r for r in downstream if r['sample_id'] == sample['sample_id']]
        item.update(state=state, updated_utc=now(),
                    result_status_counts=dict(Counter(r['status'] for r in local)))
        atomic_json(out / 'sample_progress.json', progress)
        table(out / 'sample_progress.tsv', progress)

    atomic_json(out / 'sample_progress.json', progress)
    table(out / 'sample_progress.tsv', progress)
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
    for sample_index, sample in enumerate(config['samples'], 1):
        sid, source = sample['sample_id'], sample['source_version']
        record_progress(sample, 'checking')
        print(f'[{sample_index}/{len(config["samples"])}] Input {sid}: source={source}, mode={config["mode"]}', flush=True)
        with measure('input_check', sample_id=sid, version=source) as input_clock:
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
                    if global_run and spec is sample['source_result'] and spec.get('result_dir'):
                        if not Path(spec['result_dir']).is_dir():
                            raise ValueError('Source Suite result directory is missing')
                        if not (spec.get('aa_dir') and spec.get('ac_dir')):
                            spec['_input_error'] = 'Source AA/AC files unavailable for historical replay comparison'
                            continue
                    if check_only:
                        from .concordance_import import inspect_import_format
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
                        if spec is sample['source_result'] and not global_run:
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
                if config['mode'] == 'with_bam':
                    if global_run:
                        checks[-1].update(cnv_scope='global', historical_seed_bed_supplied=bool(seeds_path),
                                          historical_result_issue=sample['source_result'].get('_input_error'))
                    checks[-1].update(source_cnvkit_version=sample.get('source_cnvkit_version'),
                                      source_cnvkit_origin=sample['source_cnvkit_evidence']['origin'],
                                      source_cnvkit_issue=sample.get('_source_cnvkit_issue'))
                record_progress(sample, 'ready')
            except Exception as exc:
                checks.append(dict(sample_id=sid, status='failed', reason=str(exc), stage='input_check'))
                sample['_input_failed'] = str(exc)
                input_clock.status = 'failed'
                record_progress(sample, 'input_failed')
    input_identity = _signature(source_files)
    identity_file = out / 'input_identity.json'
    if identity_file.exists() and json.loads(identity_file.read_text())['hash'] != input_identity:
        raise ValueError('Input BAM/seed identity changed; choose a new output directory')
    atomic_json(identity_file, dict(hash=input_identity, sources=source_files))
    atomic_json(out / 'input_check.json', checks)
    table(out / 'input_check.tsv', checks)
    aa_recovery = None
    if config.get('recover_aa_from'):
        if all(row['status'] == 'ok' for row in checks):
            from .concordance_aa_recovery import verify_recovery_source
            aa_recovery = verify_recovery_source(config, out)
        elif not check_only:
            raise ValueError('AA recovery requires all original inputs to pass validation; see input_check.tsv')
    if check_only:
        atomic_json(out / 'import_format_check.json', import_checks)
        table(out / 'import_format_check.tsv', import_checks)
        atomic_json(out / 'config.resolved.json', config)
        return dict(check_only=True, input_check=checks, import_format_check=import_checks)
    verified = {}
    def verified_runtime(vid, sid):
        runtime, entry = runtimes[vid], entries[vid]
        if vid not in verified:
            env = dict(os.environ, AA_DATA_REPO=config['reference']['aa_data_repo'], PYTHONDONTWRITEBYTECODE='1')
            try:
                if config.get('parameter_policy') == 'suite_defaults':
                    from .concordance_official_solver import official_solver
                    solver, reason = official_solver(runtime, entry, env)
                    official_solver_evidence = dict(requested=solver, reason=reason)
                    runtime = dict(runtime, solver=solver)
                    if entry['aa']['version'].startswith('1.6.'):
                        runtime['_allow_native_mosek_fallback'] = True
                    runtimes[vid] = runtime
                verified[vid] = timed_call('runtime_check', verify_runtime, runtime, entry, env,
                    timing=dict(sample_id=sid, version=vid))
                if config.get('parameter_policy') == 'suite_defaults':
                    verified[vid]['official_solver'] = official_solver_evidence
            except Exception as exc:
                verified[vid] = dict(ready=False, diagnostics=[str(exc)], id=vid)
            atomic_json(out / 'runtime_check.json', verified)
        return runtimes[vid]

    for sample_index, sample in enumerate(config['samples'], 1):
        sid, source = sample['sample_id'], sample['source_version']
        versions = list(dict.fromkeys([source, *config['target_versions']]))
        if sample.get('_input_failed'):
            if global_run:
                sample.update(upstream_seed_ids=[], selected_seed_ids=[], global_source_seed_calibration='unverified',
                              global_seed_sets={v:dict(status='failed', stage='input_check',
                              reason=sample['_input_failed'], seeds=[]) for v in versions})
                from .concordance_global_report import global_report_units
                up, down = global_report_units(sample, {}, versions)
                upstream.extend(up)
                downstream.extend(down)
            for vid in ([] if global_run else versions):
                r = empty_record(sid, '__input_check__', vid, 'failed', 'input_check', sample['_input_failed'])
                r.update(analysis_kind='aa_ac', source_version=source, comparison_eligible=False)
                downstream.append(r)
            if config['mode'] == 'with_bam' and not global_run:
                local_cnv_records.extend(timed_call('local_cnv_comparators', run_local_cnv_comparators,
                    sample, config, entries, runtimes, verified, [], out, timing=dict(sample_id=sid)))
            record_progress(sample, 'input_failed')
            atomic_json(out / 'local_cnv_comparator_records.json', local_cnv_records)
            atomic_json(out / 'downstream_records.json', downstream)
            continue
        record_progress(sample, 'running')
        print(f'[{sample_index}/{len(config["samples"])}] Analyzing {sid}; source={source}; {len(versions)} versions', flush=True)
        seeds = sample['seeds']
        sampling = config.get('sampling', {})
        if not global_run:
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
            r = timed_call('result_import', import_result, spec, sid, vid, dest,
                entry.get('supported_labels', entry.get('supported_classes', [])), lengths,
                timing=dict(sample_id=sid, version=vid))
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
                    r = timed_call('ac_reclassification', reclassify_cached, baseline, runtimes[vid], entries[vid],
                        config, out / 'reclassifications' / sid / vid, timing=dict(sample_id=sid, version=vid))
                r.update(source_version=source, seed_id='__whole_result__', source_aa_version=baseline.get('source_aa_version'))
                r['reference_scope'] = baseline.get('reference_scope', 'declared_reference_' + sample['reference_build'])
                downstream.append(r)
        elif global_run:
            from .concordance_global_workflow import run_global_sample
            observed_up, observed_down = run_global_sample(
                sample, baseline, versions, entries, config, out, runtimes, verified, verified_runtime,
                aa_recovery=aa_recovery)
            upstream.extend(observed_up)
            downstream.extend(observed_down)
        elif not chosen:
            reason = 'No supplied candidate regions; selected-region analysis is unevaluable'
            sample['source_bridge_status'] = 'unverified'
            for vid in versions:
                upstream.append(dict(sample_id=sid, seed_id='__no_seeds__', version=vid,
                    source_version=source, status='no_seeds', stage='seed_processing', reason=reason,
                    comparison_eligible=False, calibration_status='unverified', cnvkit_version=None))
                downstream.append(empty_record(sid, '__no_seeds__', vid, 'no_seeds', 'seed_processing',
                    reason, analysis_kind='aa_ac', source_version=source, comparison_eligible=False, detected=None))
            local_config, native_source_rows = config, []
        else:
            # The per-sample fields are top-level overrides. Copying every
            # manifest sample here would make cohort preparation quadratic.
            local_config = dict(config)
            local_config['_io_policy_path'] = sample['_policy_path']
            local_config['_cnv_cache_dir'] = str(out / 'cnv_cache')
            local_config['_component_cache_dir'] = str(out / '.component_cache')
            native_source_rows = []
            for vid in versions:
                print(f'  {sid} / {vid}: local upstream and selected-seed downstream', flush=True)
                runtime, entry = runtimes.get(vid), entries[vid]
                if runtime:
                    runtime = verified_runtime(vid, sid)
                    aa_verified = verified[vid].get('aa', {})
                    measurement = None
                    if aa_verified.get('version') == entry['aa']['version'] and aa_verified.get('sha256'):
                        upstream_attempt = _new_attempt(out / 'upstream_jobs' / sid / vid)
                        upstream_rows = run_local_upstream(sample, runtime, entry, local_config, upstream_attempt)
                        upstream.extend(upstream_rows)
                        if vid == source:
                            native_source_rows = upstream_rows
                        measurement = seed_stage_measurement(upstream_attempt, upstream_rows, sample, runtime,
                                                             verified[vid], local_config)
                    else:
                        for seed in seeds:
                            upstream.append(dict(sample_id=sid, seed_id=seed['seed_id'], version=vid, source_version=source,
                                status='unsupported', reason='AA seed-filter source/runtime could not be verified',
                                comparison_eligible=False, calibration_status='unverified', cnvkit_version=runtime.get('cnvkit_version')))
                    joint = timed_call('selected_seed_downstream', run_bounded_downstream,
                        sample, chosen, runtime, entry, verified[vid], local_config,
                        out / 'downstream_jobs' / sid / vid, measurement,
                        timing=dict(sample_id=sid, version=vid))
                    aa_cache = joint.get('component_cache', {}).get('aa', {})
                    if joint.get('run_dir'):
                        record_native_aa(joint['run_dir'], sid, vid,
                            role='producer_reference' if aa_cache.get('hit') else 'current_execution',
                            producer_id=aa_cache.get('producer_run_dir', ''))
                    if vid == source:
                        bridge = _bridge(baseline, joint, chosen)
                        sample['source_bridge_status'] = 'passed' if bridge and all(r['passed'] for r in bridge) else 'failed'
                        atomic_json(out / 'source_bridge' / (sid + '.json'), bridge)
                    for seed in chosen:
                        downstream.append(project_record(joint, seed))
                else:
                    if vid == source:
                        sample['source_bridge_status'] = 'unverified'
                    for seed in seeds:
                        upstream.append(dict(sample_id=sid, seed_id=seed['seed_id'], version=vid, source_version=source,
                            status='missing', reason='No configured runtime', comparison_eligible=False,
                            calibration_status='unverified', cnvkit_version=None))
                    for seed in chosen:
                        r = empty_record(sid, seed['seed_id'], vid, 'missing', 'aa', 'No configured AA+AC runtime')
                        r.update(analysis_kind='aa_ac', source_version=source, comparison_eligible=False)
                        downstream.append(r)
        if config['mode'] == 'with_bam' and not global_run:
            local_cnv_records.extend(timed_call('local_cnv_comparators', run_local_cnv_comparators,
                sample, local_config, entries, runtimes, verified, native_source_rows, out,
                timing=dict(sample_id=sid)))
        atomic_json(out / 'local_cnv_comparator_records.json', local_cnv_records)
        atomic_json(out / 'downstream_records.json', downstream)
        atomic_json(out / 'upstream_records.json', upstream)
        if global_run:
            states = sample['global_version_statuses']
            failed = [v for v, r in states.items() if r['status'] not in ('ok', 'not_detected')]
            progress[sample_index-1].update(completed_versions=len(states)-len(failed), total_versions=len(states),
                                            incomplete_versions=failed)
            state = 'failed' if len(failed) == len(states) else 'partial' if failed else 'finished'
            print(f'  {sid}: {len(states)-len(failed)}/{len(states)} Suite paths completed', flush=True)
            record_progress(sample, state)
        else:
            record_progress(sample, 'finished')
    atomic_json(out / 'config.resolved.json', config)
    atomic_json(out / 'downstream_records.json', downstream)
    atomic_json(out / 'upstream_records.json', upstream)
    table(out / 'run_metadata.tsv', [dict(program_version=PROGRAM_VERSION, created_utc=now(),
        reference_version=FIXED_REFERENCE, mode=config['mode'], manifest_hash=run_identity,
        whole_genome_bam_operations=bam_scope, **s) for s in source_files])
    atomic_json(out / 'local_cnv_comparator_records.json', local_cnv_records)
    summary = timed_call('report_generation', write_concordance_reports, out, downstream, upstream,
        config, local_cnv_records=local_cnv_records)
    atomic_json(out / 'report_inputs.sha256.json', {name: digest(out / name) for name in
        ('downstream_records.json', 'upstream_records.json', 'config.resolved.json', 'local_cnv_comparator_records.json')})
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description='AmpliconSuite Cross-Version Concordance: global CNVkit discovery and version-specific seed reconstruction')
    parser.add_argument('--version', action='version', version=PROGRAM_VERSION)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('check', 'run'):
        p = sub.add_parser(name)
        p.add_argument('--manifest', required=True)
        p.add_argument('--out', required=True)
    p = sub.add_parser('report'); p.add_argument('--out', required=True)
    p.add_argument('--from-run', help='Read verified records from a finished run and write a separate report to a new directory')
    p = sub.add_parser('demo'); p.add_argument('--out', required=True)
    p = sub.add_parser('import-samples', help='Convert a local sample TSV and shared JSON config into a batch manifest')
    p.add_argument('--samples', required=True)
    p.add_argument('--config', required=True)
    p.add_argument('--out', required=True)
    p = sub.add_parser('prepare-panel', help='Fetch verified upstream sources without reading samples')
    p.add_argument('--out', required=True)
    p.add_argument('--versions', nargs='+')
    p.add_argument('--lock')
    sub.add_parser('versions')
    args = parser.parse_args(argv)
    try:
        if args.command == 'import-samples':
            from .concordance_batch import import_samples
            config = import_samples(args.samples, args.config, args.out)
            print(f'Saved batch manifest: {Path(args.out).resolve()} ({len(config["samples"])} samples; serial analysis)')
            return 0
        if args.command == 'prepare-panel':
            from .versions import prepare_panel
            print('Saved ' + str(prepare_panel(args.out, args.versions, args.lock)))
            return 0
        if args.command == 'versions':
            from .concordance_cnv import load_builtin_map
            lock = load_lock()
            recipe_map = load_builtin_map()['versions']
            print(json.dumps(dict(reference_version=lock['reference_version'], versions=[
                {**{k: (v['version'] if isinstance(v, dict) and 'version' in v else v)
                    for k, v in e.items() if k in ('id', 'suite', 'aa', 'ac', 'bfb', 'source_status', 'runtime_status')},
                 'cnvkit_recipe': recipe_map.get(e['id'], {}).get('cnvkit_version')}
                for e in lock['versions'] if e.get('suite')]), indent=2))
            return 0
        if args.command == 'demo':
            from .concordance_demo import demo
            demo(args.out)
            return 0
        if args.command == 'report':
            from .concordance_report import write_concordance_reports
            out = Path(args.out).resolve()
            source = Path(args.from_run).resolve() if args.from_run else out
            original = json.loads((source / 'provenance.json').read_text())
            run_version = original.get('program_version')
            if run_version not in {PROGRAM_VERSION, '0.5.0.dev0', '0.4.4.dev2', '0.4.4.dev1', '0.4.4.dev0', '0.4.3', '0.4.2', '0.4.1', '0.4.0',
                                   '0.3.0', '0.3.1', '0.3.2'}:
                raise ValueError('Use the original program version to redraw this report, or run the manifest into a new output directory')
            names = ('downstream_records.json', 'upstream_records.json', 'config.resolved.json')
            if run_version in {PROGRAM_VERSION, '0.5.0.dev0', '0.4.4.dev2', '0.4.4.dev1', '0.4.4.dev0', '0.4.3', '0.4.2'}:
                names += ('local_cnv_comparator_records.json',)
            expected = json.loads((source / 'report_inputs.sha256.json').read_text())
            if set(expected) != set(names):
                raise ValueError('Incomplete report input checksum manifest')
            inputs = {}
            for name in names:
                data = (source / name).read_bytes()
                if hashlib.sha256(data).hexdigest() != expected[name]:
                    raise ValueError('Report inputs changed after execution; revalidate before replay')
                inputs[name] = json.loads(data)
            if args.from_run:
                if out == source or source in out.parents or out in source.parents:
                    raise ValueError('Choose a new report directory separate from the original run')
                out.mkdir(parents=True, exist_ok=False)
            write_concordance_reports(out, inputs['downstream_records.json'], inputs['upstream_records.json'],
                                      inputs['config.resolved.json'], save_records=False,
                                      local_cnv_records=inputs.get('local_cnv_comparator_records.json'))
            atomic_json(out / 'report_provenance.json', dict(
                source_directory=str(source), source_run_version=run_version, report_version=PROGRAM_VERSION,
                report_inputs_sha256=expected, algorithms_rerun=False, created_utc=now()))
        else:
            config, lock, entries = load_manifest(args.manifest)
            result = run_concordance(config, lock, entries, args.out, args.command == 'check')
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
