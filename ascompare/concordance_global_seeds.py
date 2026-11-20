"""Apply the exact locked Suite's whole-genome CNV seed selection."""
import json
import math
from pathlib import Path
import re
import shutil

from .common import atomic_json, digest
from .concordance_cache import signature
from .concordance_cnv import LocalFailure, _environment, _run_tool
from .concordance_component_identity import component_identity
from .concordance_official_solver import _suite_root
from .concordance_timing import measure
from .full_suite import read_cns_segment_count
from .runner import execute, local_reference
from .versions import _check_source


# Run in the locked Suite interpreter. Import the wrapper as a module to call
# only its filter functions; the full CLI would also check/scan the BAM.
_PREFILTER = """
import importlib.util
import json
import os
from pathlib import Path
import sys
script, cns, build, gain = sys.argv[1:]
sys.path.insert(0, str(Path(script).parent))
spec = importlib.util.spec_from_file_location('locked_suite', script)
suite = importlib.util.module_from_spec(spec)
spec.loader.exec_module(suite)
from paalib import cnv_prefilter
root = Path(os.environ['AA_DATA_REPO'])
suite.AA_REPO = str(root) + '/'
sizes = suite.get_ref_sizes(str(root / build / (build + '_noAlt.fa.fai')))
centromeres = suite.get_ref_centromeres(build)
bed = suite.convert_cnvkit_cns_to_bed('./', 'suite', cnsfile=cns)
filtered = cnv_prefilter.prefilter_bed(bed, build, centromeres, sizes, float(gain), '.')
Path('prefilter_result.json').write_text(json.dumps(dict(bed=Path(filtered).name)))
"""
FINAL_BED = 'suite_AA_CNV_SEEDS.bed'


def _native_seeds(path, lengths):
    from .concordance import read_seed_bed
    rows = sorted(read_seed_bed(path, lengths), key=lambda r: r['source_line'])
    lines = Path(path).read_text().splitlines()
    for row in rows:
        fields = lines[row['source_line'] - 1].split('\t')
        if len(fields) < 4 or fields[:3] != [row['chrom'], str(row['start']), str(row['end'])]:
            raise ValueError('Invalid native AA seed BED row ' + str(row['source_line']))
        row.update(source_fields=fields[3:], native_suite_bed=True)
    return rows


def _tumor_identity(sample):
    from .concordance import fast_fingerprint
    bam = fast_fingerprint(sample['alignment'])
    index = dict(path=str(Path(sample['index']).resolve()), sha256=digest(sample['index']))
    validated = {r['path']: r for r in sample['_input_identities']}
    if any(validated.get(identity['path'], {}).get(k) != v
           for identity in (bam, index) for k, v in identity.items() if k != 'method'):
        raise ValueError('Tumor BAM/index changed after preflight')
    return [bam, index]


def _cached_filter(attempt, key):
    try:
        meta = json.loads((attempt.parent / 'global_filter_cache.json').read_text())
        producer_name = meta['producer']
        if meta['key'] != key or not re.fullmatch(r'attempt_\d+', producer_name):
            return None
        producer = attempt.parent / producer_name
        if producer.is_symlink() or producer == attempt or not producer.is_dir():
            return None
        hashes = meta['sha256']
        for name in [FINAL_BED, 'upstream_provenance.json', 'seed_processing.log', 'prefilter.log']:
            if name not in hashes:
                return None
        for name, expected in hashes.items():
            path = producer / name
            if (path.resolve() != producer.resolve() / name or '..' in Path(name).parts
                    or Path(name).is_absolute() or not path.is_file() or digest(path) != expected):
                return None
        return producer
    except (KeyError, ValueError, OSError, TypeError):
        return None


def global_suite_seeds(sample, runtime, entry, config, cns, attempt):
    """Return all native, ordered seeds, including an empty successful outcome.

    The attempt directory is fresh per invocation. Only successful, unchanged
    seed filters under this same Suite job can be resumed.
    """
    attempt, cns = Path(attempt).resolve(), Path(cns).resolve()
    attempt.mkdir(parents=True, exist_ok=True)
    if any(attempt.iterdir()):
        raise ValueError('Global seed filter requires a fresh attempt directory')
    if cns.suffix != '.cns' or cns.name.endswith(('.call.cns', '_rescaled.cns', '.bintest.cns')):
        raise ValueError('Global Suite seed filtering requires native uncalled segmentation .cns')
    segment_count = read_cns_segment_count(cns)
    with cns.open() as fh:
        if fh.readline().strip().split('\t')[:5] != ['chromosome', 'start', 'end', 'gene', 'log2']:
            raise ValueError('Native Suite CNS conversion expects chromosome/start/end/gene/log2 in that order')
    suite_root = _suite_root(runtime, entry)
    suite = _check_source(suite_root, entry['suite'])
    suite_script = suite_root / entry['suite']['entrypoint']
    aa_root = Path(runtime['aa_script']).resolve()
    for _ in Path(entry['aa']['entrypoint']).parts:
        aa_root = aa_root.parent
    if aa_root / entry['aa']['entrypoint'] != Path(runtime['aa_script']).resolve():
        raise ValueError('AA entrypoint differs from the pinned seed helper checkout')
    _check_source(aa_root, entry['aa'])
    helper = Path(runtime['aa_script']).resolve().with_name('amplified_intervals.py')
    gain = float(config.get('local_cnv', {}).get('cngain', config.get('full_suite', {}).get('cngain', 4.5)))
    min_size_raw = float(config.get('local_cnv', {}).get('cnsize_min', config.get('full_suite', {}).get('cnsize_min', 50000)))
    if not math.isfinite(gain) or gain <= 0 or not min_size_raw.is_integer() or min_size_raw < 0:
        raise ValueError('Invalid Suite seed gain/minimum size')
    min_size = int(min_size_raw)
    official = config.get('parameter_policy') == 'suite_defaults'
    if official and (gain != 4.5 or min_size != 50000):
        raise ValueError('suite_defaults requires cngain=4.5 and cnsize_min=50000')
    build = sample['reference_build']
    policy = config['_io_policy_path']
    suite_python = runtime.get('suite_python') or runtime['ac_python']
    reference = config['reference']
    aa_env = _environment(runtime, runtime['aa_python'])
    suite_env = _environment(runtime, suite_python)
    verified = config.get('_seed_verified_runtime') or runtime.get('_seed_verified_runtime') or {}
    identity = dict(cns_sha256=digest(cns), sample_id=sample['sample_id'], build=build, version=runtime['id'],
        suite=suite, aa=component_identity('aa', runtime, entry, verified, aa_env),
        suite_python=str(suite_python), suite_python_sha256=digest(suite_python),
        filter_source_sha256=digest(__file__), guard_sha256=digest(Path(__file__).with_name('concordance_policy.py')),
        dependencies=signature(verified.get('ac', {}).get('packages', {})),
        input_identities=_tumor_identity(sample), reference=config.get('_reference_identity', []),
        policy_sha256=digest(policy), parameters=dict(gain=gain, min_size=min_size, official=official))
    key = signature(identity)
    # Suite's interpreter normally equals AC's (the source of its package inventory).
    reusable = (identity['reference'] and identity['input_identities']
                and str(Path(suite_python).absolute()) == str(Path(runtime['ac_python']).absolute())
                and verified.get('ac', {}).get('packages', {}).get('returncode') == 0
                and identity['aa']['packages'].get('returncode') == 0)
    producer = _cached_filter(attempt, key) if reusable else None
    if producer:
        with measure('seed_filter_cache', cache_status='hit', producer_id=str(producer)):
            shutil.copyfile(producer / FINAL_BED, attempt / FINAL_BED)
            atomic_json(attempt / 'upstream_provenance.json', dict(version=runtime['id'], key=key,
                cache_hit=True, producer=str(producer), seed_bed_sha256=digest(attempt / FINAL_BED)))
            return _native_seeds(attempt / FINAL_BED, config.get('_reference_lengths'))

    repo = str(local_reference(reference['aa_data_repo'], attempt / 'private_AA_DATA_REPO'))
    aa_env.update(AA_DATA_REPO=repo)
    suite_env.update(AA_DATA_REPO=repo, MPLCONFIGDIR=str(attempt / 'matplotlib_cache'))
    timeout = config.get('job_timeout_seconds', 3600)
    command = [str(suite_python), '-B', '-I', '-c', _PREFILTER, str(suite_script), str(cns), build, str(gain)]
    prelog = attempt / 'prefilter.log'
    with measure('suite_cnv_prefilter', measurement_kind='subprocess', log_path=str(prelog)) as timing:
        rc = execute(command, attempt, suite_env, prelog, timeout=timeout)
        timing.status = 'complete' if rc == 0 else 'failed'
    if rc:
        raise LocalFailure('failed', 'seed_processing', 'Native Suite CNV prefilter failed; see ' + str(prelog))
    filtered = json.loads((attempt / 'prefilter_result.json').read_text())['bed']
    if Path(filtered).name != filtered or not (attempt / filtered).is_file():
        raise LocalFailure('missing', 'seed_processing', 'Native Suite prefilter BED absent')
    # A stable relative filename is significant: amplified_intervals records
    # its own input filename in the last column of every native AA seed.
    args = ['--ref', build, '--bed', filtered, '--bam', str(Path(sample['alignment']).resolve()),
            '--gain', str(gain), '--cnsize_min', str(min_size), '--out', FINAL_BED[:-4]]
    helper_text = helper.read_text()
    if '--logfile' in helper_text and '--logfile' in suite_script.read_text():
        args.extend(['--logfile', 'native_seed_filter.log'])
    if not official and '--no_cstats' in helper_text:
        args.append('--no_cstats')
    logfile = _run_tool(runtime['aa_python'], helper, args, policy, attempt, aa_env, 'seed_processing', timeout)
    log = logfile.read_text(errors='replace')
    native_log = attempt / 'native_seed_filter.log'
    if native_log.is_file():
        log += '\n' + native_log.read_text(errors='replace')
    if 'Traceback (most recent call last)' in log or 'zero division error' in log.lower():
        raise LocalFailure('failed', 'seed_processing', 'Native AA seed filter logged an error; see ' + str(logfile))
    for audit in attempt.glob('bam_access.*.json'):
        blocked = json.loads(audit.read_text()).get('blocked')
        if blocked:
            raise LocalFailure('failed', 'seed_processing', 'AA seed filter BAM access incomplete: ' + str(blocked[0]))
    if not (attempt / FINAL_BED).is_file():
        raise LocalFailure('missing', 'seed_processing', 'Native AA seed filtering did not produce ' + FINAL_BED)
    seeds = _native_seeds(attempt / FINAL_BED, config.get('_reference_lengths'))
    if digest(cns) != identity['cns_sha256'] or _tumor_identity(sample) != identity['input_identities']:
        raise ValueError('Native seed inputs changed during filtering')
    provenance = attempt / 'upstream_provenance.json'
    atomic_json(provenance, dict(version=runtime['id'], key=key, identity=identity, cache_hit=False,
        cns=str(cns), segment_count=segment_count, prefilter_bed=filtered, seed_count=len(seeds),
        seed_bed_sha256=digest(attempt / FINAL_BED), io_policy=policy,
        amplified_intervals=dict(path=str(helper), sha256=digest(helper), args=args)))
    if reusable:
        files = [FINAL_BED, filtered, 'suite_CNV_CALLS.bed', 'seed_processing.log', 'prefilter.log',
                 'prefilter_result.json', provenance.name]
        if native_log.is_file():
            files.append(native_log.name)
        if (attempt / 'private_AA_DATA_REPO/coverage.stats').is_file():
            files.append('private_AA_DATA_REPO/coverage.stats')
        files.extend(p.name for p in attempt.glob('bam_access.*.json'))
        atomic_json(attempt.parent / 'global_filter_cache.json', dict(key=key, producer=attempt.name,
                    sha256={f: digest(attempt / f) for f in files}))
    return seeds


def global_seed_measurement(attempt, seeds, sample, runtime, verified, config):
    from .concordance_aa_stats import native_seed_measurement
    attempt = Path(attempt).resolve()
    if not seeds:
        return None
    provenance = json.loads((attempt / 'upstream_provenance.json').read_text())
    producer = _cached_filter(attempt, provenance['key']) if provenance.get('cache_hit') else attempt
    if (not producer or _native_seeds(producer / FINAL_BED, config.get('_reference_lengths')) != seeds
            or digest(producer / FINAL_BED) != provenance['seed_bed_sha256']):
        return None
    return native_seed_measurement(producer, sample, runtime, verified, config)
