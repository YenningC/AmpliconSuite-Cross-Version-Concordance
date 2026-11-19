"""Bounded local CNVkit and native AA seed-filter comparison.

The endpoint is retention of an existing candidate, never confirmed amplification.
Local bins and local centering change the context of CNV estimation; therefore an
uncalibrated run cannot establish equivalence to a whole-genome pipeline.
"""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil

from .concordance_policy import guarded_command
from .common import atomic_json, digest, now
from .runner import execute, local_reference


class LocalFailure(Exception):
    def __init__(self, status, stage, message):
        super().__init__(message)
        self.status, self.stage = status, stage


def load_builtin_map():
    return json.loads(Path(__file__).with_name('builtin_cnvkit_map.json').read_text())


def _merge(intervals):
    result = []
    for chrom, start, end in sorted(intervals):
        if result and result[-1][0] == chrom and start <= result[-1][2]:
            result[-1] = (chrom, result[-1][1], max(end, result[-1][2]))
        else:
            result.append((chrom, start, end))
    return result


def _subtract(intervals, excluded):
    result = []
    for chrom, start, end in intervals:
        pieces = [(start, end)]
        for c, a, b in excluded:
            if c != chrom:
                continue
            pieces = [(x, y) for s, e in pieces for x, y in
                      ((s, min(e, a)), (max(s, b), e)) if x < y]
        result.extend((chrom, s, e) for s, e in pieces)
    return _merge(result)


def _bed(path, lengths, stage='input_check'):
    result = []
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip() or line.startswith(('#', 'track', 'browser')):
            continue
        fields = line.split()
        try:
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if not 0 <= start < end <= lengths[chrom]:
                raise ValueError('coordinates outside reference')
        except (IndexError, KeyError, ValueError) as exc:
            raise LocalFailure('failed', stage, f'Invalid BED line {number}: {path}') from exc
        result.append((chrom, start, end))
    return _merge(result)


def prepare_windows(seeds, fasta, settings, out):
    """Create one finite BED for all candidates and explicit background bins."""
    try:
        lengths = {f[0]: int(f[1]) for f in
                   (line.split() for line in Path(str(fasta) + '.fai').read_text().splitlines()) if f}
        flank = int(settings.get('flank_bp', 100000))
        bin_size = int(settings.get('bin_size', 1000))
        max_bp = int(settings.get('max_window_bp', 25000000))
        if flank < 0 or bin_size <= 0 or max_bp <= 0:
            raise ValueError('invalid local window parameters')
        context = []
        for seed in seeds:
            chrom, start, end = seed['chrom'], int(seed['start']), int(seed['end'])
            if not 0 <= start < end <= lengths[chrom]:
                raise ValueError('seed outside declared reference')
            context.append((chrom, max(0, start - flank), min(lengths[chrom], end + flank)))
    except (OSError, KeyError, IndexError, ValueError) as exc:
        raise LocalFailure('failed', 'input_check', str(exc)) from exc
    context = _merge(context)
    background = _subtract(_bed(settings['background_bed'], lengths), context)
    if not background:
        raise LocalFailure('unsupported', 'input_check', 'Explicit background has no bases outside candidate windows')
    intervals = _merge(context + background)
    total_bp = sum(e - s for _, s, e in intervals)
    if total_bp > max_bp:
        raise LocalFailure('unsupported', 'input_check', f'Local BED budget exceeded: {total_bp} > {max_bp} bp')
    if any(s == 0 and e == lengths[c] for c, s, e in intervals):
        raise LocalFailure('unsupported', 'input_check', 'A local window covers an entire reference contig')
    bins = []
    for label, regions in (('LOCAL_CONTEXT', context), ('LOCAL_BACKGROUND', background)):
        for chrom, start, end in regions:
            bins.extend((chrom, pos, min(pos + bin_size, end), label)
                        for pos in range(start, end, bin_size))
    bins.sort()
    target = Path(out) / 'local.targets.bed'
    target.write_text(''.join('\t'.join(map(str, row)) + '\n' for row in bins))
    return dict(target_bed=str(target), context=context, background=background,
                intervals=intervals, total_bp=total_bp, bin_count=len(bins),
                flank_bp=flank, bin_size=bin_size, max_window_bp=max_bp)


def _sampled_identity(path):
    """At most three 64 KiB blocks, including for BAM/FASTA; never a full scan."""
    path = Path(path).resolve(strict=True)
    stat = path.stat()
    size, block = stat.st_size, 65536
    offsets = sorted(set((0, max(0, size // 2 - block // 2), max(0, size - block))))
    parts = []
    with path.open('rb') as handle:
        for offset in offsets:
            handle.seek(offset)
            parts.append({'offset': offset, 'sha256': hashlib.sha256(handle.read(block)).hexdigest()})
    return dict(path=str(path), bytes=size, mtime_ns=stat.st_mtime_ns,
                method='sampled_first_middle_last_64KiB', blocks=parts)


def _environment(runtime, python):
    env = dict(os.environ)
    env.update(runtime.get('env', {}))
    extra = runtime.get('extra_path', [])
    if not isinstance(extra, list):
        raise LocalFailure('failed', 'input_check', 'runtime.extra_path must be a list')
    env['PATH'] = os.pathsep.join([str(Path(python).absolute().parent), *extra, env.get('PATH', '')])
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['AS_NO_REPO_CHECK'] = '1'
    return env


def _run_tool(python, script, args, policy, cwd, env, stage, timeout=3600):
    command = guarded_command(python, script, args, policy)
    logfile = Path(cwd) / (stage + '.log')
    env = dict(env)
    env.setdefault('MPLCONFIGDIR', str(Path(cwd) / 'matplotlib_cache'))
    Path(env['MPLCONFIGDIR']).mkdir(parents=True, exist_ok=True)
    try:
        rc = execute(command, Path(cwd), env, logfile, timeout=timeout)
    except TimeoutError as exc:
        raise LocalFailure('failed', stage, str(exc)) from exc
    if rc:
        message = f'{stage} exited {rc}; see {logfile}'
        text = logfile.read_text(errors='replace')[-12000:]
        if 'ScanForbidden' in text or 'access prohibited' in text:
            message += '; bounded I/O policy stopped an unsupported read or an exhausted budget'
        raise LocalFailure('failed', stage, message)
    return logfile


_PROBE = '''import base64, hashlib, importlib, importlib.metadata, json, pathlib, shutil, subprocess, sys
result = {'python': sys.version, 'executable': sys.executable, 'modules': {}, 'dependencies': {}}
for name in ('cnvlib', 'skgenome'):
    module = importlib.import_module(name)
    root = pathlib.Path(module.__file__).resolve().parent
    inventory = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in sorted(root.rglob('*.py'))}
    result['modules'][name] = {'root': str(root), 'inventory': inventory,
        'sha256': hashlib.sha256(json.dumps(inventory, sort_keys=True).encode()).hexdigest()}
distribution = importlib.metadata.distribution('CNVkit')
record = {'available': bool(distribution.files), 'checked_python_files': 0, 'changed_files': [], 'missing_files': [], 'import_origins_match': True}
for name in ('cnvlib', 'skgenome'):
    expected = pathlib.Path(distribution.locate_file(name + '/__init__.py')).resolve().parent
    if str(expected) != result['modules'][name]['root']:
        record['import_origins_match'] = False
for entry in distribution.files or []:
    if str(entry).endswith('.py') and entry.hash:
        path = pathlib.Path(distribution.locate_file(entry))
        record['checked_python_files'] += 1
        if not path.is_file():
            record['missing_files'].append(str(entry))
        else:
            actual = base64.urlsafe_b64encode(hashlib.new(entry.hash.mode, path.read_bytes()).digest()).decode().rstrip('=')
            if actual != entry.hash.value: record['changed_files'].append(str(entry))
result['distribution_record'] = record
for name in ('CNVkit', 'numpy', 'scipy', 'pandas', 'pysam', 'biopython', 'pomegranate'):
    try: result['dependencies'][name] = importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError: result['dependencies'][name] = None
path = shutil.which('Rscript')
result['Rscript'] = {'path': path, 'sha256': hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest() if path else None}
if path:
    for label, args in (('version', ['--version']), ('DNAcopy_version', ['--vanilla', '-e', "cat(as.character(utils::packageVersion('DNAcopy')))"])):
        probe = subprocess.run([path] + args, capture_output=True, text=True, timeout=30)
        result['Rscript'][label] = {'returncode': probe.returncode, 'stdout': probe.stdout.strip(), 'stderr': probe.stderr.strip()}
freeze = subprocess.run([sys.executable, '-m', 'pip', 'freeze', '--all'], capture_output=True, text=True, timeout=30)
result['pip_freeze'] = {'returncode': freeze.returncode, 'packages': sorted(freeze.stdout.splitlines()), 'stderr': freeze.stderr.strip()}
pathlib.Path(sys.argv[1]).write_text(json.dumps(result, sort_keys=True, indent=2) + '\\n')
'''


def _probe_cnvkit(runtime, policy, out, env):
    python, script = runtime.get('cnvkit_python'), runtime.get('cnvkit_script')
    if not python or not script or not runtime.get('cnvkit_version'):
        raise LocalFailure('unsupported', 'version_check', 'Explicit CNVkit interpreter, script and version are required')
    for path in (python, script):
        if not Path(path).is_file():
            raise LocalFailure('missing', 'version_check', f'Missing CNVkit executable/script: {path}')
    logfile = _run_tool(python, script, ['version'], policy, out, env, 'cnvkit_version', 90)
    actual = re.findall(r'(?m)^\s*(?:CNVkit\s+)?(\d+\.\d+\.\d+[^\s]*)\s*$',
                        logfile.read_text(errors='replace'))
    if not actual or actual[-1] != str(runtime['cnvkit_version']):
        raise LocalFailure('unsupported', 'version_check',
                           f'Actual CNVkit version does not match {runtime["cnvkit_version"]}: {actual}')
    helper, result = Path(out) / 'probe_cnvkit.py', Path(out) / 'cnvkit_identity.json'
    helper.write_text(_PROBE)
    _run_tool(python, helper, [str(result)], policy, out, env, 'cnvkit_source', 90)
    if not result.is_file():
        raise LocalFailure('missing', 'version_check', 'CNVkit source identity probe produced no output')
    identity = json.loads(result.read_text())
    if identity.get('dependencies', {}).get('CNVkit') != actual[-1]:
        raise LocalFailure('unsupported', 'version_check', 'CNVkit distribution and CLI version differ')
    record = identity.get('distribution_record', {})
    if (not record.get('available') or not record.get('checked_python_files') or
            not record.get('import_origins_match') or record.get('changed_files') or record.get('missing_files')):
        raise LocalFailure('unsupported', 'version_check',
                           'CNVkit imported source cannot be verified against installed distribution RECORD; see cnvkit_identity.json')
    identity.update(version=actual[-1], script=str(Path(script).absolute()), script_sha256=digest(script))
    atomic_json(result, identity)
    return identity


# Native coverage's preliminary sort/read-length checks iterate the BAM from its
# beginning. Replace only these auxiliary checks with finite indexed probes.
# The native depth, reference, bias correction and segmentation implementations
# are called unchanged. This adaptation is explicit in every provenance record.
_CNV_LAUNCHER = '''import itertools, json, pathlib, runpy, sys
from numpy import median
import pysam
from cnvlib import samutil
script, targets = sys.argv[1:3]
def sorted_header(bam, *args, **kwargs):
    with pysam.AlignmentFile(bam, 'rb') as handle:
        return handle.header.to_dict().get('HD', {}).get('SO') == 'coordinate'
def existing_index(bam, *args, **kwargs):
    with pysam.AlignmentFile(bam, 'rb') as handle:
        if not handle.has_index(): raise ValueError('Existing BAM index required; indexing prohibited')
    return bam + '.bai'
def read_length(bam, span=1000, *args, **kwargs):
    own = isinstance(bam, (str, pathlib.Path))
    handle = pysam.AlignmentFile(str(bam), 'rb') if own else bam
    lengths = []
    try:
        for line in pathlib.Path(targets).read_text().splitlines():
            c, s, e = line.split()[:3]
            for read in itertools.islice(handle.fetch(c, int(s), int(e)), max(0, span-len(lengths))):
                if read.query_length: lengths.append(read.query_length)
            if len(lengths) >= span: break
    finally:
        if own: handle.close()
    if not lengths: raise ValueError('No reads in finite CNV coverage windows; cannot estimate auxiliary read length')
    return median(lengths)
samutil.ensure_bam_sorted = sorted_header
samutil.ensure_bam_index = existing_index
samutil.get_read_length = read_length
sys.argv = [script] + sys.argv[3:]
runpy.run_path(script, run_name='__main__')
'''


def _cached_cnv(sample, runtime, config, windows, identity, out, policy, env):
    reference = config['reference']
    reference_mode = config.get('local_cnv', {}).get('reference_mode', 'matched_normal')
    input_keys = ('alignment', 'index') + (('normal_alignment', 'normal_index')
                                         if reference_mode == 'matched_normal' else ())
    inputs = {name: _sampled_identity(sample[name]) for name in input_keys}
    inputs['fasta'] = _sampled_identity(reference['fasta'])
    inputs['fai_sha256'] = digest(str(reference['fasta']) + '.fai')
    parameters = {'windows': {k: v for k, v in windows.items() if k != 'target_bed'},
                  'coverage': 'native_depth', 'processes': 1, 'segment_method': 'cbs',
                  'edge_correction': False, 'local_cns_gap_clipping': True,
                  'reference_mode': reference_mode,
                  'runtime_variant': 'bounded_indexed_auxiliary_probes',
                  'launcher_sha256': hashlib.sha256(_CNV_LAUNCHER.encode()).hexdigest()}
    # Software environment changes invalidate a measurement without changing the
    # paired scientific design. Keep values out of the cache manifest.
    env_spec = {'env': runtime.get('cnvkit_env', runtime.get('env', {})),
                'extra_path': runtime.get('cnvkit_extra_path', runtime.get('extra_path', []))}
    env_sha = hashlib.sha256(json.dumps(env_spec, sort_keys=True).encode()).hexdigest()
    signature = {'sample_id': sample['sample_id'], 'inputs': inputs,
                 'io_contract': {'policy': json.loads(Path(policy).read_text()),
                                 'guard_sha256': digest(Path(__file__).with_name('concordance_policy.py'))},
                 'cnvkit': identity, 'cnv_environment_sha256': env_sha, 'parameters': parameters}
    key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()
    cache_root = Path(config.get('_cnv_cache_dir', Path(out) / 'cnv_cache')) / key
    manifest_path = cache_root / 'complete.json'
    def valid_output(relative, checksum):
        path = cache_root / relative
        return (not Path(relative).is_absolute() and path.resolve().is_relative_to(cache_root.resolve())
                and path.suffix in ('.cnn', '.cnr', '.cns') and path.is_file()
                and path.stat().st_size <= 64 * 1024 ** 2 and digest(path) == checksum)
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text())
            if (manifest['signature'] == json.loads(json.dumps(signature)) and
                    len(manifest['outputs']) == (5 if reference_mode == 'matched_normal' else 4)
                    and manifest['cns'] in manifest['outputs'] and all(
                    valid_output(p, sha)
                    for p, sha in manifest['outputs'].items())):
                return cache_root / manifest['cns'], dict(cache_hit=True, key=key, path=str(cache_root),
                                                         signature=signature)
        except (OSError, KeyError, ValueError):
            pass
    cache_root.mkdir(parents=True, exist_ok=True)
    attempt = cache_root / f'attempt_{1 + len(list(cache_root.glob("attempt_*"))):04d}'
    attempt.mkdir()
    targets = attempt / 'local.targets.bed'
    shutil.copyfile(windows['target_bed'], targets)
    launcher = attempt / 'bounded_cnvkit.py'
    launcher.write_text(_CNV_LAUNCHER)
    empty = attempt / 'tumor.antitargetcoverage.cnn'
    empty.write_text('chromosome\tstart\tend\tgene\tlog2\tdepth\n')
    tumor, normal = attempt / 'tumor.targetcoverage.cnn', attempt / 'normal.targetcoverage.cnn'
    ref, cnr, cns = attempt / 'reference.cnn', attempt / 'tumor.cnr', attempt / 'tumor.cns'
    commands = [('coverage_tumor', ['coverage', str(sample['alignment']), str(targets), '-p', '1', '-o', str(tumor)])]
    expected = [tumor]
    links = [(0, 'alignment', 'index', 'tumor')]
    if reference_mode == 'matched_normal':
        commands.append(('coverage_normal', ['coverage', str(sample['normal_alignment']), str(targets), '-p', '1', '-o', str(normal)]))
        expected.append(normal)
        links.append((1, 'normal_alignment', 'normal_index', 'normal'))
        ref_args = ['reference', str(normal), '-f', str(reference['fasta']), '--no-edge', '-o', str(ref)]
    else:
        ref_args = ['reference', '-t', str(targets), '-f', str(reference['fasta']), '--no-edge', '-o', str(ref)]
    commands.extend([('reference', ref_args),
                     ('fix', ['fix', str(tumor), str(empty), str(ref), '--no-edge', '-o', str(cnr)]),
                     ('segment', ['segment', str(cnr), '-m', 'cbs', '-p', '1', '-o', str(cns)])])
    expected.extend([ref, cnr, cns])
    atomic_json(attempt / 'provenance.json', dict(created=now(), signature=signature,
                commands=commands, note='Finite local bins; empty antitarget array; native WGS-style no-edge correction'))
    # Native bedcov locates indexes by filename; use private links when the
    # supplied index has a nonstandard name. No index is created or touched.
    for position, bam_key, index_key, stem in links:
        link = attempt / (stem + '.bam')
        link.symlink_to(Path(sample[bam_key]).resolve())
        suffix = '.csi' if str(sample[index_key]).lower().endswith('.csi') else '.bai'
        Path(str(link) + suffix).symlink_to(Path(sample[index_key]).resolve())
        commands[position][1][1] = str(link)
    atomic_json(attempt / 'commands.json', commands)
    for (stage, args), result in zip(commands, expected):
        _run_tool(runtime['cnvkit_python'], launcher,
                  [runtime['cnvkit_script'], str(targets), *args], policy, attempt, env,
                  stage, config.get('local_cnv', {}).get('timeout_seconds', 3600))
        if not result.is_file() or not result.stat().st_size:
            raise LocalFailure('missing', stage, f'Expected CNVkit result is absent/empty: {result}')
    outputs = {str(p.relative_to(cache_root)): digest(p) for p in expected}
    atomic_json(manifest_path, dict(signature=signature, outputs=outputs,
                                   cns=str(cns.relative_to(cache_root)), completed=now()))
    return cns, dict(cache_hit=False, key=key, path=str(cache_root), signature=signature)


def cns_to_local_bed(cns, intervals, destination):
    """Use every segment, retaining CN and clipping only unobserved inter-window gaps."""
    rows = []
    with Path(cns).open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not {'chromosome', 'start', 'end', 'log2'}.issubset(reader.fieldnames or []):
            raise LocalFailure('failed', 'result_reading', 'CNS lacks chromosome/start/end/log2 columns')
        try:
            for row in reader:
                start, end, log2 = int(row['start']), int(row['end']), float(row['log2'])
                if start < 0 or end <= start or not math.isfinite(log2):
                    raise ValueError('invalid CNS value')
                cn = 2 ** (log2 + 1)
                if not math.isfinite(cn) or cn <= 0:
                    raise ValueError('invalid CNS copy number')
                for chrom, a, b in intervals:
                    if chrom == row['chromosome'] and max(a, start) < min(b, end):
                        rows.append((chrom, max(a, start), min(b, end), cn))
        except (KeyError, ValueError, OverflowError) as exc:
            raise LocalFailure('failed', 'result_reading', f'Malformed CNS: {exc}') from exc
    if not rows:
        raise LocalFailure('missing', 'result_reading', 'CNS has no valid segments in the finite windows')
    Path(destination).write_text(''.join(f'{c}\t{s}\t{e}\t{cn:.12g}\n' for c, s, e, cn in rows))
    return rows


def retention_fraction(seed, intervals):
    pieces = [(c, max(int(seed['start']), s), min(int(seed['end']), e)) for c, s, e in intervals
              if c == seed['chrom'] and max(int(seed['start']), s) < min(int(seed['end']), e)]
    return sum(e - s for _, s, e in _merge(pieces)) / (int(seed['end']) - int(seed['start']))


def calibrate_source_retention(sample, runtime, rows, out):
    """Check an existing, originally recorded same-input positive control.

    The provenance must come from the original upstream run. This function does
    not manufacture historical input identity from today's files. Its endpoint
    is only reproduction of supplied source candidates, not global CN accuracy.
    """
    if runtime['id'] != sample.get('source_version'):
        for row in rows:
            row.update(calibration_role='candidate_observation', calibration_status='unverified',
                       calibration_reason='Candidate retention is an observation; calibration is assessed on the independent source-version replay')
        return rows
    for row in rows:
        row['calibration_role'] = 'source_control'
    required = ('baseline_cns', 'source_seed_bed', 'baseline_provenance')
    evidence = dict(scope='selected_source_candidate_retention_only', status='unverified',
                    required_source_files=list(required),
                    required_original_metadata={
                        'provenance_kind': 'original_upstream_run',
                        'cnvkit_version': 'actual original CNVkit version',
                        'reference_build': 'original declared build',
                        'reference_mode': rows[0].get('reference_mode', 'matched_normal'),
                        'inputs': {'alignment': 'original sampled identity from _sampled_identity',
                                   'normal_alignment': 'original sampled identity from _sampled_identity'},
                        'outputs': {'baseline_cns_sha256': 'original CNS SHA256',
                                    'source_seed_bed_sha256': 'original seed BED SHA256'}},
                    warning='Do not create original-run metadata retrospectively from current files. Sampled identity is not a full BAM checksum.')
    if rows[0].get('reference_mode') == 'flat':
        evidence['required_original_metadata']['inputs'].pop('normal_alignment')
    reason = 'Original-run input/source provenance is absent; a version declaration alone cannot calibrate local CNV'
    try:
        if not all(sample.get(k) and Path(sample[k]).is_file() for k in required):
            raise ValueError(reason)
        metadata = json.loads(Path(sample['baseline_provenance']).read_text())
        if metadata.get('provenance_kind') != 'original_upstream_run':
            raise ValueError('Calibration metadata is not original upstream-run provenance')
        if not (str(sample.get('source_cnvkit_version')) == str(runtime['cnvkit_version']) ==
                str(metadata.get('cnvkit_version'))):
            raise ValueError('Source CNVkit version is absent or differs from the local control runtime')
        if metadata.get('reference_build') != sample.get('reference_build'):
            raise ValueError('Original and local reference builds differ')
        reference_mode = rows[0].get('reference_mode', 'matched_normal')
        if metadata.get('reference_mode', 'matched_normal') != reference_mode:
            raise ValueError('Original and local reference strategies differ; flat is not a matched normal')
        input_keys = ('alignment', 'normal_alignment') if reference_mode == 'matched_normal' else ('alignment',)
        for key in input_keys:
            if metadata.get('inputs', {}).get(key) != _sampled_identity(sample[key]):
                raise ValueError(f'Original same-input sampled identity is missing or differs: {key}')
        for key, field in (('baseline_cns', 'baseline_cns_sha256'), ('source_seed_bed', 'source_seed_bed_sha256')):
            if metadata.get('outputs', {}).get(field) != digest(sample[key]):
                raise ValueError(f'Original output checksum differs: {key}')
        original_seeds = []
        for line in Path(sample['source_seed_bed']).read_text().splitlines():
            if line.strip() and not line.startswith(('#', 'track', 'browser')):
                f = line.split()
                if int(f[1]) < 0 or int(f[2]) <= int(f[1]):
                    raise ValueError('Invalid source seed BED')
                original_seeds.append((f[0], int(f[1]), int(f[2])))
        source_segments = cns_to_local_bed(sample['baseline_cns'],
                [(r['seed']['chrom'], int(r['seed']['start']), int(r['seed']['end'])) for r in rows],
                Path(out) / 'calibration_source_cnv.bed')
        source_intervals = [(c, s, e) for c, s, e, cn in source_segments]
        evidence.update(status='evaluated', provenance_sha256=digest(sample['baseline_provenance']),
                        baseline_cns_sha256=digest(sample['baseline_cns']),
                        source_seed_bed_sha256=digest(sample['source_seed_bed']), results=[])
        for row in rows:
            fraction = retention_fraction(row['seed'], original_seeds)
            cn_coverage = retention_fraction(row['seed'], source_intervals)
            supported = min(fraction, cn_coverage) >= row['min_seed_overlap']
            status = ('passed' if row['status'] == 'retained' else 'failed') if supported else 'unverified'
            detail = dict(seed_id=row['seed_id'], source_candidate_fraction=fraction,
                          source_cns_coverage_fraction=cn_coverage, local_retention_fraction=row['retention_fraction'],
                          status=status)
            evidence['results'].append(detail)
            row.update(calibration_status=status, calibration_scope=evidence['scope'],
                       calibration_reason=('Same-version source candidate was reproduced locally' if status == 'passed'
                                           else 'Same-version source candidate was not reproduced locally' if status == 'failed'
                                           else 'Supplied source candidate/CNS does not sufficiently cover this requested seed'))
    except (OSError, ValueError, KeyError, IndexError, LocalFailure) as exc:
        reason = str(exc)
        evidence['reason'] = reason
        for row in rows:
            row.update(calibration_status='unverified', calibration_reason=reason)
    atomic_json(Path(out) / 'calibration_evidence.json', evidence)
    for row in rows:
        row['calibration_evidence'] = str(Path(out) / 'calibration_evidence.json')
    return rows


def _previous_filter_bed(manifest, identity):
    """Validate a prior successful controlled filter before reusing its small BED."""
    if not manifest.is_file():
        return None
    try:
        saved = json.loads(manifest.read_text())
        if saved['identity'] != identity:
            return None
        attempt = saved['producer_attempt']
        if not re.fullmatch(r'attempt_[0-9]+', attempt):
            return None
        directory = manifest.parent / attempt
        bed = directory / 'local_amplified.bed'
        if (directory.is_symlink() or bed.is_symlink() or not bed.is_file()
                or bed.stat().st_size > 64 * 1024 * 1024 or digest(bed) != saved['bed_sha256']):
            return None
        for name in ('upstream_provenance.json', 'seed_processing.log'):
            artifact = directory / name
            if artifact.is_symlink() or not artifact.is_file() or digest(artifact) != saved['producer_sha256'][name]:
                return None
        return bed
    except (OSError, KeyError, TypeError, ValueError):
        return None


def run_local_upstream(sample, runtime, entry, config, out):
    """Run bounded upstream once per sample/version; return candidate-level records."""
    out = Path(out).absolute()
    out.mkdir(parents=True, exist_ok=True)
    settings = dict(config.get('local_cnv', {}))
    seeds = sample.get('seeds', [])
    target_version = runtime['id']
    source_version = sample.get('source_version')
    mapped = load_builtin_map()['versions'].get(target_version, {})
    base = dict(sample_id=sample['sample_id'], version=target_version, source_version=source_version,
                cnvkit_version=runtime.get('cnvkit_version'), endpoint='local_candidate_retention',
                calibration_role='source_control' if target_version == source_version else 'candidate_observation',
                calibration_status='unverified',
                calibration_reason='No matched same-CNVkit, same-input whole-genome calibration evidence has been evaluated; local bins, centering and segmentation context differ',
                matching_rule='union_overlap_bp / original_seed_bp >= configured min_seed_overlap',
                threshold_origin='user-configurable concordance threshold; not an official AA/CNVkit threshold',
                native_recipe_cnvkit=mapped.get('cnvkit_version'),
                native_recipe_match=runtime.get('cnvkit_version') == mapped.get('cnvkit_version'),
                comparison_design='native_recipe_cnvkit' if runtime.get('cnvkit_version') == mapped.get('cnvkit_version')
                                  else 'controlled_cnvkit',
                runtime_variant='bounded_indexed_auxiliary_probes',
                dependency_stack_kind='recorded_runtime_not_reconstructed_container',
                reference_mode=settings.get('reference_mode', 'matched_normal'),
                retention_fraction=None, amplification_confirmed=False, comparison_eligible=False)
    def records(status, stage, reason, **fields):
        rows = [dict(base, seed_id=s.get('seed_id', '__no_seeds__'), seed=s, status=status,
                     stage=stage, reason=reason, **fields) for s in (seeds or [{}])]
        atomic_json(out / 'upstream_records.json', rows)
        return rows
    if not seeds:
        return records('no_seeds', 'input_check', 'No existing candidate regions were supplied')
    try:
        threshold = float(settings.get('min_seed_overlap', 1.0))
        if not math.isfinite(threshold) or not 0 < threshold <= 1:
            raise LocalFailure('failed', 'input_check', 'min_seed_overlap must be >0 and <=1')
        base['min_seed_overlap'] = threshold
        mode = settings.get('reference_mode', 'matched_normal')
        if mode not in ('matched_normal', 'flat'):
            raise LocalFailure('unsupported', 'input_check', 'reference_mode must be matched_normal or explicitly flat')
        if not settings.get('background_bed'):
            raise LocalFailure('unsupported', 'input_check', 'An explicit reliable background BED is required for local CNV')
        if mode == 'matched_normal' and (not sample.get('normal_alignment') or not sample.get('normal_index')):
            raise LocalFailure('unsupported', 'input_check', 'Matched-normal mode requires a paired indexed normal BAM; flat is never selected silently')
        alignment_keys = ('alignment', 'normal_alignment') if mode == 'matched_normal' else ('alignment',)
        if any(not str(sample.get(k, '')).lower().endswith('.bam') for k in alignment_keys):
            raise LocalFailure('unsupported', 'input_check', 'Local CNV supports indexed BAM only; no CRAM conversion or whole-genome scan')
        required_keys = ('alignment', 'index') + (('normal_alignment', 'normal_index') if mode == 'matched_normal' else ())
        required = [sample.get(k) for k in required_keys]
        required += [settings['background_bed'], config['reference']['fasta'], str(config['reference']['fasta']) + '.fai']
        if any(not p or not Path(p).is_file() for p in required):
            raise LocalFailure('missing', 'input_check', 'Required BAM/index, background BED or indexed FASTA is missing')
        policy = config.get('_io_policy_path')
        if not policy or not Path(policy).is_file():
            raise LocalFailure('unsupported', 'input_check', 'Bounded indexed-BAM policy is required')
        windows = prepare_windows(seeds, config['reference']['fasta'], settings, out)
        cnv_environment = dict(runtime, env=runtime.get('cnvkit_env', runtime.get('env', {})),
                               extra_path=runtime.get('cnvkit_extra_path', runtime.get('extra_path', [])))
        env = _environment(cnv_environment, runtime.get('cnvkit_python', runtime['aa_python']))
        identity = _probe_cnvkit(runtime, policy, out, env)
        cns, cache = _cached_cnv(sample, runtime, config, windows, identity, out, policy, env)
        local_bed = out / 'all_local_cnv_segments.bed'
        segments = cns_to_local_bed(cns, windows['intervals'], local_bed)
        helper = Path(runtime['aa_script']).absolute().parent / 'amplified_intervals.py'
        if not helper.is_file():
            raise LocalFailure('unsupported', 'seed_processing', 'Locked AA source lacks amplified_intervals.py')
        aa_env = _environment(runtime, runtime['aa_python'])
        private_ref = local_reference(config['reference']['aa_data_repo'], out / 'private_AA_DATA_REPO')
        aa_env['AA_DATA_REPO'] = str(private_ref)
        gain = float(settings.get('cngain', config.get('full_suite', {}).get('cngain', 4.5)))
        min_size = int(settings.get('cnsize_min', config.get('full_suite', {}).get('cnsize_min', 50000)))
        if not math.isfinite(gain) or gain <= 0 or min_size < 0:
            raise LocalFailure('failed', 'input_check', 'Invalid AA seed gain/minimum size')
        prefix = out / 'local_amplified'
        args = ['--ref', config['reference']['build'], '--bed', str(local_bed),
                '--bam', str(sample['alignment']), '--gain', str(gain),
                '--cnsize_min', str(min_size), '--out', str(prefix)]
        if '--no_cstats' in helper.read_text():
            args.append('--no_cstats')
        provenance = dict(created=now(), version=target_version, source_version=source_version, cnvkit=identity,
                          reference_mode=mode,
                          windows=windows, cache=cache, all_local_segment_count=len(segments),
                          amplified_intervals=dict(path=str(helper), sha256=digest(helper), args=args),
                          io_policy=str(Path(policy).absolute()), calibration_status='unverified',
                          calibration_reason=base['calibration_reason'], baseline_inputs={
                              name: _sampled_identity(sample[name]) for name in ('baseline_cns', 'source_seed_bed')
                              if sample.get(name) and Path(sample[name]).is_file()},
                          limitations=['Local CNV is not whole-genome CNVkit', 'Finite background is user-selected and is not certified diploid',
                                       'Native CNS segments are clipped at unobserved gaps between local windows',
                                       'AA seed filtering retains its bounded indexed background estimation',
                                       'Candidate retention does not establish biological amplification or pipeline compatibility'])
        output = prefix.with_suffix('.bed')
        filter_manifest = out.parent / 'filter_cache.json'
        filter_identity, previous_bed = None, None
        if (config.get('_controlled_cnv_role') and config.get('_fixed_source_aa_identity')
                and config.get('_reference_identity')):
            filter_inputs = dict(cnv_measurement=cache['key'], local_segments_sha256=digest(local_bed),
                fixed_aa=config['_fixed_source_aa_identity'], build=config['reference']['build'],
                reference=config['_reference_identity'], helper_sha256=digest(helper),
                env=runtime.get('env', {}), extra_path=runtime.get('extra_path', []), gain=gain,
                min_size=min_size, no_cstats='--no_cstats' in args,
                candidate_ids=[s['seed_id'] for s in seeds])
            filter_identity = hashlib.sha256(json.dumps(filter_inputs, sort_keys=True).encode()).hexdigest()
            previous_bed = _previous_filter_bed(filter_manifest, filter_identity)
            provenance['seed_filter_cache'] = dict(identity=filter_identity, hit=previous_bed is not None,
                                                    producer_bed=str(previous_bed) if previous_bed else None)
        atomic_json(out / 'upstream_provenance.json', provenance)
        if previous_bed:
            shutil.copyfile(previous_bed, output)
            if digest(output) != digest(previous_bed):
                raise LocalFailure('failed', 'seed_processing', 'Seed filter cache changed during reuse')
            atomic_json(out / 'seed_processing.log', dict(cache_hit=True, producer_bed=str(previous_bed)))
        else:
            _run_tool(runtime['aa_python'], helper, args, policy, out, aa_env,
                      'seed_processing', settings.get('timeout_seconds', 3600))
            filter_log = (out / 'seed_processing.log').read_text(errors='replace')
            if 'Traceback (most recent call last)' in filter_log or 'zero division error' in filter_log:
                raise LocalFailure('failed', 'seed_processing', 'AA seed filter logged a numerical/processing error; empty output is not a valid negative')
        if not output.is_file():
            raise LocalFailure('missing', 'result_reading', 'AA seed filtering exited successfully but its BED is missing')
        lengths = {f[0]: int(f[1]) for f in
                   (line.split() for line in Path(str(config['reference']['fasta']) + '.fai').read_text().splitlines()) if f}
        intervals = _bed(output, lengths, 'result_reading')
        if filter_identity and not previous_bed:
            atomic_json(filter_manifest, dict(identity=filter_identity, producer_attempt=out.name,
                bed_sha256=digest(output), producer_sha256={n: digest(out / n)
                    for n in ('upstream_provenance.json', 'seed_processing.log')}))
        rows = []
        for seed in seeds:
            fraction = retention_fraction(seed, intervals)
            rows.append(dict(base, seed_id=seed['seed_id'], seed=seed,
                             status='retained' if fraction >= threshold else 'not_retained',
                             stage='seed_processing', reason='Candidate retention after bounded local CNV and native AA seed filtering',
                             retention_fraction=fraction, retained_intervals=intervals,
                             cnv_cache_hit=cache['cache_hit'], cnv_cache_key=cache['key'],
                             provenance=str(out / 'upstream_provenance.json')))
            rows[-1]['comparison_eligible'] = True
        if config.get('_controlled_cnv_role') == 'comparison':
            for row in rows:
                row.update(calibration_role='candidate_observation', calibration_status='unverified',
                           calibration_reason='Controlled CNVkit version compared to the sample-matched local standard')
        else:
            calibrate_source_retention(sample, runtime, rows, out)
        atomic_json(out / 'upstream_records.json', rows)
        return rows
    except LocalFailure as exc:
        return records(exc.status, exc.stage, str(exc))
    except (OSError, ValueError, KeyError, TypeError, TimeoutError) as exc:
        return records('failed', 'input_check' if not (out / 'upstream_provenance.json').exists()
                       else 'seed_processing', f'{type(exc).__name__}: {exc}')
