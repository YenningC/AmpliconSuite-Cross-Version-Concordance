"""Read a finished global run and isolate AA failures caused by the unique-base limit."""
import hashlib
import json
from pathlib import Path
import re

from .common import digest
from .concordance_global_report import global_loci, seeds_by_locus
from .concordance_global_seeds import FINAL_BED, _native_seeds

_REPORT_INPUTS = frozenset(('config.resolved.json', 'downstream_records.json',
                            'upstream_records.json', 'local_cnv_comparator_records.json'))
_RECOVERY_KEYS = {'recover_aa_from', 'aa_rescue_unique_bp'}
_SAMPLE_OUTPUTS = {'_policy_path', 'global_seed_sets', 'upstream_seed_ids',
                   'selected_seed_ids', 'global_executed_seed_ids',
                   'global_source_seed_calibration', 'source_bridge_status',
                   'global_version_statuses'}


def _read_json(path):
    return json.loads(Path(path).read_text())


def _within(value, root):
    path = Path(value)
    if not path.is_absolute() or path.resolve() != path or not path.is_relative_to(root):
        raise ValueError(f'Previous artifact path is outside its run: {path}')
    return path


def _match_file(path, checksum):
    if not path.is_file() or path.is_symlink() or digest(path) != checksum:
        raise ValueError(f'Previous artifact is missing or changed: {path}')


def _read_report_inputs(root):
    checksums = _read_json(root / 'report_inputs.sha256.json')
    if not isinstance(checksums, dict) or set(checksums) != _REPORT_INPUTS:
        raise ValueError('AA recovery requires a finished global run with a complete report checksum manifest')
    documents = {}
    for name, checksum in checksums.items():
        data = (root / name).read_bytes()
        if hashlib.sha256(data).hexdigest() != checksum:
            raise ValueError(f'Previous report input changed: {name}')
        documents[name] = json.loads(data)
    return checksums, documents


def _verified_seeds(root, version, sample, call, lengths):
    cns = Path(call['cns'])
    receipt = _read_json(_within(call['cnvkit_provenance'], root))
    if str(cns) != (receipt.get('cns') or receipt.get('file')):
        raise ValueError(f'Previous CNS receipt does not match the sealed report for {version}')
    _match_file(cns, receipt.get('sha256'))
    attempt = _within(call['seed_attempt'], root)
    job = root / 'global_upstream' / sample['sample_id'] / version
    if attempt.parent != job:
        raise ValueError(f'Previous seed attempt is outside its Suite job for {version}')
    provenance_path = attempt / 'upstream_provenance.json'
    provenance = _read_json(provenance_path)
    if (provenance['version'] != version or provenance['cns'] != str(cns)
            or provenance['io_policy'] != sample['_policy_path']):
        raise ValueError(f'Previous seed provenance does not match its input for {version}')
    bed = attempt / FINAL_BED
    _match_file(bed, provenance['seed_bed_sha256'])
    if _native_seeds(bed, lengths) != call['seeds']:
        raise ValueError(f'Previous native seed BED differs from the sealed report for {version}')
    return dict(sha256=provenance['seed_bed_sha256'], cns_sha256=receipt['sha256'])


def _verify_original_seed_input(attempt, native):
    from .concordance_components import _seed_bed_text

    bed = _within(attempt / 'selected_seeds.bed', attempt)
    if bed.read_text() != _seed_bed_text(native):
        raise ValueError(f'Previous AA seed input differs from the frozen panel: {attempt}')


def _verify_completed_job(version, job, records, sample_id, native):
    from .concordance import project_record

    cache = _read_json(job / 'cache.json')
    joint = cache['record']
    attempt = _within(joint['run_dir'], job)
    files = cache['files']
    if (joint['status'] not in ('ok', 'not_detected') or joint['sample_id'] != sample_id
            or joint['version'] != version or not joint.get('comparison_eligible')
            or any(r.get('run_dir') != str(attempt) for r in records)
            or str(attempt / 'stages.json') not in files or str(attempt / 'selected_seeds.bed') not in files
            or not any(Path(name).name.endswith('_summary.txt') for name in files)):
        raise ValueError(f'Previous completion cache is inconsistent for {version}')
    for name, checksum in files.items():
        _match_file(_within(name, attempt), checksum)
    _verify_original_seed_input(attempt, native)
    stages = _read_json(attempt / 'stages.json')
    if (stages['aa']['status'] != 'complete'
            or stages['ac']['status'] not in ('complete', 'skipped_zero_amplicons')):
        raise ValueError(f'Previous completed result lacks successful AA/AC stages for {version}')
    fields = ('status', 'labels', 'features', 'intervals', 'amplicons', 'reconstruction_intervals',
              'supported_labels', 'coverage_status')
    for row in records:
        projected = project_record(joint, row['seed'])
        if any(row.get(field) != projected.get(field) for field in fields):
            raise ValueError(f'Previous completed record differs from the sealed report for {version}')
    return joint


def _budget_failure(version, job, records, native, new_limit):
    """Use the recorded guard rejection, plus the matching failed AA attempt."""
    reasons = {r.get('reason', '') for r in records
               if 'Unique genomic coverage budget exceeded' in r.get('reason', '')}
    if not reasons:
        return None
    if len(reasons) != 1 or any(r['status'] != 'failed' or r['stage'] != 'aa' for r in records):
        raise ValueError(f'Previous AA limit failure is inconsistent for {version}')
    match = re.fullmatch(r'Unique genomic coverage budget exceeded(?:; genome tiling is prohibited)?; see (bam_access\.\d+\.json)',
                         reasons.pop())
    attempts = {r.get('run_dir') for r in records}
    if not match or len(attempts) != 1 or None in attempts:
        raise ValueError(f'Previous AA limit failure lacks a unique audit for {version}')
    attempt = _within(attempts.pop(), job)
    stages_path = _within(attempt / 'stages.json', attempt)
    stages = _read_json(stages_path)
    audit_path = _within(attempt / match[1], attempt)
    audit = _read_json(audit_path)
    rejected = audit.get('first_rejected_unique_request', {})
    old_limit = rejected.get('limit_bp', 0)
    if (stages['aa']['status'] != 'failed' or stages['aa'].get('returncode') in (None, 0)
            or 'ac' in stages or stages['aa'].get('reason') != records[0]['reason']
            or not audit.get('blocked', []) or not audit['blocked'][0].startswith('Unique genomic coverage budget exceeded')
            or audit.get('effective_max_unique_bp') != old_limit or old_limit <= 0
            or not audit['unique_requested_bp'] <= old_limit < rejected['projected_unique_bp']):
        raise ValueError(f'Previous BAM access audit does not confirm a unique-base AA rejection for {version}')
    _verify_original_seed_input(attempt, native)
    if new_limit < rejected['projected_unique_bp']:
        raise ValueError(f'aa_rescue_unique_bp must permit the previously rejected request for {version}')
    return dict(previous_audit=str(audit_path), previous_audit_sha256=digest(audit_path),
                previous_limit_bp=old_limit, first_rejected_request=rejected,
                solver=stages['aa'].get('solver_selection'))


def verify_recovery_source(config, out):
    """Validate the immutable previous run and all samples before restarting any AA."""
    root, out = Path(config['recover_aa_from']).resolve(), Path(out).resolve()
    if root == out or root.is_relative_to(out) or out.is_relative_to(root):
        raise ValueError('AA recovery needs a new output directory separate from the original')
    checksums, inputs = _read_report_inputs(root)
    prior = inputs['config.resolved.json']
    if prior['mode'] != 'with_bam' or prior.get('cnv_scope') != 'global' or prior.get('recover_aa_from'):
        raise ValueError('AA recovery requires the finished original global with_bam run')
    changed = set(k for k in set(prior) | set(config) if prior.get(k) != config.get(k)) - _RECOVERY_KEYS - {'samples'}
    if changed:
        raise ValueError('AA recovery changed original settings: ' + ', '.join(sorted(changed)))
    for name in ('input_identity.json', 'versions.lock.json'):
        if _read_json(root / name) != _read_json(out / name):
            raise ValueError(f'AA recovery changed original inputs or version lock: {name}')
    new_samples = config['samples']
    if len(new_samples) != len(prior['samples']):
        raise ValueError('AA recovery must preserve the original sample panel')
    recovered = {}
    for old, new in zip(prior['samples'], new_samples):
        unchanged = {k: v for k, v in old.items() if k not in _SAMPLE_OUTPUTS}
        if unchanged != {k: v for k, v in new.items() if k not in _SAMPLE_OUTPUTS}:
            raise ValueError(f'AA recovery changed the input for {old["sample_id"]}')
        sid, sets = old['sample_id'], old['global_seed_sets']
        previous_policy = _within(old['_policy_path'], root)
        if _read_json(previous_policy) != _read_json(Path(new['_policy_path'])):
            raise ValueError(f'AA recovery changed the original indexed BAM policy for {sid}')
        versions = old['global_version_statuses']
        discovery = root / 'global_upstream' / sid / 'discovery.json'
        sampling = _read_json(root / 'sampling' / (sid + '.json'))
        loci = global_loci(sets)
        selected = sampling['selected']
        if (sets != _read_json(discovery) or list(versions) != list(sets)
                or sampling['all_seeds'] != loci or [s['seed_id'] for s in loci] != old['upstream_seed_ids']
                or [s['seed_id'] for s in selected] != old['selected_seed_ids'] or any(s not in loci for s in selected)
                or sampling['native_selected_ids'] != old['global_executed_seed_ids']):
            raise ValueError(f'AA recovery requires the original frozen seed panel for {sid}')
        by_version = {v: [] for v in versions}
        for row in inputs['downstream_records.json']:
            if row['sample_id'] == sid:
                by_version[row['version']].append(row)
        successes, retries, seeds = {}, {}, {}
        selected_ids = set(old['selected_seed_ids'])
        panel = selected_ids or {'__no_seeds__'}
        for version, status in versions.items():
            records = by_version[version]
            if len(records) != len(panel) or {r['seed_id'] for r in records} != panel:
                raise ValueError(f'Previous result rows do not cover the frozen panel for {sid} / {version}')
            call = sets[version]
            native = []
            if call['status'] == 'ok':
                seeds[version] = _verified_seeds(root, version, old, call, config['_reference_lengths'])
                members = seeds_by_locus(loci, call['seeds'])
                seed_ids = {s['seed_id'] for locus in selected_ids for s in members.get(locus, [])}
                native = [s for s in call['seeds'] if s['seed_id'] in seed_ids]
                if [s['seed_id'] for s in native] != old['global_executed_seed_ids'].get(version, []):
                    raise ValueError(f'Previous selected seed rows disagree with the frozen panel for {version}')
            job = root / 'downstream_jobs' / sid / version
            if status['status'] in ('ok', 'not_detected') and native:
                successes[version] = _verify_completed_job(version, job, records, sid, native)
            retry = _budget_failure(version, job, records, native, config['aa_rescue_unique_bp'])
            if retry:
                if call['status'] != 'ok' or not old['global_executed_seed_ids'].get(version):
                    raise ValueError(f'Previous AA failure lacks input seeds for {sid} / {version}')
                retries[version] = retry
        recovered[sid] = dict(previous_sample=old, sampling=sampling, previous_downstream=by_version,
                              completed=successes, retries=retries, checksums=checksums, seeds=seeds)
    return recovered
