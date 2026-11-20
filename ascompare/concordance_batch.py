"""Import an explicit, local sample table into the existing serial workflow.

This adapter reads only its TSV/JSON and known result-directory names. BAM and
component checks remain the job of ``asconcordance check`` and ``run``.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile

from .common import digest
from .inputs import safe_id


REQUIRED_COLUMNS = {'sample_id', 'source_version', 'reference_build'}
SAMPLE_COLUMNS = REQUIRED_COLUMNS | {
    'patient_id', 'alignment', 'index', 'normal_alignment', 'normal_index',
    'source_cnvkit_version', 'baseline_cns', 'source_seed_bed',
    'baseline_provenance',
}
RESULT_COLUMNS = {
    'result_dir', 'run_dir', 'aa_dir', 'ac_dir', 'stages_json', 'seeds_bed',
    'aa_version', 'ac_version',
}
ALLOWED_COLUMNS = SAMPLE_COLUMNS | RESULT_COLUMNS | {'result_version'}


def _rows(path):
    """Read strict UTF-8 TSV, including a spreadsheet UTF-8 byte-order mark."""
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.reader(handle, delimiter='\t', strict=True)
        try:
            header = [value.strip() for value in next(reader)]
        except StopIteration:
            raise ValueError('Sample TSV is empty') from None
        if not all(header) or len(header) != len(set(header)):
            raise ValueError('Sample TSV has an empty or duplicate column name')
        unknown = set(header) - ALLOWED_COLUMNS
        if unknown:
            raise ValueError('Unknown sample TSV columns: ' + ', '.join(sorted(unknown)))
        missing = REQUIRED_COLUMNS - set(header)
        if missing:
            raise ValueError('Sample TSV needs columns: ' + ', '.join(sorted(missing)))
        count = 0
        for values in reader:
            if not values or not any(value.strip() for value in values):
                continue
            line = reader.line_num
            if len(values) != len(header):
                raise ValueError(f'Sample TSV line {line}: expected {len(header)} tab-separated fields, got {len(values)}')
            row = {key: value.strip() for key, value in zip(header, values) if value.strip()}
            missing = REQUIRED_COLUMNS - set(row)
            if missing:
                raise ValueError(f'Sample TSV line {line}: empty required fields: ' + ', '.join(sorted(missing)))
            count += 1
            yield line, row
        if not count:
            raise ValueError('Sample TSV has no sample rows')


def import_samples(samples_path, config_path, output_path):
    """Validate and save a new manifest without starting scientific components.

    A row describes one version's existing output. ``result_version`` defaults
    to ``source_version``. Repeated sample IDs describe extra existing versions;
    their nonempty sample metadata must agree, and each needs exactly one source
    result. Sample order follows first appearance in the TSV.
    """
    # Imported lazily to keep the command-line entry point free of import cycles.
    from .concordance import _resolve, cnv_scope, load_manifest

    samples_path = Path(samples_path).expanduser().resolve(strict=True)
    config_path = Path(config_path).expanduser().resolve(strict=True)
    output_path = Path(output_path).expanduser().absolute()
    if output_path.exists() or output_path.is_symlink():
        raise ValueError('Output manifest already exists; choose a new --out path')
    config = json.loads(config_path.read_text(encoding='utf-8'))
    if not isinstance(config, dict):
        raise ValueError('Shared config must be a JSON object')
    if config.get('samples') not in (None, []):
        raise ValueError('Shared config must omit samples or use an empty samples list; the TSV supplies all samples')
    if 'batch_import' in config:
        raise ValueError('Shared config must not contain batch_import from an earlier import')
    config = _resolve(config, config_path.parent)
    if config.get('mode') not in ('with_bam', 'results_only'):
        raise ValueError('Shared config mode must be with_bam or results_only')

    samples, results = {}, {}
    for line, row in _rows(samples_path):
        for key in ('sample_id', 'source_version', 'patient_id', 'result_version'):
            if key in row:
                try:
                    safe_id(row[key])
                except (TypeError, ValueError) as exc:
                    raise ValueError(f'Sample TSV line {line}, {key}: {exc}') from exc
        row = _resolve(row, samples_path.parent)
        sid, source = row['sample_id'], row['source_version']
        version = row.get('result_version', source)
        sample = samples.setdefault(sid, {})
        for key in sorted(SAMPLE_COLUMNS):
            if key not in row:
                continue
            if key in sample and sample[key] != row[key]:
                raise ValueError(f'Sample TSV line {line}: conflicting {key} for sample {sid}')
            sample[key] = row[key]
        by_version = results.setdefault(sid, {})
        if version in by_version:
            raise ValueError(f'Sample TSV line {line}: duplicate result_version {version} for sample {sid}')
        spec = {key: row[key] for key in sorted(RESULT_COLUMNS) if key in row}
        if not (spec.get('result_dir') or spec.get('run_dir') or
                (spec.get('aa_dir') and spec.get('ac_dir'))):
            raise ValueError(f'Sample TSV line {line}: provide result_dir, run_dir, or both aa_dir and ac_dir')
        by_version[version] = spec

    for sid, sample in samples.items():
        source = sample['source_version']
        if source not in results[sid]:
            raise ValueError(f'Sample {sid} has no result row for its source_version {source}')
        sample['source_result'] = results[sid][source]
        sample['comparison_results'] = [dict(version=version, **spec)
            for version, spec in results[sid].items() if version != source]
        if sample['source_result'].get('seeds_bed'):
            sample['seeds_bed'] = sample['source_result']['seeds_bed']
        if config['mode'] == 'with_bam':
            required = {'alignment', 'index'}
            if cnv_scope(config) == 'local':
                required.add('seeds_bed')
            missing = required - set(sample)
            if missing:
                raise ValueError(f'Sample {sid}: with_bam needs ' + ', '.join(sorted(missing)))
            if bool(sample.get('normal_alignment')) != bool(sample.get('normal_index')):
                raise ValueError(f'Sample {sid}: normal_alignment and normal_index must be supplied together')

    config['samples'] = list(samples.values())
    config['batch_import'] = dict(schema_version=1,
        samples_tsv=dict(path=str(samples_path), sha256=digest(samples_path)),
        shared_config=dict(path=str(config_path), sha256=digest(config_path)),
        sample_order=list(samples), sample_count=len(samples),
        result_row_count=sum(map(len, results.values())))
    # The normal manifest validator is authoritative. All paths are already
    # absolute, so its temporary location cannot change either input base.
    with tempfile.TemporaryDirectory(prefix='asconcordance-import-') as temporary:
        candidate = Path(temporary) / 'manifest.json'
        candidate.write_text(json.dumps(config, ensure_ascii=False, allow_nan=False), encoding='utf-8')
        validated, _, _ = load_manifest(candidate)
    payload = json.dumps(validated, indent=2, ensure_ascii=False, allow_nan=False) + '\n'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation also protects against a concurrent import using the
    # same destination. Existing manifests and original inputs are read-only.
    with output_path.open('x', encoding='utf-8') as handle:
        handle.write(payload)
    return validated
