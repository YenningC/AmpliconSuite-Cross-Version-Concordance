#!/usr/bin/env python3
"""Compare two completed, explicitly synthetic positive FULL / CYCLES control runs."""

import argparse
import csv
import hashlib
import json
from pathlib import Path


SCIENCE_FIELDS = (
    'status', 'detected', 'labels', 'intervals', 'reconstruction_intervals', 'amplicons',
    'supported_labels', 'unsupported_labels', 'category_status', 'aa_amplicon_count',
    'evaluable_amplicon_count', 'unevaluable_amplicon_count', 'comparison_eligible',
)
NATIVE_FILES = (
    'aa/*_graph.txt', 'aa/*_cycles.txt', 'aa/*_summary.txt',
    'ac/*amplicon_classification_profiles.tsv', 'ac/**/*_intervals.bed',
)
REPORTS = (
    'result.tsv', 'downstream_similarity_samples.tsv', 'compatibility_combined.tsv',
    'upstream_statistics.tsv', 'local_cnv_comparator_seeds.tsv', 'cohort_compatibility_metrics.tsv',
)


def read_json(path):
    return json.loads(path.read_text())


def hashes(directory, pattern):
    return {p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.glob(pattern)) if p.is_file()}


def records(root):
    if read_json(root / 'config.resolved.json').get('synthetic') is not True:
        raise ValueError('This control requires a synthetic run')
    return {(r['sample_id'], r['seed_id'], r['version']): r
            for r in read_json(root / 'downstream_records.json')}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def comparable_table(root, name):
    # Two provenance columns record real paths in different output directories.
    # Match relative paths and preserve every scientific column and row.
    with (root / name).open() as handle:
        table = csv.DictReader(handle, delimiter='\t')
        rows = []
        for row in table:
            for column in ('target_provenance', 'baseline_provenance'):
                if row.get(column):
                    row[column] = Path(row[column]).relative_to(root).as_posix()
            rows.append(row)
        return table.fieldnames, rows


def compare(full, cycles):
    before, after = records(full), records(cycles)
    require(before and before.keys() == after.keys(), 'Different or empty control panels')
    require(len({(sample, seed) for sample, seed, version in before}) == 1,
            'The positive control must contain one sample with one selected seed')
    suites = {}
    for key in sorted(before):
        left, right = before[key], after[key]
        require(left['status'] == right['status'] == 'ok', f'Incomplete positive control {key}')
        for field in SCIENCE_FIELDS:
            require(left.get(field) == right.get(field), f'{key} changed {field}')
        a, b = Path(left['run_dir']), Path(right['run_dir'])
        for pattern in NATIVE_FILES:
            files = hashes(a, pattern)
            require(files and files == hashes(b, pattern), f'{key} changed {pattern}')
        require(list((a / 'aa').glob('*.png')) and not list((b / 'aa').glob('*.png')),
                f'{key} unexpected native figure inventory')
        sa, sb = (read_json(d / 'stages.json') for d in (a, b))
        for stage, mode in ((sa, 'FULL'), (sb, 'CYCLES')):
            argv = stage['aa']['argv']
            require(argv[argv.index('--runmode') + 1] == mode, f'{key} unexpected run mode')
            require(argv[argv.index('--extendmode') + 1] == 'EXPLORE', f'{key} extension changed')
        suites[key[2]] = dict(aa=sb['aa']['execution'], ac=sb['ac']['execution'],
                               classification=right['labels'])
    table_checks = {}
    for name in REPORTS:
        require(comparable_table(full, name) == comparable_table(cycles, name), f'{name} changed')
        table_checks[name] = ('byte_identical' if (full / name).read_bytes() == (cycles / name).read_bytes()
                              else 'identical_after_normalizing_recorded_output_paths')
    figures = {}
    for pattern in ('*.png', 'sample_ladders/*.png'):
        original = hashes(full, pattern)
        require(original and original == hashes(cycles, pattern), f'Concordance figures changed: {pattern}')
        figures.update(original)
    aa_calls = sum(r['aa'] == 'subprocess' for r in suites.values())
    ac_calls = sum(r['ac'] == 'subprocess' for r in suites.values())
    return dict(synthetic=True, optimized_aa_calls=aa_calls, optimized_ac_calls=ac_calls, scope='Positive fixed-seed control on the supplied runs',
                reference=read_json(cycles / 'config.resolved.json')['reference_version'],
                paired_units=len(before), suite_outputs=suites, native_files_compared=NATIVE_FILES,
                existing_report_tables=table_checks, identical_report_png_sha256=figures)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('full_run', type=Path)
    parser.add_argument('cycles_run', type=Path)
    args = parser.parse_args()
    print(json.dumps(compare(args.full_run.resolve(), args.cycles_run.resolve()), indent=2))
