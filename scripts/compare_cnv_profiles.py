#!/usr/bin/env python3
"""Offline, descriptive paired comparison of existing CNVkit profiles on BED candidates."""
import argparse
import csv
import gzip
import hashlib
import math
from collections import defaultdict
from pathlib import Path


def open_text(path):
    return gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()


def read_candidates(path):
    seeds = []
    with open_text(path) as f:
        for number, line in enumerate(f, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.split()
            try:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                if not chrom or start < 0 or end <= start:
                    raise ValueError()
            except (ValueError, IndexError) as exc:
                raise ValueError(f'{path}:{number}: invalid BED interval') from exc
            seeds.append((chrom, start, end))
    if not seeds or len(seeds) != len(set(seeds)):
        raise ValueError('Candidate BED must contain one or more distinct coordinate triples')
    return seeds


def read_profile(path, seeds, stage):
    """Stream the existing table, validate per-chromosome order and retain only candidate overlaps."""
    if not path.name.endswith(('.' + stage, '.' + stage + '.gz')):
        raise ValueError(f'{path} must end in .{stage} or .{stage}.gz')
    pieces, previous, by_chrom = defaultdict(list), {}, defaultdict(list)
    for seed in seeds:
        by_chrom[seed[0]].append(seed)
    with open_text(path) as f:
        table = csv.reader(f, delimiter='\t')
        header = next(table, [])
        if len(set(header)) != len(header) or not {'chromosome', 'start', 'end', 'log2'} <= set(header):
            raise ValueError(f'{path}: invalid CNVkit header')
        c, a, b, v = (header.index(k) for k in ('chromosome', 'start', 'end', 'log2'))
        w = header.index('weight') if 'weight' in header else None
        for line, fields in enumerate(table, 2):
            try:
                if len(fields) != len(header):
                    raise ValueError('wrong number of columns')
                chrom, start, end, log2 = fields[c], int(fields[a]), int(fields[b]), float(fields[v])
                if not chrom or start < 0 or end <= start or not math.isfinite(log2):
                    raise ValueError('invalid interval or log2')
                if start < previous.get(chrom, 0):
                    raise ValueError('overlapping or unsorted intervals on chromosome')
                previous[chrom] = end
                if w is not None:
                    weight = float(fields[w])
                    if not math.isfinite(weight) or weight < 0:
                        raise ValueError('invalid weight')
                    if weight == 0:
                        continue
            except ValueError as exc:
                raise ValueError(f'{path}:{line}: {exc}') from exc
            for seed in by_chrom.get(chrom, ()):
                left, right = max(start, seed[1]), min(end, seed[2])
                if left < right:
                    pieces[seed].append((left, right, log2))
    return pieces


def intersections(left, right):
    i = j = 0
    while i < len(left) and j < len(right):
        a, b, x = left[i]
        c, d, y = right[j]
        start, end = max(a, c), min(b, d)
        if start < end:
            yield start, end, x, y
        if b <= d:
            i += 1
        if d <= b:
            j += 1


def measured_cns(segments, bins):
    return [(start, end, log2) for start, end, log2, _ in intersections(segments, bins)]


def summarize(seed, stage, source, local, gain, build):
    bp, source_total, local_total, source_gain, local_gain = 0, 0., 0., 0, 0
    cut = math.log2(gain) - 1
    for start, end, x, y in intersections(source, local):
        width = end - start
        bp += width
        source_total += width*x
        local_total += width*y
        source_gain += width*(x >= cut)
        local_gain += width*(y >= cut)
    chrom, start, end = seed
    seed_id = 'seed_' + hashlib.sha256(f'{chrom}:{start}-{end}'.encode()).hexdigest()[:12]
    return dict(seed_id=seed_id, chromosome=chrom, start=start, end=end, declared_build=build,
                stage=stage, status='observed' if bp else 'no_paired_bp', candidate_bp=end-start,
                source_observed_bp=sum(b-a for a,b,_ in source), local_observed_bp=sum(b-a for a,b,_ in local),
                paired_bp=bp, paired_candidate_fraction=bp/(end-start),
                source_mean_log2=source_total/bp if bp else None, local_mean_log2=local_total/bp if bp else None,
                paired_mean_delta_log2=(local_total-source_total)/bp if bp else None,
                gain_threshold_cn=gain, source_gain_bp_fraction=source_gain/bp if bp else None,
                local_gain_bp_fraction=local_gain/bp if bp else None)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, epilog=(
        'CNR uses measured, positive-weight bins when weights are provided. CNS is masked by the CNR '
        'support of each run. Means and gain fractions weight paired bases by their overlap lengths. '
        'Gain is 2 * 2**log2 >= the given copy-number threshold. Build identifiers are declared by the '
        'caller. This diagnostic does not authenticate the historical run or establish calibration.'))
    parser.add_argument('--seeds-bed', type=Path, required=True, help='Original candidates, BED0 (optional gzip)')
    for label in ('source', 'local'):
        for stage in ('cnr', 'cns'):
            parser.add_argument('--'+label+'-'+stage, type=Path, required=True,
                                help=f'Existing {label} .{stage} table (optional gzip)')
    for label in ('source', 'local', 'seed'):
        parser.add_argument('--'+label+'-build', required=True,
                            choices=['hg19', 'GRCh37', 'hg38', 'GRCh38', 'GRCh38_viral'],
                            help=f'Exact declared assembly of {label} input; all three must match')
    parser.add_argument('--gain-cn', required=True, type=float, help='Positive copy-number cutoff')
    parser.add_argument('--out', type=Path, required=True, help='New tab-separated output file')
    args = parser.parse_args(argv)
    try:
        if len({args.source_build, args.local_build, args.seed_build}) != 1:
            raise ValueError('Source, local and BED reference builds must match exactly')
        if not math.isfinite(args.gain_cn) or args.gain_cn <= 0:
            raise ValueError('gain-cn must be positive and finite')
        seeds = read_candidates(args.seeds_bed)
        bins_src = read_profile(args.source_cnr, seeds, 'cnr')
        bins_loc = read_profile(args.local_cnr, seeds, 'cnr')
        seg_src = read_profile(args.source_cns, seeds, 'cns')
        seg_loc = read_profile(args.local_cns, seeds, 'cns')
        rows = []
        for seed in seeds:
            for stage, source, local in [('cnr', bins_src[seed], bins_loc[seed]),
                    ('cns', measured_cns(seg_src[seed], bins_src[seed]), measured_cns(seg_loc[seed], bins_loc[seed]))]:
                rows.append(summarize(seed, stage, source, local, args.gain_cn, args.source_build))
        with args.out.open('x', newline='') as f:
            writer = csv.DictWriter(f, list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
    except (ValueError, OSError, EOFError) as exc:
        parser.exit(2, f'error: {exc}\n')


if __name__ == '__main__':
    main()
