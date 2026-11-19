#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Generate the artificial paired-BAM fixture used for local CNVkit comparisons."""

import argparse
import json
from pathlib import Path
import pysam

CHROMOSOMES = ('chr1', 'chr2')
CONTIG_LENGTH = 15_000_000


def write_reference(genome):
    """Make a tiny AA data repository using the build label accepted by native AA."""
    genome.mkdir(parents=True)
    fasta = genome / 'synthetic.fa'
    fasta.write_text(''.join('>' + chrom + '\n' + 'ACGT' * (CONTIG_LENGTH // 4) + '\n' for chrom in CHROMOSOMES))
    pysam.faidx(str(fasta))
    (genome / 'empty.bed').touch()
    (genome / 'synthetic_mappability.bed').write_text(
        ''.join(f'{chrom}\t0\t{CONTIG_LENGTH}\t1\n' for chrom in CHROMOSOMES))
    keys = ['duke35_filename', 'mapability_exclude_filename', 'gene_filename', 'exon_file',
            'oncogene_filename', 'centromere_filename', 'conserved_regions_filename', 'segdup_filename']
    (genome / 'file_list.txt').write_text('fa_file synthetic.fa\nchrLen_file synthetic.fa.fai\n' +
                                         ''.join(key + (' synthetic_mappability.bed\n' if key == 'duke35_filename'
                                                        else ' empty.bed\n') for key in keys))
    return fasta


def write_bams(inputs):
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': chrom, 'LN': CONTIG_LENGTH} for chrom in CHROMOSOMES]}
    counts = {}
    for sample in ('normal', 'tumor'):
        destination = inputs / ('SYNTHETIC_' + sample + '.bam')
        count = 0
        with pysam.AlignmentFile(destination, 'wb', header=header) as output:
            for contig in range(len(CHROMOSOMES)):
                positions = [(start + offset, start, offset, flag)
                             for start in range(1000, CONTIG_LENGTH - 1000, 400)
                             for offset, flag in [(0, 99), (180, 147)]]
                if sample == 'tumor' and contig == 0:
                    positions.extend((start + offset, start, offset, flag)
                                     for start in range(5_090_000, 5_220_000, 60)
                                     for offset, flag in [(0, 99), (180, 147)])
                for position, start, offset, flag in sorted(positions):
                    read = pysam.AlignedSegment()
                    read.query_name = f'SYNTHETIC_{sample}_{contig}_{start}'
                    read.flag = flag
                    read.reference_id = read.next_reference_id = contig
                    read.reference_start = position
                    read.next_reference_start = start + 180 - offset
                    read.query_sequence = 'A' * 100
                    read.query_qualities = pysam.qualitystring_to_array('I' * 100)
                    read.cigarstring = '100M'
                    read.mapping_quality = 60
                    read.template_length = 280 if offset == 0 else -280
                    output.write(read)
                    count += 1
        pysam.index(str(destination))
        counts[sample] = count
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True, help='New directory for synthetic inputs')
    args = parser.parse_args()
    out = args.out.resolve()
    try:
        out.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        parser.error('Output directory already exists; choose a new --out path')
    inputs = out / 'inputs'
    repo = inputs / 'synthetic_reference'
    fasta = write_reference(repo/'hg19')
    counts = write_bams(inputs)
    background = inputs/'SYNTHETIC_background.bed'
    background.write_text('chr1\t1000000\t1400000\nchr2\t1000000\t1400000\n')
    seeds = inputs/'SYNTHETIC_AA_CNV_SEEDS.bed'
    seeds.write_text('chr1\t5100000\t5210000\n')

    record = dict(label='SYNTHETIC_ONLY', synthetic=True,
                  description='Artificial BAM/reference/candidates for local CNVkit comparison. '
                              'The hg19 directory label is required by AA; all sequences are synthetic.',
                  generated_read_counts=counts, build_label='hg19', coordinate_system='BED0_half_open',
                  aa_data_repo=str(repo.relative_to(out)), fasta=str(fasta.relative_to(out)),
                  seeds_bed=str(seeds.relative_to(out)), background_bed=str(background.relative_to(out)))
    for name in ('tumor', 'normal'):
        bam = (inputs/f'SYNTHETIC_{name}.bam').relative_to(out)
        record[name + '_bam'] = str(bam)
        record[name + '_index'] = str(bam) + '.bai'
    (out/'fixture.json').write_text(json.dumps(record, indent=2) + '\n')
    print(f'SYNTHETIC_ONLY fixture: {out / "fixture.json"}')


if __name__ == '__main__':
    main()
