"""SYNTHETIC ONLY: emulates output contracts to test orchestration, never science."""
import argparse
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--bed');p.add_argument('--out');p.add_argument('-o');p.add_argument('-i');p.add_argument('--bam');p.add_argument('--ref');p.add_argument('--runmode');p.add_argument('--extendmode');p.add_argument('--downsample');p.add_argument('--min_size');p.add_argument('--min_flow'); a=p.parse_args()
prefix=Path(a.out or a.o); prefix.parent.mkdir(parents=True,exist_ok=True)
old='SYNTHETIC_history' in str(prefix)
if a.bed:
 chrom,start,end,seed=Path(a.bed).read_text().split(); start=int(start);end=int(end)
 if seed=='negative':
  Path(str(prefix)+'_summary.txt').write_text('#Amplicons = 0\n')
 else:
  Path(str(prefix)+'_summary.txt').write_text('#Amplicons = 1\n')
  Path(str(prefix)+'_amplicon1_graph.txt').write_text(f'sequence\t{chrom}:{start}-\t{chrom}:{end-1}+\t8\t20\t{end-start}\t10\n')
  Path(str(prefix)+'_amplicon1_cycles.txt').write_text(f'Segment\t1\t{chrom}\t{start}\t{end-1}\nCycle=1;Copy_count=8;Segments=1+\n')
else:
 fields=Path(a.i).read_text().splitlines()[0].split('\t'); seg=Path(fields[1]).read_text().splitlines()[0].split(); chrom,start,end=seg[2],int(seg[3]),int(seg[4])+1
 labels=['BFB'] if old else ['ecDNA','BFB']
 header='sample_name\tamplicon_number\tamplicon_decomposition_class\tecDNA+\tBFB+\tFAN+\n'
 row='reconstruction\tamplicon1\tCyclic\t'+('None detected' if old else 'Positive')+'\tPositive\tNone detected\n'
 Path(str(prefix)+'_amplicon_classification_profiles.tsv').write_text(header+row)
 beddir=prefix.parent/(prefix.name+'_classification_bed_files');beddir.mkdir(exist_ok=True)
 for label in labels:
  (beddir/f'reconstruction_amplicon1_{label}_1_intervals.bed').write_text(f'{chrom}\t{start}\t{end}\n')
