"""Tiny generated BAM + explicit fixture caller; no actual AA reconstruction."""
import sys
from pathlib import Path
from .common import atomic_json,table
from .inputs import validate
from .runner import run_jobs


def demo(directory):
    import pysam
    from .reporting import write_reports
    out=Path(directory).resolve();out.mkdir(parents=True,exist_ok=True)
    data=out/'synthetic_inputs';data.mkdir(exist_ok=True)
    ref=data/'reference'/'hg19';ref.mkdir(parents=True,exist_ok=True)
    fasta=ref/'toy.fa';fasta.write_text('>chr1\n'+'ACGT'*25000+'\n');pysam.faidx(str(fasta))
    (ref/'file_list.txt').write_text('fa_file toy.fa\nchrLen_file toy.fa.fai\n')
    bed=data/'seeds.bed';bed.write_text('chr1\t100\t1000\tpositive\nchr1\t2000\t3000\tnegative\n')
    rows=[]
    for n in range(1,4):
        bam=data/f'SYNTHETIC_{n}.bam'
        with pysam.AlignmentFile(str(bam),'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':100000}]}) as f:
            read=pysam.AlignedSegment();read.query_name='synthetic_read';read.query_sequence='A'*100;read.flag=0;read.reference_id=0;read.reference_start=150;read.mapping_quality=60;read.cigar=[(0,100)];read.query_qualities=pysam.qualitystring_to_array('I'*100);f.write(read)
        pysam.index(str(bam))
        rows.append(dict(sample_id=f'SYNTHETIC_{n}',alignment=str(bam),index=str(bam)+'.bai',seeds_bed=str(bed),reference_build='hg19'))
    fixture=str(Path(__file__).with_name('fixture_caller.py'))
    versions=[dict(id=name,aa_python=sys.executable,ac_python=sys.executable,aa_script=fixture,ac_script=fixture) for name in ('SYNTHETIC_latest','SYNTHETIC_history')]
    config={'synthetic':True,'reference_version':'SYNTHETIC_latest','versions':versions,'reference':{'build':'hg19','fasta':str(fasta),'aa_data_repo':str(ref.parent)},'statistics':{'min_samples':20,'bootstrap':500,'permutations':1000}}
    checked,sources=validate(config,rows);atomic_json(out/'input_check.json',checked);table(data/'samples.tsv',rows)
    if any(x['status']!='ok' for x in checked):
        raise RuntimeError('Synthetic preflight failed')
    atomic_json(out/'config.resolved.json',config);atomic_json(out/'provenance.json',{'synthetic':True,'scope':'fixture only; no biological validation','inputs':sources})
    runtimes={x['id']:{'ready':True,'aa':{'help':''},'ac':{'version':'fixture'},'supported_labels':['ecDNA','BFB','Complex-non-cyclic','Linear']} for x in versions}
    records=run_jobs(config,checked,runtimes,out);atomic_json(out/'records.json',records)
    if any(x['status'] not in ('ok','not_detected') for x in records):
        raise RuntimeError('Synthetic runner/parser did not pass; inspect records.json')
    write_reports(out,records,config['reference_version'],options=config['statistics'],synthetic=True)
    print('SYNTHETIC ONLY demo complete: '+str(out))
