"""Exercise the actual locked legacy AA background handoff on synthetic inputs.

Runs the native seed filter plus cold/reuse AA on empty and one-seed BEDs.
Checks native acquisition, graph and cycles equality with no patient inputs.
"""
import argparse
import json
import os
from pathlib import Path
import pysam

from ascompare.common import atomic_json, digest
from ascompare.concordance_aa_stats import stage_measurement, validated_values
from ascompare.concordance_policy import guarded_command
from ascompare.runner import execute, local_reference
from ascompare.versions import _check_source, load_lock


def make_fixture(root):
    repo = root/'synthetic_reference'
    genome = repo/'hg19'  # Native CLI requires this build label; all sequences are artificial.
    genome.mkdir(parents=True)
    length = 15000000
    chroms = ['chr1','chr2']
    fasta = genome/'synthetic.fa'
    fasta.write_text(''.join('>'+chrom+'\n'+'A'*length+'\n' for chrom in chroms))
    pysam.faidx(str(fasta))
    (genome/'empty.bed').touch()
    names = ['duke35_filename', 'mapability_exclude_filename', 'gene_filename', 'exon_file',
             'oncogene_filename', 'centromere_filename', 'conserved_regions_filename', 'segdup_filename']
    (genome/'file_list.txt').write_text('fa_file synthetic.fa\nchrLen_file synthetic.fa.fai\n'+
                                        ''.join(key+' empty.bed\n' for key in names))
    bam = root/'SYNTHETIC tumor.bam'
    with pysam.AlignmentFile(bam,'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':chrom,'LN':length} for chrom in chroms]}) as output:
        for contig in range(len(chroms)):
            for start in range(1000,length-1000,400):
                for offset,flag in [(0,99),(180,147)]:
                    read=pysam.AlignedSegment()
                    read.query_name=f'SYNTHETIC_{contig}_{start}'
                    read.flag=flag
                    read.reference_id=read.next_reference_id=contig
                    read.reference_start=start+offset
                    read.next_reference_start=start+180-offset
                    read.query_sequence='A'*100
                    read.cigarstring='100M'
                    read.mapping_quality=60
                    read.template_length=280 if offset==0 else -280
                    output.write(read)
    pysam.index(str(bam))
    return repo,bam


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtimes',required=True,help='local JSON runtime list or config containing runtimes')
    parser.add_argument('--out',required=True,help='fresh, whitespace-free output directory')
    parser.add_argument('--downsample',type=float,default=10,help='AA coverage target for the synthetic reconstruction')
    args=parser.parse_args()
    root=Path(args.out).resolve()
    root.mkdir(parents=True,exist_ok=False)
    repo,bam=make_fixture(root)
    bed=root/'SYNTHETIC_empty_AA_CNV_SEEDS.bed'
    bed.touch()
    policy=root/'synthetic_policy.json'
    atomic_json(policy,dict(alignments={str(bam):str(bam)+'.bai'}))
    runtimes=json.loads(Path(args.runtimes).read_text())
    if isinstance(runtimes,dict):
        runtimes=runtimes['runtimes']
    indexed={r['id']:r for r in runtimes}
    entries={e['id']:e for e in load_lock()['versions']}
    results=[]
    for version in ('suite-1.2.2','suite-1.3.3','suite-1.3.4','suite-1.3.7'):
        runtime=indexed[version]
        aa=Path(runtime['aa_script']).resolve()
        entry=entries[version]['aa']
        lock=_check_source(runtime['aa_source'],entry)
        if aa!=(Path(lock['source_root'])/entry['entrypoint']).resolve():
            raise ValueError('AA entrypoint differs from checked source')
        work=root/version
        work.mkdir()
        env=dict(os.environ,**runtime.get('env',{}))
        env.update(MPLCONFIGDIR=str(root/'mplcache'),PYTHONDONTWRITEBYTECODE='1',PYTHONHASHSEED='0',
                   OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
        result=dict(version=version,aa_version=entry['version'],aa_commit=lock['commit'],synthetic=True,
                    python_hash_seed='0',downsample=args.downsample,scope='background-only and single fixed seed with native AA FULL/EXPLORE; synthetic, no AC')
        for mode in ('producer','cold','reuse','seeded_cold','seeded_reuse'):
            attempt=work/mode
            attempt.mkdir()
            reference=local_reference(repo,attempt/'reference')
            env['AA_DATA_REPO']=str(reference)
            script=aa.with_name('amplified_intervals.py') if mode=='producer' else aa
            alignment=str(bam)
            if mode.endswith('reuse'):
                stats=work/'producer'/'reference'/'coverage.stats'
                measured=dict(path=str(stats),sha256=digest(stats),fields=validated_values(stats,bam))
                alignment=stage_measurement(measured,attempt,bam)
            selected_bed=bed
            if mode.startswith('seeded'):
                selected_bed=attempt/'SYNTHETIC_AA_CNV_SEEDS.bed'
                selected_bed.write_text('chr1\t5100000\t5210000\n')
            argv=['--bam',alignment,'--bed',str(selected_bed),'--out',str(attempt/'synthetic'),'--ref','hg19']
            if not mode.endswith('reuse'):
                argv+=['--no_cstats']
            if mode!='producer':
                argv+=['--runmode','FULL','--extendmode','EXPLORE','--downsample',str(args.downsample)]
            log=attempt/'component.log'
            rc=execute(guarded_command(runtime['aa_python'],script,argv,policy),attempt,env,log,timeout=240)
            audit=json.loads(next(attempt.glob('bam_access.*.json')).read_text())
            result[mode]=dict(returncode=rc,unique_bp=audit['unique_requested_bp'],
                              fetches=audit['queries'],blocked=audit['blocked'])
            if rc or audit['blocked']:
                raise RuntimeError('Native smoke failed: '+str(log))
            if mode.endswith('reuse'):
                result[mode]['reuse_confirmed']='Reusing cstats from '+str(reference/'coverage.stats') in log.read_text()
                if Path(alignment).is_symlink():
                    Path(alignment).unlink()
        result['native_stats_identical']=((work/'producer/reference/coverage.stats').read_bytes()==
                                          (work/'cold/reference/coverage.stats').read_bytes())
        result['seeded_graph_cycles_identical']=all(
            (work/'seeded_cold'/name).stat().st_size>0 and
            digest(work/'seeded_cold'/name)==digest(work/'seeded_reuse'/name)
            for name in ('synthetic_amplicon1_graph.txt','synthetic_amplicon1_cycles.txt'))
        if not (result['reuse']['reuse_confirmed'] and result['seeded_reuse']['reuse_confirmed']
                and result['seeded_graph_cycles_identical'] and result['native_stats_identical'] and
                result['reuse']['unique_bp']<result['cold']['unique_bp']):
            raise RuntimeError('Native smoke assertions failed for '+version)
        results.append(result)
        atomic_json(root/'results.json',results)
        print(version,result['cold']['unique_bp'],result['reuse']['unique_bp'],flush=True)


if __name__=='__main__':
    main()
