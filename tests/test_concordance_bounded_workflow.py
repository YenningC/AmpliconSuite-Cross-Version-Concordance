"""Actual guarded synthetic AA/AC subprocesses on newly-created toy BAMs.

Component verification and optional local CNV are explicit SYNTHETIC stubs.
No official component execution or patient input is used in these tests.
"""
import copy
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam

from ascompare.concordance import (_bridge,fast_fingerprint,load_manifest,project_record,
                                 run_concordance,run_bounded_downstream)
from ascompare.common import atomic_json,digest
from ascompare.full_suite import run_full_suite
from ascompare.parsing import parse_run
from ascompare.statistics import AGREEMENT_METRICS
from ascompare.versions import load_lock
from test_parsing import write_fixture_run


class BoundedConcordanceWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name).resolve()
        self.bam=self.root/'SYNTHETIC_tumor.bam'
        header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':10000},{'SN':'chr2','LN':20000}]}
        with pysam.AlignmentFile(self.bam,'wb',header=header) as handle:
            for number,start in enumerate((105,125,155,800)):
                read=pysam.AlignedSegment();read.query_name=f'SYNTHETIC_{number}';read.query_sequence='A'*30
                read.flag=0;read.reference_id=0;read.reference_start=start;read.mapping_quality=60
                read.cigar=[(0,30)];read.query_qualities=pysam.qualitystring_to_array('I'*30);handle.write(read)
        pysam.index(str(self.bam))  # Create the NEW toy index, never a user input index.
        self.normal=self.root/'SYNTHETIC_normal.bam';shutil.copyfile(self.bam,self.normal)
        shutil.copyfile(str(self.bam)+'.bai',str(self.normal)+'.bai')
        self.fasta=self.root/'SYNTHETIC.fa'
        self.fasta.write_text('>chr1\n'+'A'*10000+'\n>chr2\n'+'A'*20000+'\n');pysam.faidx(str(self.fasta))
        self.repo=self.root/'reference';(self.repo/'hg19').mkdir(parents=True)
        (self.repo/'hg19/file_list.txt').write_text('SYNTHETIC reference only\n')
        self.template=write_fixture_run(self.root/'template')
        self.source=self.root/'source';shutil.copytree(self.template,self.source)
        self.seeds=[dict(chrom='chr1',start=110,end=130,seed_id='shared_1'),
                    dict(chrom='chr1',start=150,end=180,seed_id='shared_2')]
        self.seedfile=self.root/'SYNTHETIC_seeds.bed';self.seedfile.write_text('chr1\t110\t130\nchr1\t150\t180\n')
        self.counter=self.root/'synthetic_execution_trace.txt'
        self.aa=self.root/'synthetic_aa.py';self.ac=self.root/'synthetic_ac.py'
        self.aa.write_text('''import os,pathlib,shutil,sys,pysam
a=sys.argv
with open(os.environ['SYNTHETIC_TRACE'],'a') as f:f.write('AA\\n')
bam=pysam.AlignmentFile(a[a.index('--bam')+1],'rb')
mode=os.environ.get('SYNTHETIC_FETCH_MODE','finite')
prefix=pathlib.Path(a[a.index('--out')+1])
if mode!='finite':
 try:
  if mode=='unbounded':list(bam.fetch())
  else:
   for _ in range(50):list(bam.fetch('chr1',90,250))
 except RuntimeError:pass
 # Deliberately attempt to disguise a caught budget failure as zero amplicons.
 pathlib.Path(str(prefix)+'_summary.txt').write_text('#Amplicons = 0\\n')
else:
 assert sum(1 for _ in bam.fetch('chr1',90,250))==3
 assert list(bam.fetch('chr2',990,1110))==[]  # Finite connected-region query is permitted.
 for p in pathlib.Path('''+repr(str(self.template/'aa'))+''').glob('*.txt'):
  shutil.copyfile(p,prefix.parent/p.name)
''')
        self.ac.write_text('''import os,pathlib,shutil,sys
a=sys.argv
with open(os.environ['SYNTHETIC_TRACE'],'a') as f:f.write('AC\\n')
rows=pathlib.Path(a[a.index('-i')+1]).read_text().splitlines()
assert len(rows)==1 and all(pathlib.Path(p).is_file() for p in rows[0].split()[1:])
prefix=pathlib.Path(a[a.index('-o')+1])
for p in pathlib.Path('''+repr(str(self.template/'ac'))+''').iterdir():
 if p.is_dir():shutil.copytree(p,prefix.parent/p.name)
 else:shutil.copyfile(p,prefix.parent/p.name)
pathlib.Path(str(prefix)+'.log').write_text('SYNTHETIC - INFO - Complete\\n')
print('BFBArchitect SYNTHETIC integration enabled')
''')
        self.entries={e['id']:e for e in load_lock()['versions']}
        self.entry=self.entries['suite-1.6.0']
        self.runtime=self.make_runtime('suite-1.6.0')
        self.verified=self.make_verified(self.runtime,self.entry)
        self.policy=self.root/'policy.json'
        self.budgets={'max_fetch_bp':1000,'max_total_fetch_bp':3000,'max_total_reads':100}
        atomic_json(self.policy,dict(self.budgets,alignments={str(self.bam):str(self.bam)+'.bai'}))
        self.sample=dict(sample_id='SYNTHETIC_BOUNDED',source_version='suite-1.6.0',reference_build='hg19',
            alignment=str(self.bam),index=str(self.bam)+'.bai',normal_alignment=str(self.normal),normal_index=str(self.normal)+'.bai',
            seeds=self.seeds,_input_identities=[fast_fingerprint(self.bam),{'path':str(self.bam)+'.bai','sha256':digest(str(self.bam)+'.bai')}])
        self.config=dict(reference={'build':'hg19','fasta':str(self.fasta),'aa_data_repo':str(self.repo)},
            _reference_lengths={'chr1':10000,'chr2':20000},_io_policy_path=str(self.policy),
            io_policy=self.budgets,aa_parameters={},ac_parameters={'min_size':1,'min_flow':1},job_timeout_seconds=30)

    def tearDown(self):self.temp.cleanup()

    def make_runtime(self,version):
        return dict(id=version,aa_python=sys.executable,aa_script=str(self.aa),ac_python=sys.executable,ac_script=str(self.ac),
            solver='clarabel',env={'SYNTHETIC_TRACE':str(self.counter),'SYNTHETIC_FETCH_MODE':'finite'})

    @staticmethod
    def make_verified(runtime,entry,env=None):
        return dict(ready=True,diagnostics=[],validation_level='SYNTHETIC_VERIFICATION_STUB_NOT_OFFICIAL_COMPONENTS',
            aa={'version':entry['aa']['version'],'help':'--solver --no_cstats','manifest_sha256':'SYNTHETIC_AA'},
            ac={'version':entry['ac']['version'],'manifest_sha256':'SYNTHETIC_AC'})

    def run_downstream(self,name='bounded_job'):
        return run_bounded_downstream(self.sample,self.seeds,self.runtime,self.entry,self.verified,self.config,self.root/name)

    def trace(self):return self.counter.read_text().splitlines() if self.counter.exists() else []

    def test_actual_guarded_call_resume_and_cache_labels_are_reparsed_from_raw(self):
        record=self.run_downstream();self.assertEqual(record['status'],'ok',record.get('reason'))
        self.assertEqual(self.trace(),['AA','AC'])
        self.assertIn(['chr2',1000,1100],record['reconstruction_intervals'])
        audits=list(Path(record['run_dir']).glob('bam_access.*.json'));self.assertEqual(len(audits),1)
        audit=json.loads(audits[0].read_text());self.assertEqual(audit['queries'],2);self.assertEqual(audit['yielded_reads'],3)
        self.assertEqual(audit['blocked'],[])
        cached=self.run_downstream();self.assertEqual(cached['run_dir'],record['run_dir'])
        cache=self.root/'bounded_job/cache.json';data=json.loads(cache.read_text())
        data['record']['labels']=['BFB'];data['record']['intervals']=[];atomic_json(cache,data)
        corrected=self.run_downstream();self.assertEqual(corrected['labels'],['ecDNA'])
        self.assertEqual(corrected['intervals'],record['intervals']);self.assertEqual(self.trace(),['AA','AC'])

    def test_truncated_graph_creates_new_attempt_and_does_not_reuse_cache(self):
        first=self.run_downstream();graph=next((Path(first['run_dir'])/'aa').glob('*_graph.txt'));graph.write_text('')
        second=self.run_downstream();self.assertEqual(second['status'],'ok')
        self.assertTrue(second['run_dir'].endswith('attempt_0002'));self.assertEqual(self.trace(),['AA','AC','AA','AC'])
        self.assertEqual(graph.read_text(),'')  # Failed prior evidence is retained.

    def test_caught_unbounded_and_budget_failures_cannot_become_negative(self):
        for mode in ('unbounded','budget'):
            with self.subTest(mode=mode):
                self.runtime['env']['SYNTHETIC_FETCH_MODE']=mode
                record=self.run_downstream(mode)
                self.assertEqual(record['status'],'failed');self.assertEqual(record['stage'],'aa')
                self.assertFalse(record['comparison_eligible'])
                audit=json.loads(next(Path(record['run_dir']).glob('bam_access.*.json')).read_text())
                self.assertTrue(audit['blocked'])
                self.assertEqual((Path(record['run_dir'])/'aa/reconstruction_summary.txt').read_text(),'#Amplicons = 0\n')
        self.assertEqual(self.trace(),['AA','AA'])

    def test_ac_exit_zero_with_traceback_or_missing_final_marker_is_not_success(self):
        original=self.ac.read_text()
        cases={'traceback':original+"\nprint('Traceback (most recent call last): SYNTHETIC hidden failure')\n",
               'no_completion':original.replace(' - INFO - Complete',' - INFO - Still running')}
        for name,body in cases.items():
            with self.subTest(name=name):
                self.ac.write_text(body)
                record=self.run_downstream(name)
                self.assertEqual(record['status'],'failed')
                self.assertEqual(record['stage'],'ac')
                self.assertFalse(record['comparison_eligible'])

    def test_source_bridge_uses_four_metrics_and_keeps_shared_remote_structure(self):
        replay=self.run_downstream()
        original=parse_run(self.source,self.sample['sample_id'],'__source__',self.runtime['id'],self.entry['supported_classes'])
        original['comparison_eligible']=True
        rows=_bridge(original,replay,self.seeds);self.assertEqual(len(rows),2)
        self.assertTrue(all(r['passed'] and all(r[m]==0 for m in AGREEMENT_METRICS) for r in rows))
        for seed in self.seeds:
            projected=project_record(replay,seed)
            self.assertEqual(projected['reconstruction_intervals'],replay['reconstruction_intervals'])
            self.assertEqual(len(projected['amplicons']),1)
            self.assertIn(['chr2',1000,1100],projected['intervals'])
        altered=copy.deepcopy(original)
        altered['labels']=['BFB']
        for amp in altered['amplicons']:amp['labels']=['BFB']
        for feature in altered['features']:feature['label']='BFB'
        rows=_bridge(altered,replay,self.seeds)
        self.assertTrue(all(not r['passed'] and r['label_discordance']==1 and r['class_footprint_distance']==1 for r in rows))
        self.assertTrue(all(r['detection_discordance']==0 and r['interval_distance']==0 for r in rows))

    def test_legacy_cli_and_full_suite_good_sample_entry_are_disabled(self):
        out=self.root/'legacy_output'
        result=subprocess.run([sys.executable,'-m','ascompare','run','--config','missing.json','--samples','missing.tsv','--out',str(out)],
                              capture_output=True,text=True,cwd=Path(__file__).resolve().parents[1])
        self.assertEqual(result.returncode,2);self.assertIn('Legacy BAM workflows are disabled',result.stderr)
        self.assertFalse(out.exists())
        with self.assertRaisesRegex(ValueError,'Whole-genome Suite execution is disabled'):
            run_full_suite({'versions':[]},[{'sample_id':'SYNTHETIC','status':'ok'}],{},self.root/'full_suite')
        self.assertFalse((self.root/'full_suite').exists())

    def run_public_workflow(self,name='public',alter_source=False):
        if alter_source:
            profile=next((self.source/'ac').glob('*profiles.tsv'))
            profile.write_text(profile.read_text().replace('Cyclic\tPositive\tNone detected','Cyclic\tNone detected\tPositive'))
            for bed in (self.source/'ac/classification_classification_bed_files').glob('*_ecDNA_1_intervals.bed'):
                bed.rename(bed.with_name(bed.name.replace('_ecDNA_','_BFB_')))
        manifest=self.root/(name+'.json')
        sample={k:v for k,v in self.sample.items() if not k.startswith('_') and k!='seeds'}
        sample.update(seeds_bed=str(self.seedfile),source_result={'run_dir':str(self.source),
            'component_versions':{k:self.entry[k]['version'] for k in ('aa','ac')}})
        config=dict(mode='with_bam',synthetic=True,reference_version='suite-1.6.2',target_versions=['suite-1.6.2','suite-1.6.0'],
            samples=[sample],runtimes=[self.make_runtime('suite-1.6.0'),self.make_runtime('suite-1.6.2')],
            reference=self.config['reference'],io_policy=self.budgets,ac_parameters={'min_size':1,'min_flow':1},
            sampling={'max_seeds':2,'random_seed':42},statistics={'min_samples':20},threads=1)
        atomic_json(manifest,config);config,lock,entries=load_manifest(manifest)
        inputs=[self.bam,Path(str(self.bam)+'.bai'),self.normal,Path(str(self.normal)+'.bai'),self.seedfile]
        inputs += [p for p in self.source.rglob('*') if p.is_file()]
        before={str(p):digest(p) for p in inputs}
        def cnv(sample,runtime,entry,config,out):
            return [dict(sample_id=sample['sample_id'],seed_id=s['seed_id'],version=runtime['id'],source_version=sample['source_version'],
                status='retained',retention_fraction=1.,comparison_eligible=True,calibration_status='passed',native_recipe_match=True,
                matching_rule='SYNTHETIC_TEST_exact_interval',cnvkit_version='SYNTHETIC_CNV_STUB',synthetic=True) for s in sample['seeds']]
        with patch('ascompare.concordance.verify_runtime',side_effect=self.make_verified), \
             patch('ascompare.concordance_cnv.run_local_upstream',side_effect=cnv),patch('ascompare.concordance_report._plots'):
            result=run_concordance(config,lock,entries,self.root/name)
        self.assertEqual(before,{str(p):digest(p) for p in inputs})
        return result,self.root/name

    def test_public_with_bam_workflow_keeps_real_bounded_calls_and_shared_seed_pairing(self):
        result,out=self.run_public_workflow()
        self.assertEqual(self.trace(),['AA','AC','AA','AC'])
        cfg=json.loads((out/'config.resolved.json').read_text());self.assertEqual(cfg['samples'][0]['source_bridge_status'],'passed')
        records=json.loads((out/'downstream_records.json').read_text());self.assertEqual(len(records),4)
        self.assertTrue(all(r['status']=='ok' for r in records))
        self.assertEqual(len({r['sample_id'] for r in records}),1)
        beds=list((out/'downstream_jobs').glob('*/*/attempt_*/selected_seeds.bed'))
        self.assertEqual(len(beds),2);self.assertEqual(beds[0].read_text(),beds[1].read_text())
        self.assertTrue(all(['chr2',1000,1100] in r['reconstruction_intervals'] for r in records))
        self.assertTrue(all(r['n_selected_samples']==1 for r in result['compatibility_combined']))
        with (out/'sample_similarity.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        self.assertTrue(all(r['n_units']=='2' for r in rows))
        with next((out/'downstream/aa_ac').glob('*/statistics.tsv')).open() as handle:
            stats=list(csv.DictReader(handle,delimiter='\t'))
        signal=[r for r in stats if r['subset']=='signal_positive' and r['metric'] in AGREEMENT_METRICS]
        self.assertTrue(signal);self.assertTrue(all(r['n_samples']=='1' for r in signal))
        self.assertTrue(all(r['decision']=='insufficient_evidence' for r in result['compatibility_combined']))

    def test_public_source_bridge_failure_is_retained_in_final_report(self):
        result,out=self.run_public_workflow(name='bridge_failure',alter_source=True)
        cfg=json.loads((out/'config.resolved.json').read_text());self.assertEqual(cfg['samples'][0]['source_bridge_status'],'failed')
        self.assertTrue(all('input_source_to_rerun_bridge_not_passed' in r['reasons'] for r in result['compatibility_combined']))
        self.assertTrue(all(r['decision']!='within_margin' for r in result['compatibility_combined']))


if __name__=='__main__':unittest.main()
