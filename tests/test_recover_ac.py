"""Synthetic recovery only: no patient input or official runtime is used."""
import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

from ascompare.common import atomic_json, digest
from ascompare.full_suite import build_command
from ascompare.parsing import LABELS
from ascompare.recover_ac import recover_ac, build_recovery_command, command_signature, _verify_recovery_runtime
from test_parsing import write_fixture_run


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name).resolve()
        self.attempt=self.root/'job/attempt_0001';self.attempt.mkdir(parents=True)
        source=write_fixture_run(self.root/'synthetic_native')
        pipeline=self.attempt/'pipeline';pipeline.mkdir()
        shutil.copytree(source/'aa',pipeline/'sample_AA_results')
        (pipeline/'sample_AA_CNV_SEEDS.bed').write_text('chr1\t100\t200\n')
        cnv=pipeline/'sample_cnvkit_output';cnv.mkdir()
        (cnv/'tumor.cns').write_text('chromosome\tstart\tend\tlog2\nchr1\t0\t1000\t2.5\n')
        (self.attempt/'suite.log').write_text('Running CNVKit batch\nRunning amplified_intervals\nRunning AA\nRunning AC\ndeclare: -A: invalid option\n')
        metadata=dict(AmpliconSuite_pipeline_version='unused',AA_version='1.6.r2',AC_version='2.0.1',cnvkit_version='0.9.14',
            cnvkit_cmd='cnvkit.py batch tumor.bam --normal normal.bam',amplified_intervals_cmd='amplified_intervals.py',AA_cmd='AA',AC_cmd='AC')
        metadata['AmpliconSuite-pipeline_version']='1.6.2'
        atomic_json(pipeline/'sample_run_metadata.json',metadata)
        self.script=self.root/'synthetic_wrapper.py'
        self.script.write_text('''import pathlib,sys,shutil
a=sys.argv;out=pathlib.Path(a[a.index('-o')+1]);aa=pathlib.Path(a[a.index('--completed_AA_runs')+1])
out.mkdir();ac=out/'sample_classification'
shutil.copytree('''+repr(str(source/'ac'))+''',ac)
rows=[]
for c in aa.glob('*_cycles.txt'):
 rows.append('reconstruction\\t'+str(c)+'\\t'+str(c.with_name(c.name.replace('_cycles.txt','_graph.txt'))))
(ac/'sample.input').write_text('\\n'.join(rows)+'\\n')
(ac/'sample.log').write_text('2026 - INFO - Complete\\n')
(out/'sample_finish_flag.txt').write_text('All stages completed\\n')
print('BFBArchitect 1.0.0 integration enabled')
''')
        ref=self.root/'reference';(ref/'hg19').mkdir(parents=True)
        f=ref/'hg19/file_list.txt';f.write_text('synthetic reference\n')
        st=f.stat();atomic_json(self.attempt/'reference_audit_before.json',[dict(path=str(f),size=st.st_size,mtime_ns=st.st_mtime_ns,sha256=digest(f))])
        self.runtime=dict(id='suite-1.6.2',suite_python=sys.executable,suite_script=str(self.script),aa_python=sys.executable,
            ac_python=sys.executable,aa_script=str(self.script),ac_script=str(self.script),cnvkit_script=str(self.script))
        self.verified=dict(ready=True,suite=dict(version='1.6.2',help='--completed_AA_runs --completed_run_metadata --run_AC --no_repo_check --AA_solver'),
            aa=dict(version='1.6.r2'),ac=dict(version='2.0.1'),cnvkit=dict(version='0.9.14'),supported_labels=list(LABELS))
        self.config=dict(versions=[self.runtime],reference=dict(build='hg19',aa_data_repo=str(ref)),threads=1)
        self.sample=dict(sample_id='SYNTHETIC_RECOVERY',status='ok',reference_lengths={'chr1':10000,'chr2':10000})
        atomic_json(self.root/'config.resolved.json',self.config)
        atomic_json(self.root/'runtime_check.json',{'suite-1.6.2':self.verified})
        atomic_json(self.root/'input_check.json',[self.sample])
        argv=build_command(self.config,self.runtime,self.verified,self.attempt,str(self.attempt/'alignment/tumor.bam'),str(self.attempt/'alignment/normal.bam'))
        atomic_json(self.attempt/'execution.json',dict(endpoint='full_pipeline',argv=argv,seed_thresholds={'cngain':4.5}))
        self.rehash()

    def tearDown(self):self.temp.cleanup()

    def rehash(self):
        self.hashes={str(p):digest(p) for p in self.attempt.rglob('*') if p.is_file()}
        atomic_json(self.attempt.parent/'result.json',dict(status='failed',stage='ac',endpoint='full_pipeline',
            sample_id=self.sample['sample_id'],version=self.runtime['id'],run_dir=str(self.attempt),output_sha256=self.hashes))

    def recover(self):
        with patch('ascompare.recover_ac._verify_recovery_runtime',return_value=self.verified):
            return recover_ac(self.attempt,self.root/'config.resolved.json',self.root/'runtime_check.json',
                self.root/'input_check.json',self.root/'recovery',extra_path=['/synthetic/bash5/bin'])

    def test_subprocess_recovery_inherits_and_keeps_original_failure(self):
        before=(self.attempt.parent/'result.json').read_bytes()
        record=self.recover();self.assertEqual(record['status'],'ok')
        self.assertEqual(record['recovery_kind'],'recovered_ac');self.assertEqual(record['endpoint'],'full_pipeline')
        self.assertEqual(record['rerun_stages'],['ac']);self.assertEqual(record['labels'],['ecDNA'])
        self.assertEqual((self.attempt.parent/'result.json').read_bytes(),before)
        self.assertTrue(all(digest(p)==sha for p,sha in self.hashes.items()))
        self.assertFalse((self.root/'recovery/aa').is_symlink())
        self.assertFalse((self.root/'recovery/pipeline/sample_run_metadata.json').exists())
        manifest=json.loads((self.root/'recovery/recovery_manifest.json').read_text())
        self.assertNotIn('--bam',manifest['argv']);self.assertNotIn('--run_AA',manifest['argv'])
        self.assertEqual(manifest['inherited_cns_segment_count'],1)

    def test_changed_inherited_graph_refused(self):
        p=next((self.attempt/'pipeline/sample_AA_results').glob('*_graph.txt'));p.write_text(p.read_text()+'changed\n')
        with self.assertRaisesRegex(ValueError,'changed'):self.recover()
        self.assertFalse((self.root/'recovery').exists())

    def test_missing_checksum_refused(self):
        data=json.loads((self.attempt.parent/'result.json').read_text());data['output_sha256']={}
        atomic_json(self.attempt.parent/'result.json',data)
        with self.assertRaisesRegex(ValueError,'checksum'):self.recover()

    def test_malformed_nonempty_cns_refused(self):
        (self.attempt/'pipeline/sample_cnvkit_output/tumor.cns').write_text('chromosome\tstart\tend\tlog2\nchr1\t0\t100\tnan\n')
        self.rehash()
        with self.assertRaisesRegex(ValueError,'malformed'):self.recover()

    def test_zero_ac_input_is_failed_not_negative(self):
        self.script.write_text(self.script.read_text().replace("(ac/'sample.input').write_text('\\n'.join(rows)+'\\n')","(ac/'sample.input').write_text('')"))
        record=self.recover();self.assertEqual(record['status'],'failed');self.assertEqual(record['stage'],'ac')
        self.assertIn('every inherited AA amplicon',record['reason'])

    def test_recovery_uses_only_offline_completed_aa_branch(self):
        cmd=build_recovery_command(self.config,self.runtime,self.verified,self.root/'new')
        self.assertEqual(cmd[cmd.index('--completed_AA_runs')+1],str(self.root/'new/aa'))
        self.assertIn('--run_AC',cmd);self.assertIn('--no_repo_check',cmd)
        with self.assertRaisesRegex(ValueError,'offline'):
            build_recovery_command(self.config,self.runtime,{'suite':{'help':'--completed_AA_runs --run_AC'}},self.root/'new')

    def test_original_flag_order_can_change_without_changing_scope(self):
        execution=self.attempt/'execution.json';data=json.loads(execution.read_text())
        data['argv'].remove('--no_repo_check');data['argv'].insert(data['argv'].index('--run_AC')+1,'--no_repo_check')
        atomic_json(execution,data);self.rehash()
        self.assertEqual(self.recover()['status'],'ok')

    def test_command_normalization_rejects_unknown_duplicate_missing_and_changed_values(self):
        original=['/env/python','/source/Suite.py','--run_AA','--ref','hg19','--downsample','-1']
        self.assertEqual(command_signature(original),command_signature(original[:2]+['--downsample','-1','--ref','hg19','--run_AA']))
        for suffix in (['--run_AA'],['--upload'],['--cngain']):
            with self.subTest(suffix=suffix):
                with self.assertRaises(ValueError):command_signature(original+suffix)
        changed=list(original);changed[changed.index('hg19')]='GRCh38'
        self.assertNotEqual(command_signature(original),command_signature(changed))
        changed=list(original);changed[0]='/other/python'
        self.assertNotEqual(command_signature(original),command_signature(changed))

    def test_recovery_categories_come_from_lock_not_parser_default(self):
        baseline=copy.deepcopy(self.verified)
        baseline['supported_labels']=['ecDNA','BFB','Linear','Complex-non-cyclic','FAN']
        fresh=copy.deepcopy(baseline);fresh.pop('supported_labels')
        entry={'supported_classes':list(baseline['supported_labels'])}
        with patch('ascompare.recover_ac.verify_runtime',return_value=fresh),patch('ascompare.recover_ac.verify_full_runtime',return_value=fresh):
            actual=_verify_recovery_runtime(self.runtime,baseline,entry,self.config)
        self.assertIn('FAN',actual['supported_labels']);self.assertNotIn('Virus',actual['supported_labels'])
        baseline['supported_labels'].remove('FAN')
        with self.assertRaisesRegex(ValueError,'categories differ'):
            _verify_recovery_runtime(self.runtime,baseline,entry,self.config)


if __name__=='__main__':unittest.main()
