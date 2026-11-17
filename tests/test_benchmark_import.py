"""Synthetic artifact imports and fake AC subprocesses only; no patient data."""
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

from ascompare.benchmark_import import import_result, reclassify_result, reclassify_cached, _verify_ac_runtime, _probe_json
from ascompare.common import atomic_json, digest
from ascompare.parsing import LABELS
from ascompare.versions import _digest
from test_parsing import write_fixture_run


class BenchmarkImportTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name).resolve()
        self.source=write_fixture_run(self.root/'source')
        self.spec={'run_dir':str(self.source),'aa_version':'1.6.r2','ac_version':'2.0.1',
                   'reference_build':'hg19','validated_declared_components':True}
        self.ref=self.root/'ref';(self.ref/'hg19').mkdir(parents=True)
        (self.ref/'hg19/file_list.txt').write_text('SYNTHETIC reference\n')
        self.config={'reference':{'build':'hg19','aa_data_repo':str(self.ref)},'ac_min_size':0,'ac_min_flow':0}
        self.script=self.root/'fake_ac.py'
        self.script.write_text('''import pathlib,sys,shutil
a=sys.argv;out=pathlib.Path(a[a.index('-o')+1]);out.parent.mkdir(exist_ok=True)
for p in pathlib.Path('''+repr(str(self.source/'ac'))+''').iterdir():
 if p.is_dir():shutil.copytree(p,out.parent/p.name)
 else:shutil.copyfile(p,out.parent/p.name)
print('SYNTHETIC AC completed')
''')
        # Intentionally no AA interpreter, AA script or BAM in this runtime.
        self.runtime={'id':'TARGET','ac_python':sys.executable,'ac_script':str(self.script)}
        self.entry={'id':'TARGET','ac':{'version':'1.2.1','entrypoint':'fake_ac.py'},
                    'supported_classes':['ecDNA','BFB','Linear','Complex-non-cyclic']}
        self.verified={'ac':{'version':'1.2.1','help':'--min_size --min_flow'},'supported_labels':self.entry['supported_classes']}

    def tearDown(self):self.temp.cleanup()

    def imported(self,spec=None,support=LABELS,name='import'):
        return import_result(spec or self.spec,'SYNTHETIC','SOURCE',self.root/name,support,{'chr1':10000,'chr2':10000})

    def reclassify(self,source=None):
        source=source or self.imported()
        with patch('ascompare.benchmark_import._verify_ac_runtime',return_value=self.verified):
            return reclassify_result(source,self.runtime,self.entry,self.config,self.root/'reclassified')

    def test_complete_import_copies_only_selected_small_files_and_marks_declarations(self):
        bam=self.source/'aa/patient.bam';bam.write_bytes(b'BAM must never be opened')
        before=digest(bam)
        original_digest=digest
        def guarded_digest(path):
            if Path(path).suffix in ('.bam','.cram'):raise AssertionError('BAM opened')
            return original_digest(path)
        with patch('ascompare.benchmark_import.digest',side_effect=guarded_digest):record=self.imported()
        self.assertEqual(record['status'],'ok');self.assertEqual(record['analysis_kind'],'aa_ac')
        self.assertEqual(record['provenance_kind'],'imported_declared');self.assertTrue(record['comparison_eligible'])
        self.assertFalse(record['original_execution_verified']);self.assertEqual(digest(bam),before)
        self.assertFalse(list((self.root/'import').rglob('*.bam')))
        self.assertEqual(len(record['aa_artifact_sha256']),3)

    def test_missing_profile_is_missing_never_negative(self):
        next((self.source/'ac').glob('*profiles.tsv')).unlink()
        record=self.imported();self.assertEqual(record['status'],'missing');self.assertIsNone(record['detected'])
        self.assertFalse(record['comparison_eligible'])

    def test_missing_graph_and_summary_are_missing(self):
        next((self.source/'aa').glob('*_graph.txt')).unlink()
        self.assertEqual(self.imported()['status'],'missing')
        next((self.source/'aa').glob('*_summary.txt')).unlink()
        self.assertEqual(self.imported(name='import2')['status'],'missing')

    def test_zero_needs_explicit_completion_and_skip_evidence(self):
        zero=write_fixture_run(self.root/'zero',count=0)
        spec={**self.spec,'run_dir':str(zero)}
        self.assertEqual(self.imported(spec)['status'],'not_detected')
        spec.pop('run_dir');spec.update(aa_dir=str(zero/'aa'),ac_dir=str(zero/'ac'))
        record=self.imported(spec,name='no_stages');self.assertEqual(record['status'],'missing')
        self.assertIn('explicit AA completion',record['reason'])

    def test_unknown_component_or_category_support_is_descriptive_only(self):
        spec=dict(self.spec);spec.pop('aa_version')
        record=self.imported(spec);self.assertEqual(record['status'],'ok');self.assertFalse(record['comparison_eligible'])
        record=self.imported(support=None,name='unknown_support')
        self.assertEqual(record['labels'],['ecDNA']);self.assertFalse(record['comparison_eligible'])
        self.assertEqual(record['category_status']['FAN'],'support_unknown')

    def test_failed_component_validation_mapping_is_not_truthy_approval(self):
        spec={**self.spec,'validated_declared_components':{'aa':False,'ac':False}}
        record=self.imported(spec);self.assertFalse(record['comparison_eligible'])
        spec['validated_declared_components']={'aa':'1.6.r2','ac':'2.0.1'}
        self.assertTrue(self.imported(spec,name='validated_mapping')['comparison_eligible'])

    def test_known_failure_stage_is_preserved(self):
        (self.source/'stages.json').write_text(json.dumps({'aa':{'status':'complete'},'ac':{'status':'failed','reason':'SYNTHETIC failure'}}))
        record=self.imported();self.assertEqual(record['status'],'failed');self.assertEqual(record['stage'],'ac')

    def test_malformed_stage_mapping_is_invalid_not_an_unhandled_exception(self):
        (self.source/'stages.json').write_text('[]')
        record=self.imported();self.assertEqual(record['status'],'invalid');self.assertFalse(record['comparison_eligible'])

    def test_artifact_consistency_does_not_invent_execution_exit_status(self):
        spec={**self.spec,'aa_dir':str(self.source/'aa'),'ac_dir':str(self.source/'ac')};spec.pop('run_dir')
        record=self.imported(spec);self.assertEqual(record['status'],'ok')
        self.assertEqual(record['completion_evidence'],'artifact_consistency_only')
        stages=json.loads((self.root/'import/stages.json').read_text())
        self.assertTrue(all('returncode' not in x for x in stages.values()))
        spec['stages_json']=str(self.root/'missing_stages.json')
        self.assertEqual(self.imported(spec,name='missing_stages')['status'],'missing')

    def test_symlink_cannot_smuggle_unselected_data(self):
        target=self.source/'aa/not_an_artifact.bam';target.write_bytes(b'protected')
        (self.source/'aa/bad_amplicon2_graph.txt').symlink_to(target)
        record=self.imported();self.assertEqual(record['status'],'invalid');self.assertIn('Linked',record['reason'])

    def test_size_limit_and_output_overlap_are_refused(self):
        with patch('ascompare.benchmark_import.MAX_FILE_BYTES',1):record=self.imported()
        self.assertEqual(record['status'],'invalid');self.assertIn('limit',record['reason'])
        with self.assertRaises(ValueError):self.imported(name='source/nested_output')

    def test_ac_only_runs_real_synthetic_subprocess_without_aa_environment(self):
        source=self.imported();before=dict(source['aa_artifact_sha256'])
        record=self.reclassify(source)
        self.assertEqual(record['status'],'ok');self.assertEqual(record['analysis_kind'],'ac_only')
        self.assertFalse(record['aa_executed']);self.assertTrue(record['ac_executed']);self.assertTrue(record['ac_runtime_verified'])
        self.assertEqual(record['source_aa_version'],'1.6.r2');self.assertEqual(record['aa_artifact_sha256'],before)
        self.assertEqual(record['source_graph_identity'],source['source_graph_identity'])
        argv=json.loads((self.root/'reclassified/execution.json').read_text())['argv']
        self.assertNotIn('--bam',argv);self.assertNotIn('--run_AA',argv)

    def test_ac_only_source_hash_change_prevents_launch(self):
        source=self.imported();path=next((self.root/'import/aa').glob('*_graph.txt'));path.write_text('changed')
        with patch('ascompare.benchmark_import.execute') as execute:
            record=self.reclassify(source);execute.assert_not_called()
        self.assertEqual(record['status'],'unsupported');self.assertIn('checksum',record['reason'])

    def test_nested_benchmark_ac_parameters_are_used_and_recorded(self):
        self.config['ac_parameters']={'min_size':120000,'min_flow':3.5}
        record=self.reclassify();self.assertEqual(record['status'],'ok')
        execution=json.loads((self.root/'reclassified/execution.json').read_text())
        self.assertEqual(execution['parameters'],self.config['ac_parameters'])
        argv=execution['argv']
        self.assertEqual(argv[argv.index('--min_size')+1],'120000')
        self.assertEqual(argv[argv.index('--min_flow')+1],'3.5')

    def test_zero_graph_is_not_an_ac_negative(self):
        zero=write_fixture_run(self.root/'zero',count=0);source=self.imported({**self.spec,'run_dir':str(zero)})
        record=self.reclassify(source);self.assertEqual(record['status'],'unsupported');self.assertFalse(record['ac_executed'])

    def test_unverified_ac_is_unsupported_and_not_executed(self):
        source=self.imported()
        with patch('ascompare.benchmark_import._verify_ac_runtime',side_effect=ValueError('AC source mismatch')):
            record=reclassify_result(source,self.runtime,self.entry,self.config,self.root/'reclassified')
        self.assertEqual(record['status'],'unsupported');self.assertFalse(record['ac_executed']);self.assertFalse(record['ac_runtime_verified'])

    def test_ac_process_failure_remains_failure(self):
        self.script.write_text('import sys;print("SYNTHETIC crash");sys.exit(3)')
        record=self.reclassify();self.assertEqual(record['status'],'failed');self.assertFalse(record['comparison_eligible'])

    def test_cancellation_records_interruption_and_propagates(self):
        source=self.imported()
        with patch('ascompare.benchmark_import._verify_ac_runtime',return_value=self.verified),patch('ascompare.benchmark_import.execute',side_effect=KeyboardInterrupt):
            with self.assertRaises(KeyboardInterrupt):
                reclassify_result(source,self.runtime,self.entry,self.config,self.root/'reclassified')
        record=json.loads((self.root/'reclassified/result.json').read_text())
        self.assertEqual(record['status'],'failed');self.assertFalse(record['comparison_eligible'])

    def test_exit_zero_with_disabled_required_bfb_is_failure(self):
        self.verified['ac']['version']='2.0.1';self.entry['ac']['version']='2.0.1'
        self.script.write_text(self.script.read_text()+"\nprint('Skipping BFBArchitect integration')\n")
        record=self.reclassify();self.assertEqual(record['status'],'failed');self.assertIn('BFBArchitect',record['reason'])

    def test_ac_verifier_rejects_source_and_version_mismatch_without_aa_keys(self):
        with self.assertRaisesRegex(ValueError,'source manifest'):_verify_ac_runtime(self.runtime,self.entry,os.environ.copy())
        fake={'returncode':0,'stdout':'AC 9.9.9','stderr':''}
        with patch('ascompare.benchmark_import._check_source',return_value={}),patch('ascompare.benchmark_import._run',return_value=fake):
            with self.assertRaisesRegex(ValueError,'version differs'):_verify_ac_runtime(self.runtime,self.entry,os.environ.copy())


class ACOnlyCacheTests(unittest.TestCase):
    """Actual pinned SYNTHETIC source, actual version/package probes and subprocess."""
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name).resolve()
        self.original=write_fixture_run(self.root/'original')
        spec=dict(run_dir=str(self.original),aa_version='SYNTHETIC_AA',ac_version='SYNTHETIC_SOURCE_AC',
                  reference_build='hg19',validated_declared_components=True)
        self.source=import_result(spec,'SYNTHETIC_CACHE','SOURCE',self.root/'import',LABELS,{'chr1':10000,'chr2':10000})
        ac_source=self.root/'pinned_ac_source';ac_source.mkdir()
        self.script=ac_source/'classifier.py';self.counter=self.root/'classification_calls.txt'
        self.script.write_text('''# Explicit SYNTHETIC classifier fixture, not official AC.
import pathlib,sys,shutil
if '--version' in sys.argv:print('AmpliconClassifier 1.2.1');raise SystemExit(0)
if '--help' in sys.argv:print('-i --ref -o --min_size --min_flow');raise SystemExit(0)
with open('''+repr(str(self.counter))+''','a') as handle:handle.write('CLASSIFY\\n')
a=sys.argv;prefix=pathlib.Path(a[a.index('-o')+1])
for p in pathlib.Path('''+repr(str(self.original/'ac'))+''').iterdir():
 if p.is_dir():shutil.copytree(p,prefix.parent/p.name)
 else:shutil.copyfile(p,prefix.parent/p.name)
print('SYNTHETIC classification complete')
''')
        files={'classifier.py':digest(self.script)}
        self.entry={'id':'SYNTHETIC_TARGET','ac':{'entrypoint':'classifier.py','version':'1.2.1',
            'repo':'https://example.invalid/SYNTHETIC_LOCK_ONLY','commit':'SYNTHETIC_PINNED_SOURCE','manifest_sha256':_digest(files)},
            'supported_classes':['ecDNA','BFB','Linear','Complex-non-cyclic']}
        atomic_json(ac_source/'.ascompare-source.json',{'repo':self.entry['ac']['repo'],'commit':self.entry['ac']['commit'],'files':files})
        self.runtime={'id':self.entry['id'],'ac_python':sys.executable,'ac_script':str(self.script)}
        self.reference=self.root/'reference';(self.reference/'hg19').mkdir(parents=True)
        (self.reference/'hg19/file_list.txt').write_text('SYNTHETIC local reference\n')
        self.config={'reference':{'build':'hg19','aa_data_repo':str(self.reference)},'ac_parameters':{'min_size':1,'min_flow':1}}
        self.job=self.root/'job'

    def tearDown(self):self.temp.cleanup()

    def run_cached(self):return reclassify_cached(self.source,self.runtime,self.entry,self.config,self.job)

    def calls(self):return self.counter.read_text().splitlines() if self.counter.exists() else []

    def test_actual_runtime_probe_and_resume_reparse_ignore_normalized_tampering(self):
        first=self.run_cached();self.assertEqual(first['status'],'ok',first.get('reason'));self.assertFalse(first['cache_hit'])
        self.assertEqual(self.calls(),['CLASSIFY'])
        cache=self.job/'cache.json';data=json.loads(cache.read_text())
        data['record'].update(labels=['BFB'],version='FORGED',aa_version='FORGED',source_aa_version='FORGED',
            component_versions={'aa':'FORGED','ac':'FORGED'},source_metadata={'source_aa_version':'FORGED'})
        atomic_json(cache,data)
        result_path=Path(first['run_dir'])/'result.json';raw=json.loads(result_path.read_text());raw['labels']=['FAN'];atomic_json(result_path,raw)
        second=self.run_cached();self.assertTrue(second['cache_hit']);self.assertEqual(self.calls(),['CLASSIFY'])
        self.assertEqual(second['labels'],['ecDNA']);self.assertEqual(second['version'],self.entry['id'])
        self.assertEqual(second['component_versions'],{'aa':'SYNTHETIC_AA','ac':'1.2.1'})
        self.assertEqual(second['source_metadata']['source_aa_version'],'SYNTHETIC_AA')
        self.assertFalse(second['ac_executed_this_invocation']);self.assertTrue(second['ac_runtime_verified'])
        # Source imports get a fresh directory on public workflow resume; identical
        # graph content must still reuse the previous AC execution.
        moved=self.root/'new_import';shutil.copytree(Path(self.source['run_dir']),moved)
        self.source['run_dir']=str(moved)
        third=self.run_cached();self.assertTrue(third['cache_hit']);self.assertEqual(self.calls(),['CLASSIFY'])
        self.assertEqual(third['source_record_dir'],str(moved))

    def test_source_graph_tampering_fails_instead_of_reusing(self):
        self.assertEqual(self.run_cached()['status'],'ok')
        p=next((Path(self.source['run_dir'])/'aa').glob('*_graph.txt'));p.write_text('changed')
        result=self.run_cached();self.assertEqual(result['status'],'unsupported');self.assertFalse(result['cache_hit'])
        self.assertEqual(self.calls(),['CLASSIFY']);self.assertIn('checksum',result['reason'])

    def test_changed_pinned_ac_source_cannot_be_bypassed_by_manifest_verified_flag(self):
        self.assertEqual(self.run_cached()['status'],'ok')
        self.script.write_text(self.script.read_text()+'\n# changed source\n')
        self.runtime['verified']=True;self.config['verified']={'ready':True}
        result=self.run_cached();self.assertEqual(result['status'],'unsupported');self.assertFalse(result['ac_executed'])
        self.assertEqual(self.calls(),['CLASSIFY']);self.assertIn('changed source',result['reason'])

    def test_parameters_and_reference_identity_require_new_classification(self):
        first=self.run_cached();self.assertEqual(first['status'],'ok')
        self.config['ac_parameters']['min_flow']=2
        second=self.run_cached();self.assertEqual(second['status'],'ok');self.assertFalse(second['cache_hit'])
        self.assertNotEqual(first['run_dir'],second['run_dir'])
        (self.reference/'hg19/file_list.txt').write_text('SYNTHETIC reference changed\n')
        third=self.run_cached();self.assertEqual(third['status'],'ok');self.assertFalse(third['cache_hit'])
        self.assertEqual(self.calls(),['CLASSIFY']*3)

    def test_truncated_raw_output_retries_and_retains_prior_attempt(self):
        first=self.run_cached();p=next((Path(first['run_dir'])/'aa').glob('*_graph.txt'));p.write_text('')
        second=self.run_cached();self.assertEqual(second['status'],'ok');self.assertFalse(second['cache_hit'])
        self.assertNotEqual(second['run_dir'],first['run_dir']);self.assertEqual(p.read_text(),'')
        self.assertEqual(self.calls(),['CLASSIFY']*2)

    def test_solver_status_text_around_json_is_accepted_without_guessing(self):
        self.assertEqual(_probe_json('License status message\n{"version":"1.0.1"}\nSolver shutdown message\n'),{'version':'1.0.1'})
        with self.assertRaises(ValueError):_probe_json('No JSON; dependency failed')


if __name__=='__main__':unittest.main()
