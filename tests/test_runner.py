"""Actual fixture subprocesses test orchestration, never biological accuracy."""
import json
import os
import sys
import tempfile
import unittest
import venv
from pathlib import Path
from unittest.mock import patch

from ascompare.runner import execute, run_jobs, check_ac_runtime


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.ref = self.root / 'reference'
        self.ref.mkdir()
        self.out = self.root / 'output'
        fixture = Path(__file__).resolve().parents[1] / 'ascompare/fixture_caller.py'
        self.runtime = dict(id='SYNTHETIC_test', aa_python=sys.executable,
                            ac_python=sys.executable, aa_script=str(fixture), ac_script=str(fixture))
        self.config = dict(synthetic=True, versions=[self.runtime], reference_version='SYNTHETIC_test',
                           reference=dict(aa_data_repo=str(self.ref), build='hg19', fasta='unused'))
        self.samples = [dict(sample_id='TEST_SAMPLE', status='ok', alignment='unused.bam', seeds_bed='unused.bed',
                             reference_lengths={'chr1':10000}, seeds=[dict(seed_id='positive',chrom='chr1',start=100,end=1000)])]
        self.verified = {'SYNTHETIC_test':dict(ready=True, aa={'help':''}, ac={'version':'fixture'},
                                             supported_labels=['ecDNA','BFB','Linear','Complex-non-cyclic'])}
        # Fixture processes do not read the alignment. Production preflight and BAM/CRAM
        # staging are covered separately; this patch does not replace process execution.
        self.stage = patch('ascompare.runner.prepare_alignment', return_value='unused.bam')
        self.stage.start()

    def tearDown(self):
        self.stage.stop()
        self.temp.cleanup()

    def run_jobs(self):
        return run_jobs(self.config,self.samples,self.verified,self.out)

    def test_subprocess_success_and_resume(self):
        first = self.run_jobs()[0]
        self.assertEqual(first['status'],'ok')
        self.assertEqual(set(first['labels']),{'ecDNA','BFB'})
        second = self.run_jobs()[0]
        self.assertEqual(first['run_dir'],second['run_dir'])
        self.assertEqual(len(list(Path(first['run_dir']).parent.glob('attempt_*'))),1)

    def test_edited_cached_classification_is_reparsed(self):
        first = self.run_jobs()[0]
        saved = Path(first['run_dir']).parent / 'result.json'
        first['labels'] = ['FAN']
        saved.write_text(json.dumps(first))
        second = self.run_jobs()[0]
        self.assertEqual(set(second['labels']),{'ecDNA','BFB'})

    def test_truncated_output_retries(self):
        first = self.run_jobs()[0]
        graph = next((Path(first['run_dir'])/'aa').glob('*_graph.txt'))
        graph.write_text('')
        second = self.run_jobs()[0]
        self.assertEqual(second['status'],'ok')
        self.assertNotEqual(first['run_dir'],second['run_dir'])

    def test_corrupt_result_json_retries(self):
        first = self.run_jobs()[0]
        (Path(first['run_dir']).parent/'result.json').write_text('{')
        second = self.run_jobs()[0]
        self.assertNotEqual(first['run_dir'],second['run_dir'])

    def test_failed_job_retries_with_preserved_attempt(self):
        flag = self.root/'fail_flag';flag.touch()
        wrapper = self.root/'aa_wrapper.py'
        wrapper.write_text('import pathlib,runpy,sys\n'
                           f'if pathlib.Path({str(flag)!r}).exists(): sys.exit(7)\n'
                           f'runpy.run_path({self.runtime["aa_script"]!r},run_name="__main__")\n')
        self.runtime['aa_script'] = str(wrapper)
        first = self.run_jobs()[0]
        self.assertEqual((first['status'],first['stage']),('failed','aa'))
        flag.unlink()
        second = self.run_jobs()[0]
        self.assertEqual(second['status'],'ok')
        self.assertNotEqual(first['run_dir'],second['run_dir'])
        self.assertTrue(Path(first['run_dir']).is_dir())

    def test_zero_amplicons_skips_ac(self):
        self.samples[0]['seeds'][0]['seed_id']='negative'
        record = self.run_jobs()[0]
        self.assertEqual(record['status'],'not_detected')
        self.assertFalse((Path(record['run_dir'])/'ac.log').exists())

    def test_empty_successful_aa_not_negative(self):
        script = self.root/'empty.py';script.write_text('pass\n')
        self.runtime['aa_script']=str(script)
        record = self.run_jobs()[0]
        self.assertEqual(record['status'],'missing')
        self.assertFalse((Path(record['run_dir'])/'ac.log').exists())

    def test_missing_reference_extension_not_success(self):
        self.samples[0]['reference_lengths']={'chr1':500}
        record = self.run_jobs()[0]
        self.assertEqual((record['status'],record['stage']),('invalid','aa'))

    def test_preflight_failure_preserved(self):
        self.samples[0].update(status='failed',reason='Reference mismatch')
        record = self.run_jobs()[0]
        self.assertEqual((record['status'],record['stage'],record['seed_id']),('failed','input_check','positive'))
        self.assertTrue((self.out/'anomalies.tsv').is_file())

    def test_timeout_terminates_real_subprocess(self):
        pidfile=self.root/'pid.txt'
        script=f'import os,time,pathlib; pathlib.Path({str(pidfile)!r}).write_text(str(os.getpid())); time.sleep(10)'
        with self.assertRaises(TimeoutError):
            execute([sys.executable,'-c',script],self.root,os.environ.copy(),self.root/'timeout.log',timeout=1.0)
        self.assertTrue(pidfile.is_file())
        with self.assertRaises(ProcessLookupError):
            os.kill(int(pidfile.read_text()),0)

    def test_venv_python_and_child_python_preserved_for_each_component(self):
        # Real isolated venvs with symlinked interpreters catch accidental resolve().
        fixture = self.runtime['aa_script']
        for component in ('aa','ac'):
            env_path = self.root / (component + '_venv')
            venv.EnvBuilder(with_pip=False,symlinks=True).create(env_path)
            wrapper = self.root / (component + '_wrapper.py')
            wrapper.write_text('import sys,subprocess,os,runpy\n'
                               f'assert os.path.realpath(sys.prefix) == {str(env_path.resolve())!r}, (sys.prefix, sys.executable)\n'
                               f'assert os.path.realpath(subprocess.check_output(["python","-c","import sys;print(sys.prefix)"],text=True).strip()) == {str(env_path.resolve())!r}\n'
                               'assert os.environ["ASCOMPARE_TEST_ENV"] == "passed"\n'
                               f'runpy.run_path({fixture!r},run_name="__main__")\n')
            self.runtime[component+'_python'] = str(env_path/'bin/python')
            self.runtime[component+'_script'] = str(wrapper)
        self.runtime['env']={'ASCOMPARE_TEST_ENV':'passed'}
        record = self.run_jobs()[0]
        self.assertEqual(record['status'],'ok',record.get('reason'))

    def test_bfb_normal_candidate_skip_not_integration_failure(self):
        log='BFBArchitect 1.0.1 integration enabled.\nSkipping BFBArchitect whole-graph retry because graph preflight failed\n'
        check_ac_runtime(log,'2.0.0')

    def test_bfb_integration_must_be_confirmed(self):
        for log in ('', 'BFBArchitect is installed but centromeres could not be loaded. Skipping BFBArchitect.',
                    'BFBArchitect 1.0.1 integration enabled.\nBFBArchitect failed for graph: error'):
            with self.subTest(log=log),self.assertRaises(RuntimeError):
                check_ac_runtime(log,'2.0.0')


if __name__=='__main__':
    unittest.main()
