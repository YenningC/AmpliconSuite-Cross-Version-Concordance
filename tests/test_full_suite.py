import json
from pathlib import Path
import tempfile
import unittest
from ascompare.full_suite import build_command, parse_full_result, read_seed_bed, private_reference, audit_reference, infer_stage, partial_seed_evidence, run_full_suite, offline_policy, require_ac_shell, PipelineStageError


class FullSuiteTests(unittest.TestCase):
    def fixture(self, p):
        pipeline=p/'pipeline';pipeline.mkdir()
        (pipeline/'sample_finish_flag.txt').write_text('All stages completed\n')
        (pipeline/'sample_AA_CNV_SEEDS.bed').write_text('')
        c=pipeline/'sample_cnvkit_output';c.mkdir();(c/'tumor.cns').write_text('chromosome\tstart\tend\tlog2\nchr1\t0\t1000\t0\n')
        aa=pipeline/'sample_AA_results';aa.mkdir();(aa/'sample_summary.txt').write_text('#Amplicons = 0\n')
        (p/'suite.log').write_text('Running CNVKit batch\nRunning amplified_intervals\nRunning AA\n')
        metadata={'AmpliconSuite-pipeline_version':'1.6.2','AA_version':'1.6.r2','AC_version':'2.0.1','cnvkit_version':'0.9.14','cnvkit_cmd':'python cnvkit.py batch tumor.bam --normal normal.bam','amplified_intervals_cmd':'amplified_intervals.py','AA_cmd':'AmpliconArchitect.py','AC_cmd':'NA'}
        (pipeline/'sample_run_metadata.json').write_text(json.dumps(metadata))
        verified={'suite':{'version':'1.6.2'},'aa':{'version':'1.6.r2'},'ac':{'version':'2.0.1'},'cnvkit':{'version':'0.9.14'},'supported_labels':['ecDNA','BFB','Linear','Complex-non-cyclic','FAN']}
        return pipeline,verified

    def test_completed_no_seed_is_only_full_endpoint_negative(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);pipeline,v=self.fixture(p)
            r=parse_full_result(p,{'sample_id':'TEST'},{'id':'suite-1.6.2'},v)
            self.assertEqual(r['status'],'not_detected')
            self.assertEqual(r['endpoint'],'full_pipeline')
            self.assertEqual(r['seed_id'],'__whole_sample__')
            self.assertEqual(r['detection_stage'],'seed_processing')
            self.assertTrue(r['aa_executed'])
            self.assertFalse(r['aa_reconstruction_performed'])
            self.assertFalse(r['ac_executed'])

    def test_missing_cnv_cannot_be_negative(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);pipeline,v=self.fixture(p)
            (pipeline/'sample_cnvkit_output'/'tumor.cns').unlink()
            with self.assertRaises(FileNotFoundError):parse_full_result(p,{'sample_id':'TEST'},{'id':'suite-1.6.2'},v)

    def test_malformed_cnv_cannot_be_negative_or_seed_evaluable(self):
        for segment in ('chr1\t100\t0\t0', 'chr1\t0\t100\tnan', 'chr1\t0\t100\tbroken'):
            with self.subTest(segment=segment),tempfile.TemporaryDirectory() as d:
                p=Path(d);pipeline,v=self.fixture(p)
                (pipeline/'sample_cnvkit_output'/'tumor.cns').write_text('chromosome\tstart\tend\tlog2\n'+segment+'\n')
                with self.assertRaises(ValueError):parse_full_result(p,{'sample_id':'TEST'},{'id':'suite-1.6.2'},v)
                self.assertFalse(partial_seed_evidence(p)['seed_stage_evaluable'])

    def test_exit_zero_without_finish_flag_is_failure(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);pipeline,v=self.fixture(p)
            (pipeline/'sample_finish_flag.txt').write_text('UNSUCCESSFUL\n')
            with self.assertRaises(RuntimeError):parse_full_result(p,{'sample_id':'TEST'},{'id':'suite-1.6.2'},v)

    def test_normal_not_silently_ignored(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);pipeline,v=self.fixture(p)
            f=pipeline/'sample_run_metadata.json';m=json.loads(f.read_text());m['cnvkit_cmd']='python cnvkit.py batch tumor.bam';f.write_text(json.dumps(m))
            with self.assertRaises(ValueError):parse_full_result(p,{'sample_id':'TEST'},{'id':'suite-1.6.2'},v)

    def test_command_never_reuses_external_seeds_or_uploads(self):
        runtime={'suite_python':'/env/python','suite_script':'/src/Suite.py','aa_python':'/env/python','ac_python':'/env/python','aa_script':'/aa/AmpliconArchitect.py','cnvkit_script':'/env/cnvkit.py'}
        command=build_command({'reference':{'build':'hg19'},'solver':'clarabel'},runtime,{'suite':{'help':'--no_repo_check --no_cstats --AA_solver'}},'/safe/run','/safe/tumor.bam','/safe/normal.bam')
        for flag in ('--cnv_bed','--no_filter','--completed_AA_runs','--upload','--fetch_cram_refs'):
            self.assertNotIn(flag,command)
        self.assertIn('--normal_bam',command)
        self.assertIn('--no_repo_check',command)
        with self.assertRaises(ValueError):build_command({'reference':{'build':'hg19'}},runtime,{'suite':{'help':'--no_repo_check --AA_solver'}},'/unsafe path','/tumor','/normal')

    def test_reference_control_files_are_private(self):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d);source=base/'source';(source/'hg19').mkdir(parents=True)
            (source/'hg19'/'file_list.txt').write_text('fa_file x.fa\n');(source/'hg19'/'x.fa').write_text('>chr1\nAAAA\n');(source/'hg19'/'x.fa.fai').write_text('chr1\t4\n')
            attempt=base/'attempt';attempt.mkdir();dest=private_reference(source,attempt/'reference')
            self.assertFalse((dest/'hg19').is_symlink());self.assertFalse((dest/'hg19'/'x.fa.fai').is_symlink())
            (dest/'hg19'/'x.fa.fai').write_text('private change')
            audit_reference(attempt)
            self.assertEqual((source/'hg19'/'x.fa.fai').read_text(),'chr1\t4\n')

    def test_stage_evidence(self):
        self.assertEqual(infer_stage('Running CNVKit batch\nRunning amplified_intervals\nRunning AA\ncrash'),'aa')

    def test_partial_seed_needs_completed_upstream_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p)
            self.assertTrue(partial_seed_evidence(p)['seed_stage_evaluable'])
            (p/'suite.log').write_text('Running CNVKit batch\nRunning amplified_intervals\ncrash')
            result=partial_seed_evidence(p)
            self.assertFalse(result['seed_stage_evaluable'])
            self.assertNotIn('seed_count',result)

    def test_malformed_cached_result_is_not_reused_or_fatal(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d);job=out/'jobs'/'TEST'/'__whole_sample__'/'v';job.mkdir(parents=True)
            (job/'result.json').write_text('{broken')
            records=run_full_suite({'versions':[{'id':'v'}]},[{'sample_id':'TEST','status':'failed','reason':'fixture invalid input'}],{'v':{'ready':False}},out)
            self.assertEqual(records[0]['status'],'failed')
            self.assertEqual(records[0]['reference_scope'],'complete_reference')

    def test_legacy_offline_policy_requires_exact_audited_commit(self):
        self.assertEqual(offline_policy('543b7ebc8b0d251eefbbcd14e5c1e27635e09a6a',''),'audited_intrinsically_offline_BAM_flow')
        self.assertEqual(offline_policy('unlisted','--no_repo_check'),'explicit_no_repo_check')
        with self.assertRaises(ValueError):offline_policy('unlisted','')

    def test_per_version_solver_is_explicit_and_not_a_silent_fallback(self):
        runtime={'suite_python':'/env/python','suite_script':'/src/Suite.py','aa_python':'/env/python','ac_python':'/env/python','aa_script':'/aa/AmpliconArchitect.py','cnvkit_script':'/env/cnvkit.py'}
        verified={'suite':{'help':'','offline_policy':'audited_intrinsically_offline_BAM_flow'}}
        config={'reference':{'build':'hg19'},'solver':'clarabel'}
        with self.assertRaisesRegex(ValueError,'no substitution'):
            build_command(config,runtime,verified,'/run','/tumor','/normal')
        command=build_command(config,{**runtime,'solver':'mosek'},verified,'/run','/tumor','/normal')
        self.assertNotIn('--AA_solver',command)

    def test_ac_shell_capability_is_verified_before_expensive_pipeline(self):
        with tempfile.TemporaryDirectory() as d:
            source=Path(d);(source/'make_input.sh').write_text('#!/usr/bin/env bash\ndeclare -A graph_lookup\n')
            with self.assertRaisesRegex(ValueError,'Bash associative arrays'):
                require_ac_shell(source,{'bash':{'available':True,'associative_array_probe':{'returncode':1}}})
            require_ac_shell(source,{'bash':{'available':True,'associative_array_probe':{'returncode':0}}})

    def test_ac_shell_failure_is_not_a_pipeline_negative(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);_,v=self.fixture(p)
            (p/'suite.log').write_text('Running AC\nmake_input.sh: line 126: declare: -A: invalid option\nAll stages appear to have completed successfully.\n')
            with self.assertRaises(PipelineStageError) as caught:
                parse_full_result(p,{'sample_id':'TEST'},{'id':'suite-1.6.2'},v)
            self.assertEqual(caught.exception.stage,'ac')
