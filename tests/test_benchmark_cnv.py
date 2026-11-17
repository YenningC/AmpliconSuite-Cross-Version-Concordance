import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ascompare import benchmark_cnv as cnv


class BoundedLocalCnvTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fasta = self.root / 'test.fa'
        self.fasta.write_text('>chr1\nACGT\n')
        Path(str(self.fasta) + '.fai').write_text('chr1\t1000000\t6\t4\t5\n')
        self.background = self.root / 'background.bed'
        self.background.write_text('chr1\t100000\t120000\n')
        self.repo = self.root / 'reference'
        self.repo.mkdir()
        (self.repo / 'file_list.txt').write_text('')
        self.src = self.root / 'aa'
        self.src.mkdir()
        (self.src / 'AmpliconArchitect.py').write_text('# synthetic source\n')
        (self.src / 'amplified_intervals.py').write_text('# synthetic supports --no_cstats\n')
        self.sample = dict(sample_id='SYNTHETIC', source_version='suite-1.3.4', reference_build='hg19', seeds=[
            dict(chrom='chr1', start=10000, end=20000, seed_id='s1')])
        for key, filename in (('alignment', 'tumor.bam'), ('index', 'custom_tumor.bai'),
                              ('normal_alignment', 'normal.bam'), ('normal_index', 'custom_normal.bai')):
            p = self.root / filename
            p.write_bytes(b'SYNTHETIC MOCK ONLY')
            self.sample[key] = str(p)
        self.policy = self.root / 'io_policy.json'
        self.policy.write_text(json.dumps({'alignments': {
            self.sample['alignment']: self.sample['index'],
            self.sample['normal_alignment']: self.sample['normal_index']}}))
        self.config = dict(reference=dict(fasta=str(self.fasta), build='hg19', aa_data_repo=str(self.repo)),
            local_cnv=dict(flank_bp=1000, bin_size=1000, max_window_bp=60000, min_seed_overlap=.8,
                           background_bed=str(self.background)),
            _io_policy_path=str(self.policy), _cnv_cache_dir=str(self.root / 'cache'))
        self.runtime = dict(id='suite-1.6.2', aa_python='/mock/python',
                            aa_script=str(self.src / 'AmpliconArchitect.py'),
                            cnvkit_python='/mock/python', cnvkit_script='/mock/cnvkit.py', cnvkit_version='0.9.13')
        self.identity = dict(version='0.9.13', modules={'cnvlib': {'sha256': 'mock-source'}})
        self.commands = []

    def fake_run(self, python, script, args, policy, cwd, env, stage, timeout=3600):
        self.commands.append((str(script), list(args), stage))
        log = Path(cwd) / (stage + '.log')
        log.write_text('SYNTHETIC MOCK CALLER\n')
        if stage == 'seed_processing':
            prefix = args[args.index('--out') + 1]
            Path(prefix + '.bed').write_text('chr1\t10000\t19000\t8\n')
        elif '-o' in args:
            output = Path(args[args.index('-o') + 1])
            if stage == 'segment':
                output.write_text('chromosome\tstart\tend\tlog2\n'
                                  'chr1\t9000\t21000\t2\nchr1\t100000\t120000\t0\n')
            else:
                output.write_text('chromosome\tstart\tend\tgene\tlog2\tdepth\nchr1\t10000\t20000\t-\t0\t30\n')
        return log

    def run_mock(self, out='result', runtime=None):
        with patch.object(cnv, '_probe_cnvkit', return_value=self.identity), \
             patch.object(cnv, '_run_tool', side_effect=self.fake_run):
            return cnv.run_local_upstream(self.sample, runtime or self.runtime, {}, self.config, self.root / out)

    def test_real_stage_contract_and_retention_endpoint(self):
        result = self.run_mock()[0]
        self.assertEqual(result['status'], 'retained')
        self.assertAlmostEqual(result['retention_fraction'], .9)
        self.assertFalse(result['amplification_confirmed'])
        self.assertEqual([x[2] for x in self.commands],
                         ['coverage_tumor', 'coverage_normal', 'reference', 'fix', 'segment', 'seed_processing'])
        self.assertTrue(all('batch' not in args and 'autobin' not in args and 'access' not in args
                            for _, args, _ in self.commands))
        bed = (self.root / 'result/all_local_cnv_segments.bed').read_text()
        self.assertIn('\t2\n', bed)  # CN-neutral segment was not pre-filtered.
        self.assertIn('--no_cstats', self.commands[-1][1])
        self.assertTrue(result['native_recipe_match'])
        self.assertEqual(result['source_version'], 'suite-1.3.4')
        self.assertEqual(result['version'], 'suite-1.6.2')

    def test_missing_normal_or_background_is_unsupported_before_tools(self):
        self.sample.pop('normal_alignment')
        with patch.object(cnv, '_probe_cnvkit') as probe:
            result = cnv.run_local_upstream(self.sample, self.runtime, {}, self.config, self.root / 'missing')
        self.assertEqual(result[0]['status'], 'unsupported')
        probe.assert_not_called()

    def test_explicit_flat_mode_runs_without_normal_and_is_labelled(self):
        self.config['local_cnv']['reference_mode'] = 'flat'
        self.sample.pop('normal_alignment')
        self.sample.pop('normal_index')
        result = self.run_mock()[0]
        self.assertEqual(result['status'], 'retained')
        self.assertEqual(result['reference_mode'], 'flat')
        self.assertEqual(result['calibration_status'], 'unverified')
        self.assertEqual([x[2] for x in self.commands],
                         ['coverage_tumor', 'reference', 'fix', 'segment', 'seed_processing'])
        ref_args = self.commands[1][1]
        self.assertIn('-t', ref_args)
        self.assertNotIn('normal.targetcoverage.cnn', ' '.join(ref_args))

    def test_reference_strategy_changes_cache_and_flat_needs_background(self):
        first = self.run_mock('normal')[0]
        self.config['local_cnv']['reference_mode'] = 'flat'
        second = self.run_mock('flat')[0]
        self.assertNotEqual(first['cnv_cache_key'], second['cnv_cache_key'])
        self.config['local_cnv'].pop('background_bed')
        self.assertEqual(self.run_mock('no_background')[0]['status'], 'unsupported')

    def test_budget_prevents_any_tool_invocation(self):
        self.config['local_cnv']['max_window_bp'] = 100
        with patch.object(cnv, '_probe_cnvkit') as probe:
            result = cnv.run_local_upstream(self.sample, self.runtime, {}, self.config, self.root / 'over')
        self.assertEqual(result[0]['status'], 'unsupported')
        self.assertIn('budget exceeded', result[0]['reason'])
        probe.assert_not_called()

    def test_different_suite_same_cnv_uses_validated_cache(self):
        first = self.run_mock('v1')[0]
        self.commands.clear()
        second = self.run_mock('v2', dict(self.runtime, id='suite-1.6.1'))[0]
        self.assertFalse(first['cnv_cache_hit'])
        self.assertTrue(second['cnv_cache_hit'])
        self.assertEqual([x[2] for x in self.commands], ['seed_processing'])
        self.assertEqual(first['cnv_cache_key'], second['cnv_cache_key'])

    def test_cache_tamper_and_input_changes_force_recalculation(self):
        self.run_mock('first')
        cached = next((self.root / 'cache').rglob('tumor.cns'))
        cached.write_text('tampered')
        self.commands.clear()
        self.assertFalse(self.run_mock('second')[0]['cnv_cache_hit'])
        Path(self.sample['normal_alignment']).write_bytes(b'changed synthetic input')
        self.assertFalse(self.run_mock('third')[0]['cnv_cache_hit'])

    def test_cache_identity_includes_io_policy_limits(self):
        self.run_mock('first')
        policy = json.loads(self.policy.read_text())
        policy['max_total_fetch_bp'] = 1000
        self.policy.write_text(json.dumps(policy))
        self.assertFalse(self.run_mock('changed_policy')[0]['cnv_cache_hit'])

    def test_calibration_boolean_cannot_create_validation(self):
        self.sample['source_version'] = self.runtime['id']
        self.config['local_cnv']['calibrated'] = True
        self.config['local_cnv']['calibration_status'] = 'verified'
        result = self.run_mock()[0]
        self.assertEqual(result['calibration_status'], 'unverified')
        self.assertIn('provenance', result['calibration_reason'])

    def test_failed_or_missing_filter_is_never_nonretained(self):
        original = self.fake_run
        def failing(*args, **kwargs):
            if args[6] == 'seed_processing':
                raise cnv.LocalFailure('failed', 'seed_processing', 'read budget exhausted')
            return original(*args, **kwargs)
        with patch.object(cnv, '_probe_cnvkit', return_value=self.identity), \
             patch.object(cnv, '_run_tool', side_effect=failing):
            result = cnv.run_local_upstream(self.sample, self.runtime, {}, self.config, self.root / 'failed')[0]
        self.assertEqual(result['status'], 'failed')
        self.assertIsNone(result['retention_fraction'])

    def test_unobserved_gaps_are_not_invented_as_cnv(self):
        source = self.root / 'input.cns'
        source.write_text('chromosome\tstart\tend\tlog2\nchr1\t10\t1000\t2\n')
        result = cnv.cns_to_local_bed(source, [('chr1', 10, 20), ('chr1', 100, 200)], self.root / 'all.bed')
        self.assertEqual(result, [('chr1', 10, 20, 8.0), ('chr1', 100, 200, 8.0)])

    def test_overlap_union_deduplicates_and_threshold_configurable(self):
        self.assertEqual(cnv.retention_fraction(self.sample['seeds'][0],
            [('chr1', 10000, 15000), ('chr1', 12000, 17000)]), .7)
        self.config['local_cnv']['min_seed_overlap'] = .95
        self.assertEqual(self.run_mock()[0]['status'], 'not_retained')

    def test_default_requires_complete_seed_containment(self):
        self.config['local_cnv'].pop('min_seed_overlap')
        result = self.run_mock()[0]
        self.assertEqual(result['min_seed_overlap'], 1.0)
        self.assertEqual(result['status'], 'not_retained')

    def test_native_recipe_mapping_is_not_controlled_0914(self):
        mapped = cnv.load_builtin_map()['versions']
        self.assertEqual(len(mapped), 11)
        self.assertEqual(mapped['suite-1.3.4']['cnvkit_version'], '0.9.10')
        self.assertEqual(mapped['suite-1.5.3']['cnvkit_version'], '0.9.13')
        self.assertEqual(mapped['suite-1.6.2']['cnvkit_version'], '0.9.13')
        result = self.run_mock(runtime=dict(self.runtime, cnvkit_version='0.9.14'))[0]
        self.assertFalse(result['native_recipe_match'])
        self.assertEqual(result['comparison_design'], 'controlled_cnvkit')

    def test_no_seeds_and_cram_are_explicit_nonresults(self):
        self.sample['alignment'] = '/never/read.cram'
        result = cnv.run_local_upstream(self.sample, self.runtime, {}, self.config, self.root / 'cram')
        self.assertEqual(result[0]['status'], 'unsupported')
        self.sample['seeds'] = []
        result = cnv.run_local_upstream(self.sample, self.runtime, {}, self.config, self.root / 'empty')
        self.assertEqual(result[0]['status'], 'no_seeds')

    def test_runner_wraps_python_with_policy_and_retains_error_stage(self):
        with patch.object(cnv, 'guarded_command', return_value=['guarded']) as guard, \
             patch.object(cnv, 'execute', side_effect=TimeoutError('time budget')):
            with self.assertRaises(cnv.LocalFailure) as caught:
                cnv._run_tool('python', 'tool.py', ['coverage'], self.policy,
                              self.root, {}, 'coverage_tumor', 1)
        guard.assert_called_once_with('python', 'tool.py', ['coverage'], self.policy)
        self.assertEqual(caught.exception.stage, 'coverage_tumor')

    def test_public_environment_helper_refuses_existing_environment(self):
        existing = self.root / 'cnvkit-0.9.10-env'
        existing.mkdir()
        marker = existing / 'keep.txt'
        marker.write_text('do not change')
        helper = Path(cnv.__file__).parents[1] / 'scripts/prepare_cnvkit_smoke_envs.py'
        result = subprocess.run([sys.executable, str(helper), '--base-python', sys.executable,
                                 '--out', str(self.root), '--version', '0.9.10'], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Refusing to overwrite', result.stderr)
        self.assertEqual(marker.read_text(), 'do not change')
        self.assertFalse((self.root / 'commands-0.9.10.json').exists())

    def test_calibration_requires_original_input_and_output_identity(self):
        self.sample['source_version'] = self.runtime['id']
        cns, seeds, metadata = (self.root / p for p in ('original.cns', 'original.bed', 'original_run.json'))
        cns.write_text('chromosome\tstart\tend\tlog2\nchr1\t10000\t20000\t2\n')
        seeds.write_text('chr1\t10000\t20000\n')
        original = dict(provenance_kind='original_upstream_run', cnvkit_version='0.9.13', reference_build='hg19',
                        inputs={k: cnv._sampled_identity(self.sample[k]) for k in ('alignment', 'normal_alignment')},
                        outputs={'baseline_cns_sha256': cnv.digest(cns), 'source_seed_bed_sha256': cnv.digest(seeds)})
        metadata.write_text(json.dumps(original))
        self.sample.update(baseline_cns=str(cns), source_seed_bed=str(seeds), baseline_provenance=str(metadata),
                           source_cnvkit_version='0.9.13')
        result = self.run_mock('matched')[0]
        self.assertEqual(result['calibration_status'], 'passed')
        self.assertEqual(result['calibration_role'], 'source_control')
        self.assertEqual(result['calibration_scope'], 'selected_source_candidate_retention_only')
        candidate = dict(self.runtime, id='suite-1.6.1')
        observation = self.run_mock('candidate', candidate)[0]
        self.assertEqual(observation['status'], 'retained')
        self.assertEqual(observation['calibration_status'], 'unverified')
        self.assertEqual(observation['calibration_role'], 'candidate_observation')
        self.assertFalse((self.root / 'candidate/calibration_evidence.json').exists())
        original['inputs']['alignment']['bytes'] += 1
        metadata.write_text(json.dumps(original))
        self.assertEqual(self.run_mock('mismatch')[0]['calibration_status'], 'unverified')

    @unittest.skipUnless(os.environ.get('CNVKIT_SMOKE_PYTHON'), 'Optional actual CNVkit synthetic integration test')
    def test_actual_cnvkit_finite_synthetic_bam(self):
        import pysam
        if os.environ.get('CNVKIT_SMOKE_ARTIFACTS'):
            self.addCleanup(shutil.copytree, self.root, os.environ['CNVKIT_SMOKE_ARTIFACTS'])
        # This test creates its own tiny reference and BAMs. No patient path is read.
        self.fasta.write_text('>chr1\n' + 'ACGT' * 250000 + '\n')
        Path(str(self.fasta) + '.fai').unlink()
        pysam.faidx(str(self.fasta))
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 1000000}]}
        for key, index_key, tumor in (('alignment', 'index', True), ('normal_alignment', 'normal_index', False)):
            path = self.sample[key]
            with pysam.AlignmentFile(path, 'wb', header=header) as handle:
                n = 0
                for start, end, step in ((9000, 21000, 10 if tumor else 30), (100000, 120000, 30)):
                    for position in range(start, end - 150, step):
                        read = pysam.AlignedSegment()
                        read.query_name = f'SYNTHETIC_{n}'
                        read.query_sequence = 'A' * 150
                        read.flag = 0
                        read.reference_id = 0
                        read.reference_start = position
                        read.mapping_quality = 60
                        read.cigar = [(0, 150)]
                        read.query_qualities = pysam.qualitystring_to_array('I' * 150)
                        handle.write(read)
                        n += 1
            pysam.index(path)
            self.sample[index_key] = path + '.bai'
        self.policy.write_text(json.dumps({'alignments': {self.sample['alignment']: self.sample['index'],
                                                         self.sample['normal_alignment']: self.sample['normal_index']}}))
        python = os.environ['CNVKIT_SMOKE_PYTHON']
        self.config['local_cnv']['reference_mode'] = os.environ.get('CNVKIT_SMOKE_REFERENCE_MODE', 'matched_normal')
        if self.config['local_cnv']['reference_mode'] == 'flat':
            self.sample.pop('normal_alignment')
            self.sample.pop('normal_index')
        runtime = dict(self.runtime, cnvkit_python=python, cnvkit_script=str(Path(python).parent / 'cnvkit.py'),
                       cnvkit_version=os.environ.get('CNVKIT_SMOKE_VERSION', '0.9.14'),
                       extra_path=[p for p in os.environ.get('CNVKIT_SMOKE_EXTRA_PATH', '').split(os.pathsep) if p])
        out = self.root / 'actual_cnv_smoke'
        out.mkdir()
        env = cnv._environment(runtime, python)
        windows = cnv.prepare_windows(self.sample['seeds'], self.fasta, self.config['local_cnv'], out)
        identity = cnv._probe_cnvkit(runtime, self.policy, out, env)
        try:
            cns, cache = cnv._cached_cnv(self.sample, runtime, self.config, windows, identity, out, self.policy, env)
        except Exception:
            for log in (self.root / 'cache').rglob('*.log'):
                print(log.name, log.read_text()[-8000:])
            raise
        self.assertTrue(cns.is_file())
        self.assertGreater(len(cnv.cns_to_local_bed(cns, windows['intervals'], out / 'cnv.bed')), 0)
        for audit in (self.root / 'cache').rglob('bam_access.*.json'):
            state = json.loads(audit.read_text())
            self.assertFalse(state['blocked'])
            self.assertLessEqual(state['requested_bp'], 100000)
        audits = [json.loads(p.read_text()) for p in self.root.rglob('bam_access.*.json')]
        if self.config['local_cnv']['reference_mode'] == 'flat':
            self.assertFalse(any(Path(path).name == 'normal.bam' for state in audits for path in state['opened']))
        summary = {'label': 'SYNTHETIC_ONLY', 'cnvkit_version': identity['version'],
                   'reference_mode': self.config['local_cnv']['reference_mode'],
                   'guarded_processes': len(audits), 'blocked': sum(len(x['blocked']) for x in audits),
                   'queries_per_process': [x['queries'] for x in audits],
                   'requested_bp_per_process': [x['requested_bp'] for x in audits],
                   'native_bedcov_read_count_observable': False}
        (out / 'smoke_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        print(json.dumps(summary))


if __name__ == '__main__':
    unittest.main()
