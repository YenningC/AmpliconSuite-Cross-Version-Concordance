"""Timing records describe execution cost, never statistical compatibility."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ascompare.concordance_timing import (
    TimingSession, measure, native_aa_intervals, record_native_aa, timed_call,
)


class ConcordanceTimingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def rows(self):
        with (self.root / 'stage_timing.tsv').open() as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def log(self, lines):
        path = self.root / 'native.log'
        path.write_text('\n'.join(lines) + '\n')
        return path

    def test_nested_wall_subprocess_and_path_metadata_are_saved_without_double_counting(self):
        with patch('ascompare.concordance_timing.time.perf_counter', side_effect=[0, 1, 2, 5, 7, 10]):
            with TimingSession(self.root):
                with measure('local_upstream', sample_id='SYNTHETIC', version='suite-test'):
                    with measure('segment', measurement_kind='subprocess', log_path=self.root / 'segment.log'):
                        pass
        rows = {row['stage']: row for row in self.rows()}
        self.assertEqual(float(rows['workflow']['duration_seconds']), 10)
        self.assertEqual(float(rows['local_upstream']['duration_seconds']), 6)
        self.assertEqual(float(rows['segment']['duration_seconds']), 3)
        self.assertEqual(rows['segment']['sample_id'], 'SYNTHETIC')
        self.assertEqual(rows['segment']['version'], 'suite-test')
        self.assertEqual(rows['segment']['log_path'], str(self.root / 'segment.log'))
        self.assertEqual(rows['segment']['parent_event_id'], rows['local_upstream']['event_id'])
        self.assertEqual([r['stage'] for r in rows.values() if r['accounting'] == 'invocation_total'], ['workflow'])
        self.assertEqual(rows['segment']['measurement_kind'], 'subprocess')

    def test_failures_and_interruptions_keep_actual_elapsed_time_and_resume_history(self):
        for error, status in ((RuntimeError('synthetic failure'), 'failed'), (KeyboardInterrupt(), 'interrupted')):
            with self.assertRaises(type(error)):
                with patch('ascompare.concordance_timing.time.perf_counter', side_effect=[0, 1, 4, 5]):
                    with TimingSession(self.root):
                        with measure('aa', measurement_kind='subprocess'):
                            raise error
            rows = self.rows()[-2:]
            self.assertTrue(all(r['status'] == status for r in rows))
            self.assertEqual(float(rows[0]['duration_seconds']), 3)
            self.assertEqual(float(rows[1]['duration_seconds']), 5)
        self.assertEqual(len({r['invocation_id'] for r in self.rows()}), 2)

    def test_native_markers_use_qualified_intervals_and_prefer_coverage_counters(self):
        path = self.log([
            '#TIME 1 Initiating bam_to_breakpoint object',
            '#TIME 3 Exploring interval: synthetic',
            '#TIME 4 Fetching coverage for synthetic',
            '#TIME 5 Obtained coverage list for synthetic',
            '#TIME 6 Interval sets for amplicons determined:',
            '#TIME 7 Reconstructing amplicon1',
            '#TIME 11 Completed decomposition',
            '#TIME 12 Plotting SV View for amplicon1',
            '#TIME 15 Reconstructing amplicon2',
            '#TIME 17 Completed decomposition',
            '#TIME 18 Plotting SV View for amplicon2',
            '#TIME 20 Total Runtime',
            'Time spent counting coverage (slow, fast): 0.250 s, 0.750 s',
        ])
        rows = native_aa_intervals(path)
        by_stage = {row['stage']: row for row in rows}
        self.assertEqual(by_stage['aa_background_and_input_setup']['duration_seconds'], 2)
        self.assertTrue(by_stage['aa_background_and_input_setup']['approximate'])
        self.assertEqual(by_stage['aa_seed_extension']['duration_seconds'], 3)
        self.assertFalse(by_stage['aa_seed_extension']['approximate'])
        self.assertEqual([r['duration_seconds'] for r in rows if r['stage'] == 'aa_reconstruction'], [4, 2])
        self.assertEqual([r['duration_seconds'] for r in rows if r['stage'] == 'aa_plotting'], [3, 2])
        self.assertTrue(all(r['approximate'] for r in rows if r['stage'] == 'aa_plotting'))
        self.assertEqual(by_stage['aa_coverage_counting']['duration_seconds'], 1)
        self.assertEqual(by_stage['aa_coverage_counting']['method'], 'native cumulative coverage counters')

    def test_large_coverage_debug_log_has_constant_size_timing_output(self):
        path = self.root / 'native.log'
        with path.open('w') as handle:
            for number in range(20_000):
                handle.write(f'#TIME {number * 2} Fetching coverage for interval{number}\n')
                handle.write(f'#TIME {number * 2 + 1} Obtained coverage list for interval{number}\n')
        rows = native_aa_intervals(path)
        self.assertLessEqual(len(rows), 5)
        coverage = next(r for r in rows if r['stage'] == 'aa_coverage_counting')
        self.assertEqual(coverage['duration_seconds'], 20_000)
        self.assertIn('completed_queries=20000', coverage['notes'])

    def test_unfinished_and_nonmonotonic_native_phases_are_not_invented(self):
        rows = native_aa_intervals(self.log(['#TIME 1 Reconstructing amplicon1', '#TIME 2 Fetching coverage for x']))
        self.assertTrue(all(r['duration_seconds'] is None for r in rows))
        rows = native_aa_intervals(self.log([
            '#TIME 3 Reconstructing amplicon1', '#TIME 2 Completed decomposition',
            'Time spent counting coverage (slow, fast): 0.100 s, 0.200 s',
        ]))
        self.assertTrue(all(r['duration_seconds'] is None for r in rows if r['stage'] != 'aa_coverage_counting'))
        # The independently recorded cumulative counter is still usable.
        self.assertAlmostEqual(next(r['duration_seconds'] for r in rows if r['stage'] == 'aa_coverage_counting'), .3)

    def test_cache_lookup_and_producer_native_time_are_separate(self):
        run = self.root / 'attempt'
        (run / 'aa').mkdir(parents=True)
        (run / 'aa/reconstruction.log').write_text('#TIME 1 Reconstructing amplicon1\n#TIME 101 Completed decomposition\n')
        with patch('ascompare.concordance_timing.time.perf_counter', side_effect=[0, 1, 2, 3]):
            with TimingSession(self.root):
                with measure('aa_cache', cache_status='hit', producer_id=run):
                    pass
                record_native_aa(run, 'SYNTHETIC', 'suite-test', role='producer_reference', producer_id=str(run))
        rows = self.rows()
        lookup = next(r for r in rows if r['stage'] == 'aa_cache')
        native = next(r for r in rows if r['stage'] == 'aa_reconstruction')
        self.assertEqual(float(lookup['duration_seconds']), 1)
        self.assertEqual(float(native['duration_seconds']), 100)
        self.assertEqual(native['accounting'], 'producer_reference')
        self.assertEqual(native['cache_status'], 'hit')
        self.assertEqual(native['producer_id'], str(run))
        self.assertIn('excluded from current execution cost', native['notes'])

    def test_missing_native_log_is_explicitly_unavailable(self):
        with TimingSession(self.root):
            record_native_aa(self.root / 'missing_attempt', 'SYNTHETIC', 'suite-test')
        row = next(r for r in self.rows() if r['stage'] == 'aa_native_phases')
        self.assertEqual(row['status'], 'unavailable')
        self.assertEqual(row['duration_seconds'], '')
        self.assertIn('absent', row['notes'])

    def test_resume_after_partial_timing_event_preserves_new_valid_events(self):
        (self.root / '.stage_timing_events.jsonl').write_text('{"incomplete":')
        with self.assertWarnsRegex(UserWarning, 'incomplete timing event'):
            with TimingSession(self.root):
                with measure('input_check'):
                    pass
        self.assertEqual({r['stage'] for r in self.rows()}, {'input_check', 'workflow'})

    def test_cnv_failure_and_timeout_are_measured_as_failed_subprocesses(self):
        from ascompare.concordance_cnv import LocalFailure, _run_tool
        for outcome in (4, TimeoutError('synthetic timeout')):
            def execute(*args, **kwargs):
                Path(args[3]).write_text('SYNTHETIC failure')
                if isinstance(outcome, Exception):
                    raise outcome
                return outcome
            with patch('ascompare.concordance_cnv.guarded_command', return_value=['synthetic']), \
                 patch('ascompare.concordance_cnv.execute', side_effect=execute), \
                 patch('ascompare.concordance_timing.time.perf_counter', side_effect=[0, 1, 4, 5]):
                with TimingSession(self.root):
                    with self.assertRaises(LocalFailure):
                        _run_tool('python', 'synthetic.py', [], self.root / 'policy', self.root, {}, 'segment')
        segment = [r for r in self.rows() if r['stage'] == 'segment']
        self.assertEqual(len(segment), 2)
        self.assertTrue(all(r['status'] == 'failed' and float(r['duration_seconds']) == 3 for r in segment))

    def test_result_cache_metadata_describes_current_call_only(self):
        with TimingSession(self.root):
            result = timed_call('ac_reclassification', lambda: {'status': 'ok', 'cache_hit': True, 'run_dir': str(self.root / 'producer')})
        self.assertTrue(result['cache_hit'])
        row = next(r for r in self.rows() if r['stage'] == 'ac_reclassification')
        self.assertEqual(row['cache_status'], 'hit')
        self.assertEqual(row['producer_id'], str(self.root / 'producer'))
        self.assertIn('not repeated', row['notes'])

    def test_normal_cli_demo_and_resume_write_per_invocation_stage_timing(self):
        from ascompare.concordance import main
        out = self.root / 'workflow'
        with patch('ascompare.concordance_report._plots'):
            self.assertEqual(main(['demo', '--out', str(out)]), 0)
            provenance = json.loads((out / 'demo_provenance.json').read_text())
            self.assertEqual(main(['run', '--manifest', provenance['manifest'], '--out', str(out)]), 0)
        with (out / 'stage_timing.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        invocations = {r['invocation_id'] for r in rows}
        self.assertEqual(len(invocations), 2)
        for invocation in invocations:
            local = [r for r in rows if r['invocation_id'] == invocation]
            for stage, count in (('workflow', 1), ('input_check', 3), ('result_import', 6), ('report_generation', 1)):
                self.assertEqual(sum(r['stage'] == stage for r in local), count)
            self.assertFalse(any(r['measurement_kind'] == 'subprocess' for r in local))


if __name__ == '__main__':
    unittest.main()
