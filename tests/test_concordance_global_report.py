"""Synthetic whole-genome seed comparisons; no patient data or native executables."""
import tempfile
import unittest
from pathlib import Path

from ascompare.concordance_global_report import global_loci, global_report_units, GLOBAL_METRIC
from ascompare.concordance_report import write_concordance_reports


def seed(a, b):
    return dict(chrom='chr1', start=a, end=b, seed_id=f'native_{a}_{b}')


def call(*spans, status='ok'):
    return dict(status=status, seeds=[seed(a, b) for a, b in spans], cnvkit_version='test',
                reference_mode='matched_normal', supported_labels=['ecDNA'], native_recipe_match=True)


def sample(calls):
    ids = [row['seed_id'] for row in global_loci(calls)]
    return dict(sample_id='SYNTHETIC', source_version='S', global_seed_sets=calls,
                selected_seed_ids=ids, upstream_seed_ids=ids, global_source_seed_calibration='passed',
                source_bridge_status='passed')


def report(folder, spec, joints):
    targets = list(spec['global_seed_sets'])
    up, down = global_report_units(spec, joints, targets)
    cfg = dict(mode='with_bam', cnv_scope='global', synthetic=True, reference_version='T',
               samples=[spec], target_versions=targets, statistics=dict(bootstrap=100, permutations=100))
    return write_concordance_reports(folder, down, up, cfg)


class GlobalReportTests(unittest.TestCase):
    def test_pair_counts_and_sampled_downstream(self):
        spec = sample({'S': call((1, 6), (7, 10), (20, 23)), 'T': call((3, 8), (30, 50)),
                       'OTHER': call((80, 90)), 'FAIL': call(status='failed')})
        self.assertEqual(spec['upstream_seed_ids'], ['chr1:1-10', 'chr1:20-23', 'chr1:30-50', 'chr1:80-90'])
        spec['selected_seed_ids'] = ['chr1:1-10', 'chr1:30-50']
        joint = dict(sample_id='SYNTHETIC', version='S', status='ok', comparison_eligible=True,
                     labels=['ecDNA'], supported_labels=['ecDNA'], aa_amplicon_count=1,
                     coverage_status='complete', amplicons=[dict(amplicon_id=1, labels=['ecDNA'],
                     reconstruction_intervals=[('chr1', 1, 100)])], features=[dict(amplicon_id=1,
                     label='ecDNA', intervals=[('chr1', 1, 100)])])
        up, down = global_report_units(spec, {'S':joint}, list(spec['global_seed_sets']))
        spec['global_seed_sets']['T']['historical_provenance_limited'] = True
        self.assertEqual(len(up), 16)
        self.assertEqual(len(down), 8)
        self.assertTrue(all(r['detected'] for r in down if r['version'] == 'S'))
        self.assertEqual({r['status'] for r in down if r['version'] == 'OTHER'}, {'not_detected'})
        self.assertEqual({r['status'] for r in down if r['version'] == 'T'}, {'missing'})
        self.assertEqual({r['status'] for r in down if r['version'] == 'FAIL'}, {'failed'})
        with tempfile.TemporaryDirectory() as folder:
            summary = report(folder, spec, {'S':joint})
            stats = next(r for r in summary['upstream'] if r['target_version'] == 'T')
            self.assertEqual(stats['metric'], GLOBAL_METRIC)
            self.assertEqual(stats['estimate'], 1.)
            self.assertTrue(stats['calibration_passed'])
            result = next(r for r in summary['result'] if r['target_version']=='T' and r['comparator_version']=='S')
            self.assertEqual(result['target_native_seed_count'], 2)
            self.assertEqual(result['comparator_native_seed_count'], 3)
            self.assertEqual(result['comparator_native_seeds_retained_count'], 0)
            self.assertEqual(result['comparator_native_seeds_retained_ratio'], 0.)
            combined = next(r for r in summary['compatibility_combined'] if r['target_version'] == 'T')
            self.assertIn('historical_global_cnv_computation_provenance_limited', combined['reasons'])
            self.assertEqual(result['upstream_extra_loci'], 1)
            self.assertEqual(result['upstream_missed_loci'], 1)
            self.assertAlmostEqual(result['locus_discordance'], 2/3)
            self.assertEqual(result['upstream_pair_positive_loci'], 3)
            self.assertTrue((Path(folder) / 'upstream_retention_ladder.png').is_file())
            self.assertNotIn('局部重跑', (Path(folder) / 'report.md').read_text())
            self.assertNotIn('局部对照', (Path(folder) / 'result.md').read_text())

    def test_complete_empty_and_failure_do_not_establish_agreement(self):
        spec = sample({'S':call(), 'T':call(), 'F':call(status='failed')})
        with tempfile.TemporaryDirectory() as folder:
            summary = report(folder, spec, {})
            rows = {r['target_version']:r for r in summary['result'] if r['comparator_version'] == 'S'}
            self.assertIsNone(rows['T']['upstream_interval_jaccard'])
            self.assertEqual(rows['T']['upstream_similarity_status'], 'both_empty')
            self.assertIsNone(rows['F']['target_native_seed_count'])
            stats = {r['target_version']:r for r in summary['upstream']}
            self.assertEqual(stats['F']['no_seed_samples'], 0)
            self.assertEqual(stats['T']['no_seed_samples'], 1)
            self.assertEqual(rows['F']['upstream_similarity_status'], 'unavailable')
            self.assertTrue(all(r['decision']=='insufficient_evidence' for r in summary['compatibility_combined']))

    def test_native_retention_multiple_pieces_and_threshold(self):
        from ascompare.concordance_upstream_similarity import _native_retention
        source = {'l': {'native_seeds': [seed(10, 20), seed(0, 5)]}}
        target = {'l': {'retained_intervals': [['chr1', 10, 13], ['chr1', 15, 20]]}}
        self.assertEqual(_native_retention(source, target, .8), (1, .5))
        self.assertEqual(_native_retention(source, target, 1.), (0, 0.))

    def test_bridge_version_cannot_hide_nonoverlap_or_empty_native_set(self):
        spec = sample({'S':call((0, 5)), 'T':call((100, 200)), 'BRIDGE':call((0, 200)), 'EMPTY':call()})
        with tempfile.TemporaryDirectory() as folder:
            summary = report(folder, spec, {})
            rows = {r['target_version']:r for r in summary['result'] if r['comparator_version'] == 'S'}
            self.assertEqual(rows['T']['locus_discordance'], 0.)
            self.assertEqual(rows['T'][GLOBAL_METRIC], 1.)
            self.assertEqual(rows['EMPTY'][GLOBAL_METRIC], 1.)
            self.assertIsNone(rows['EMPTY']['target_native_seeds_retained_ratio'])
            self.assertEqual(next(r for r in summary['upstream'] if r['target_version']=='T')['estimate'], 1.)
            self.assertEqual(rows['BRIDGE']['comparator_native_seeds_retained_ratio'], 1.)
            self.assertEqual(rows['BRIDGE']['target_native_seeds_retained_ratio'], 0.)
            self.assertEqual(rows['BRIDGE'][GLOBAL_METRIC], 1.)

    def test_adjacent_are_distinct(self):
        self.assertEqual(len(global_loci({'A':call((0, 9), (9, 10))})), 2)
