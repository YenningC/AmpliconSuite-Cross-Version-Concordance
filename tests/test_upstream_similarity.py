import unittest

from ascompare.concordance_report import _upstream
from ascompare.concordance_upstream_similarity import upstream_similarity_cohort_rows, upstream_similarity_rows


class UpstreamSimilarityTests(unittest.TestCase):
    def _rows(self, calls, panel, omit=(), edit=None, declare=True):
        """Two synthetic Suite releases, overlapping original candidate intervals."""
        samples = {'synthetic': {'source_version': 'source'}}
        if declare:
            samples['synthetic']['upstream_seed_ids'] = list(panel)
        records = []
        for version, intervals in calls.items():
            for seed, (start, end) in panel.items():
                if (version, seed) in omit:
                    continue
                record = dict(sample_id='synthetic', source_version='source', version=version, seed_id=seed,
                              seed=dict(chrom='chr1', start=start, end=end), cnvkit_version='cnv-'+version,
                              status='not_retained', comparison_eligible=True, reference_mode='flat',
                              matching_rule='overlap >= minimum', min_seed_overlap=1.,
                              retained_intervals=intervals)
                if edit:
                    edit(record)
                records.append(record)
        seeds, samples, _ = _upstream(records, samples, ['target', 'source'], 'target', .05, .025, 2, None)
        pairs = upstream_similarity_rows(seeds, samples, 'target')
        pair = next(p for p in pairs if p['target_version'] == 'target' and p['comparator_version'] == 'source')
        return pair, samples

    def test_nonempty_jaccard_uses_base_pairs_and_unions_overlapping_original_candidates(self):
        panel = {'a': [0, 10], 'b': [5, 15]}
        # Calls beyond the input panel and other chromosomes must not contribute.
        pair, samples = self._rows({'target': [['chr1', 4, 18], ['chrX', 0, 100]],
                                    'source': [['chr1', 2, 6], ['chr1', 2, 6]]}, panel,
                                   edit=lambda r: r.update(status='retained') if
                                   (r['version'], r['seed_id']) == ('target', 'b') else None)
        self.assertEqual(pair['upstream_similarity_status'], 'ok')
        self.assertEqual(pair['upstream_intersection_bp'], 2)
        self.assertEqual(pair['upstream_union_bp'], 13)
        self.assertAlmostEqual(pair['upstream_interval_jaccard'], 2/13.)
        self.assertEqual(pair['upstream_retained_count'], 1)
        self.assertEqual(pair['upstream_retained_ratio'], .5)
        self.assertEqual(samples[0]['upstream_retained_ratio'], .5)
        self.assertEqual(pair['upstream_panel_retained_ratio'], .5)
        self.assertEqual(pair['upstream_unevaluable_count'], 0)
        self.assertEqual(pair['upstream_target_cnvkit_versions'], ['cnv-target'])

    def test_mutual_absence_is_undefined_and_one_sided_signal_is_zero(self):
        panel = {'a': [1, 10]}
        row, _ = self._rows({'source': [], 'target': []}, panel)
        self.assertIsNone(row['upstream_interval_jaccard'])
        self.assertEqual(row['upstream_similarity_status'], 'both_empty')
        self.assertEqual(row['upstream_union_bp'], 0)
        row, _ = self._rows({'source': [['chr1', 5, 6]], 'target': []}, panel)
        self.assertEqual(row['upstream_interval_jaccard'], 0)

    def test_zero_candidates_and_missing_rows_do_not_imply_empty_calls(self):
        pair, _ = self._rows({}, {})
        self.assertEqual(pair['upstream_similarity_status'], 'no_candidates')
        self.assertIsNone(pair['upstream_retained_ratio'])
        pair, _ = self._rows({'source': []}, {'a': [1, 9]})
        self.assertEqual(pair['upstream_similarity_status'], 'unavailable_seed_result')
        self.assertEqual(pair['upstream_count_status'], 'unavailable')
        self.assertIsNone(pair['upstream_panel_retained_ratio'])
        self.assertEqual(pair['upstream_unevaluable_count'], 1)
        self.assertIsNone(pair['upstream_intersection_bp'])

    def test_partial_pairs_have_observed_only_score_and_evaluable_count_denominator(self):
        pair, _ = self._rows({'source': [['chr1', 1, 2]], 'target': [['chr1', 1, 2]]},
                             {'a': [0, 2], 'b': [3, 7]}, omit=[('target', 'b')])
        self.assertIsNone(pair['upstream_interval_jaccard'])
        self.assertEqual(pair['upstream_observed_interval_jaccard'], 1.)
        self.assertEqual(pair['upstream_similarity_status'], 'partial_candidate_evidence')
        self.assertEqual(pair['upstream_count_status'], 'partial')
        self.assertEqual(pair['upstream_evaluable_count'], 1)
        self.assertEqual(pair['upstream_expected_seeds'], 2)
        self.assertIsNone(pair['upstream_panel_retained_ratio'])
        self.assertEqual(pair['upstream_unevaluable_count'], 1)
        self.assertEqual(pair['upstream_retained_ratio'], 0.)

    def test_missing_suite_bed_evidence_and_invalid_bed_are_not_a_negative(self):
        for intervals in [None, [['chr1', 9, 2]], [['chr1', 1.8, 2]], [1]]:
            with self.subTest(intervals=intervals):
                pair, _ = self._rows({'source': [], 'target': intervals}, {'a': [1, 9]})
                self.assertEqual(pair['upstream_similarity_status'], 'missing_or_invalid_final_seed_coordinates')
                self.assertIsNone(pair['upstream_interval_jaccard'])
                self.assertIsNone(pair['upstream_target_bp'])

    def test_nonmatching_context_undeclared_universe_and_ineligible_positive(self):
        pair, _ = self._rows({'source': [], 'target': []}, {'a': [1, 9]},
                             edit=lambda r: r.update(seed=dict(chrom='chr1', start=2, end=9))
                             if r['version'] == 'target' else None)
        self.assertEqual(pair['upstream_similarity_status'], 'mismatched_candidate_context')
        pair, _ = self._rows({'source': [], 'target': []}, {'a': [1, 9]}, declare=False)
        self.assertEqual(pair['upstream_similarity_status'], 'candidate_universe_undeclared')
        self.assertIsNone(pair['upstream_interval_jaccard'])
        pair, _ = self._rows({'source': [], 'target': []}, {'a': [1, 9]},
                             edit=lambda r: r.update(status='retained', comparison_eligible=False))
        self.assertEqual(pair['upstream_retained_count'], 0)
        self.assertIsNone(pair['upstream_retained_ratio'])

    def test_cohort_weights_patients_equally_and_keeps_count_denominators(self):
        row, _ = self._rows({'source': [['chr1', 1, 2]], 'target': [['chr1', 1, 2]]}, {'a': [0, 4]})
        def case(sample, score, n, retained, similarity_status='ok', count_status='complete'):
            return dict(row, sample_id=sample, upstream_interval_jaccard=score,
                        upstream_similarity_status=similarity_status, upstream_count_status=count_status,
                        upstream_retained_count=retained, upstream_evaluable_count=n,
                        upstream_expected_seeds=n + int(count_status == 'partial'),
                        upstream_unevaluable_count=int(count_status == 'partial'),
                        upstream_retained_ratio=retained/n)
        rows = [case('a1', 1., 20, 20), case('a2', 1., 15, 15), case('b', 0., 1, 0),
                case('c', None, 1, 0, 'partial_candidate_evidence', 'partial')]
        patient_ids = {'a1': 'A', 'a2': 'A', 'b': 'B', 'c': 'C'}
        cohort, = upstream_similarity_cohort_rows(rows, patient_ids)
        self.assertEqual(cohort['upstream_interval_jaccard'], .5)
        self.assertEqual(cohort['upstream_retained_ratio'], .5)
        self.assertEqual(cohort['upstream_interval_jaccard_n_samples'], 3)
        self.assertEqual(cohort['upstream_interval_jaccard_n_independent_units'], 2)
        self.assertEqual(cohort['n_selected_independent_units'], 3)
        self.assertEqual(cohort['upstream_retained_count'], 35)
        self.assertEqual(cohort['upstream_evaluable_count'], 37)
        self.assertEqual(cohort['upstream_expected_seeds'], 38)
        self.assertEqual(cohort['upstream_unevaluable_count'], 1)
        self.assertEqual(cohort['upstream_similarity_status_counts']['partial_candidate_evidence'], 1)
        self.assertEqual(cohort['upstream_interval_jaccard_ci_low'], 0.)
        self.assertEqual(cohort['upstream_interval_jaccard_ci_high'], 1.)
        by_sample, = upstream_similarity_cohort_rows(rows, None)
        self.assertAlmostEqual(by_sample['upstream_retained_ratio'], 2/3.)
        identity, = upstream_similarity_cohort_rows([dict(rows[0], comparator_version='target')], None)
        self.assertIsNone(identity['upstream_interval_jaccard_ci_high'])


if __name__ == '__main__':
    unittest.main()
