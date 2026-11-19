"""Synthetic paired CNVkit controls; no BAM access or Suite component execution."""
import csv
import tempfile
import unittest
from pathlib import Path

from ascompare.concordance_local_cnv_report import write_local_cnv_comparator

SRC = 'suite-1.6.1'
REF = 'suite-1.6.2'
OLD = 'suite-1.3.3'
BASE = '0.9.14'
TARGET = '0.9.10'
PANEL = {'a':(0, 10), 'b':(5, 15)}


def sample(sid='SYNTHETIC', source=BASE):
    return dict(sample_id=sid, source_version=SRC, source_cnvkit_version=source, reference_build='hg19',
                upstream_seed_ids=list(PANEL))


def observation(version, seed, intervals, sid='SYNTHETIC'):
    a, b = PANEL[seed]
    covered = {x for chrom, start, end in intervals if chrom == 'chr1' for x in range(max(start,a), min(end,b))}
    ratio = len(covered)/(b-a)
    baseline = version == BASE
    return dict(sample_id=sid, source_version=SRC, version=SRC, cnvkit_version=version,
                seed_id=seed, seed=dict(chrom='chr1', start=a, end=b), matching_rule='fraction >= 1',
                min_seed_overlap=1., reference_mode='flat', retained_intervals=intervals,
                comparison_eligible=True, status='retained' if ratio == 1 else 'not_retained',
                source_cnvkit_version=BASE, actual_cnvkit_version=version,
                source_filter_suite_version=SRC, source_filter_aa_version='1.6.r1', local_design_id=sid+'-D',
                measurement_role='source_local_standard' if baseline else 'comparison',
                measurement_origin='controlled_local_job', calibration_role='source_control' if baseline else 'candidate_observation',
                calibration_status='unverified', calibration_reason='synthetic historical provenance missing')


def collected(version, intervals, sid='SYNTHETIC'):
    return [observation(version, seed, intervals, sid) for seed in PANEL]


class LocalCnvReportTests(unittest.TestCase):
    def render(self, records, samples=None, mode='with_bam', patients=None):
        samples = [sample()] if samples is None else samples
        config = dict(mode=mode, reference_version=REF, target_versions=[OLD, REF, SRC],
                      samples=samples, runtimes=[], reference=dict(build='hg19'), local_cnv=dict(reference_mode='flat'))
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            result = write_local_cnv_comparator(out, records, config, patients, .05, True)
            with (out/result['files']['seeds']).open() as f:
                seeds = list(csv.DictReader(f, delimiter='\t'))
        return result['samples'], result['cohort'], seeds

    def test_distinct_input_version_anchor_and_paired_effects(self):
        samples, cohort, seeds = self.render(collected(BASE, [['chr1',0,10]])+collected(TARGET,[['chr1',5,15]]))
        by_version = {r['target_cnvkit_version']:r for r in samples}
        baseline = by_version[BASE]
        self.assertTrue(baseline['identity_anchor'])
        self.assertEqual(baseline['target_suite_versions'], [])
        row = by_version[TARGET]
        self.assertEqual(row['target_suite_versions'], [OLD])
        self.assertEqual((row['extra_count'],row['missed_count'],row['retention_delta']), (1,1,0))
        self.assertEqual((row['intersection_bp'],row['union_bp']), (5,15))
        self.assertAlmostEqual(row['interval_jaccard'], 1/3)
        self.assertEqual(row['baseline_panel_retention'], .5)
        self.assertEqual(row['baseline_original_calibration_status'], 'unverified')
        self.assertEqual(row['status'], 'ok')
        self.assertEqual([int(r['retention_delta']) for r in seeds if r['target_cnvkit_version']==TARGET], [-1,1])
        identity = next(x for x in cohort if x['identity_anchor'])
        self.assertIsNone(identity['interval_jaccard_ci_high'])
        self.assertIsNotNone(identity['baseline_panel_retention_ci_high'])

    def test_missing_target_keeps_baseline_and_partial_delta_is_only_descriptive(self):
        records=collected(BASE,[['chr1',0,15]])+[observation(TARGET,'a',[])]
        for r in records[:2]:
            r.update(calibration_status='passed')
        samples, _, _ = self.render(records)
        row = next(r for r in samples if r['target_cnvkit_version']==TARGET)
        self.assertEqual(row['baseline_original_calibration_status'],'passed')
        self.assertEqual(row['status'],'partial')
        self.assertEqual((row['paired_seed_count'],row['expected_seed_count']), (1,2))
        self.assertIsNone(row['retention_delta'])
        self.assertEqual(row['observed_retention_delta'],-1.)
        self.assertIsNone(row['interval_jaccard'])
        self.assertEqual(row['observed_interval_jaccard'],0.)
        self.assertEqual((row['target_evaluable_count'],row['baseline_retained_count']),(1,2))
        missing=next(r for r in samples if r['target_cnvkit_version']=='0.9.13')
        self.assertEqual(missing['status'],'target_unavailable')
        self.assertIsNone(missing['target_retained_count'])
        self.assertEqual(missing['baseline_panel_retention'],1.)

    def test_missing_baseline_and_both_empty_are_distinct(self):
        target=collected(TARGET,[])
        samples, _, _=self.render(target)
        row=next(r for r in samples if r['target_cnvkit_version']==TARGET)
        self.assertEqual(row['status'],'source_baseline_unavailable')
        self.assertEqual(row['target_retained_count'],0)
        self.assertIsNone(row['baseline_retained_count'])
        self.assertIsNone(row['interval_jaccard'])
        samples,_,_=self.render(target+collected(BASE,[]))
        row=next(r for r in samples if r['target_cnvkit_version']==TARGET)
        self.assertEqual(row['status'],'both_empty')
        self.assertEqual(row['retention_delta'],0.)
        self.assertEqual(row['union_bp'],0)
        self.assertIsNone(row['interval_jaccard'])

    def test_context_and_measured_version_must_match(self):
        baseline=collected(BASE,[['chr1',0,10]])
        target=collected(TARGET,[])
        for r in target:
            r['local_design_id']='different-reference-or-windows'
        samples, _, _=self.render(baseline+target)
        row=next(r for r in samples if r['target_cnvkit_version']==TARGET)
        self.assertEqual(row['status'],'pair_context_mismatch')
        self.assertIsNone(row['retention_delta'])
        for r in target:
            r.update(actual_cnvkit_version=None)
        samples, _, _=self.render(baseline+target)
        row=next(r for r in samples if r['target_cnvkit_version']==TARGET)
        self.assertEqual(row['status'],'target_unavailable')

    def test_repeated_seeds_require_one_actual_run_design_for_each_component(self):
        target = collected(TARGET, [])
        target[0]['local_design_id'] = 'different_windows_across_rows'
        with self.assertRaisesRegex(ValueError, 'Inconsistent source filter or local design'):
            self.render(collected(BASE, []) + target)

    def test_missing_reference_version_historical_and_bamless(self):
        for records, mode, source, expected in ((None,'with_bam',BASE,'not_recorded_in_run'),
                                               (None,'results_only',BASE,'not_applicable_results_only'),
                                               ([],'with_bam',None,'source_cnvkit_version_unknown')):
            with self.subTest(expected=expected):
                rows, _, _ = self.render(records,[sample(source=source)],mode)
                self.assertTrue(rows)
                self.assertTrue(all(r['status']==expected and r['interval_jaccard'] is None for r in rows))

    def test_cohort_counts_independent_patients_and_no_duplicate_suite_observations(self):
        records=[]; samples=[]
        for sid,keep in [('a',True),('b',True),('c',False)]:
            samples.append(sample(sid))
            records.extend(collected(BASE,[['chr1',0,15]],sid))
            records.extend(collected(TARGET,[['chr1',0,15]] if keep else [],sid))
        _, cohort, _=self.render(records,samples,patients={'a':'p','b':'p','c':'q'})
        row=next(r for r in cohort if r['target_cnvkit_version']==TARGET)
        self.assertEqual(row['interval_jaccard'],.5)
        self.assertEqual(row['retention_delta'],-.5)
        self.assertEqual(row['interval_jaccard_n_samples'],3)
        self.assertEqual(row['retention_delta_n_independent_units'],2)
        self.assertEqual(row['missed_count'],2)
        self.assertEqual(row['paired_seed_count'],6)


if __name__ == '__main__':
    unittest.main()
