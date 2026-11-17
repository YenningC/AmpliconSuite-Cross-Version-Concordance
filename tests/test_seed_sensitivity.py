import copy
import tempfile
import unittest
from pathlib import Path

from ascompare.seed_sensitivity import analyze_seed_sensitivity, write_sensitivity_report

ALL_TOLERANCES=('detection_loss','category_loss','interval_loss','comparison_error','comparison_event_loss')


def sync_features(records):
    for r in records:
        r['features']=[{'label':label,'intervals':copy.deepcopy(r['intervals'])} for label in r['labels']]


def fixture():
    rows=[]
    for sample in ('S1','S2','S3'):
        for version in ('old','ref'):
            for seed in ('a','b','c'):
                # Independent seeds all recover the same structure, so one is sufficient
                # only for this synthetic cached baseline; safe recommendation is size 2.
                rows.append(dict(sample_id=sample,seed_id=seed,version=version,status='ok',
                    labels=['ecDNA'],intervals=[['chr1',100,200]],supported_labels=['ecDNA']))
    options=dict(versions=['old','ref'],reference_version='ref',candidate_seeds={s:['a','b','c'] for s in ('S1','S2','S3')},
        baseline_design='independent_fixed_seed_runs',sizes=[1,2,3],repeats=20,bootstrap=100,seed=7)
    sync_features(rows)
    return rows,options


class SensitivityTests(unittest.TestCase):
    def test_incomplete_baseline_blocks_selection(self):
        records,options=fixture()
        result=analyze_seed_sensitivity(records[:-1],options)
        self.assertFalse(result['ready'])
        self.assertIsNone(result['recommendation'])
        self.assertIn('Missing baseline record',result['baseline_errors'][0])

    def test_failed_baseline_is_not_a_negative(self):
        records,options=fixture();records[0]['status']='failed'
        result=analyze_seed_sensitivity(records,options)
        self.assertFalse(result['ready'])

    def test_empty_candidate_samples_need_sentinel_and_do_not_inflate_sample_size(self):
        records,options=fixture()
        options['candidate_seeds']['S4']=[]
        records += [dict(sample_id='S4',seed_id='__no_seeds__',version=v,status='no_seeds') for v in options['versions']]
        result=analyze_seed_sensitivity(records,options)
        self.assertTrue(result['ready'])
        self.assertEqual(result['n_samples'],3)
        self.assertEqual(result['selected_sample_count'],4)
        self.assertEqual(result['no_seed_samples'],['S4'])
        self.assertTrue(all(r['sample_id']!='S4' for r in result['subset_rows']))

    def test_no_tolerance_no_recommendation(self):
        records,options=fixture()
        result=analyze_seed_sensitivity(records,options)
        self.assertTrue(result['ready'])
        self.assertIsNone(result['recommendation'])
        self.assertTrue(all(r['mean_loss'] in (None,0) for r in result['summary_rows']))

    def test_next_tested_size_above_stable_minimum(self):
        records,options=fixture();options['tolerances']={m:.1 for m in ALL_TOLERANCES}
        for r in records:
            r['supported_labels']=['ecDNA','BFB']
            if r['version']=='old':r['labels']=['BFB']
        sync_features(records)
        result=analyze_seed_sensitivity(records,options)
        self.assertEqual(result['minimum_observed_size'],1)
        self.assertEqual(result['recommendation']['seed_count_per_sample_cap'],2)
        self.assertEqual(len(result['leave_one_sample_out']),3)
        self.assertEqual(result['inference'],'exploratory_only_not_population_validated')

    def test_identical_seed_subsets_across_versions_and_nested_sizes(self):
        records,options=fixture();result=analyze_seed_sensitivity(records,options)
        index={(r['size'],r['repeat'],r['sample_id'],r['version']):r for r in result['subset_rows']}
        for repeat in range(1,21):
            for sample in ('S1','S2','S3'):
                small=index[1,repeat,sample,'old']['selected_seed_ids']
                large=index[2,repeat,sample,'old']['selected_seed_ids']
                self.assertEqual(small,large[:1])
                self.assertEqual(large,index[2,repeat,sample,'ref']['selected_seed_ids'])

    def test_loss_of_distinct_categories_and_coverage(self):
        records,options=fixture()
        for r in records:
            r['supported_labels']=['ecDNA','BFB']
            if r['seed_id']=='b':
                r['labels']=['BFB'];r['intervals']=[['chr1',200,300]]
            if r['seed_id']=='c' and r['version']=='ref':r['labels']=['BFB']
        sync_features(records)
        options['tolerances']={m:0 for m in ALL_TOLERANCES}
        result=analyze_seed_sensitivity(records,options)
        self.assertEqual(result['minimum_observed_size'],3)
        losses=[r for r in result['subset_rows'] if r['size']==1 and r['selected_seed_ids']==['a']]
        self.assertEqual(losses[0]['lost_categories'],['BFB'])
        self.assertEqual(losses[0]['lost_interval_bp'],100)
        self.assertEqual(losses[0]['category_loss'],.5)
        self.assertEqual(losses[0]['interval_loss'],.5)

    def test_all_negative_no_false_adequacy(self):
        records,options=fixture()
        for r in records:r.update(status='not_detected',labels=[],intervals=[])
        sync_features(records)
        options['tolerances']={m:.5 for m in ALL_TOLERANCES}
        result=analyze_seed_sensitivity(records,options)
        self.assertIsNone(result['recommendation'])
        self.assertTrue(all(r['n_informative_samples']==0 for r in result['summary_rows']))

    def test_single_prespecified_order_not_random_stability(self):
        records,options=fixture()
        options.update(mode='prespecified_nested',seed_orders=options['candidate_seeds'],tolerances={m:.1 for m in ALL_TOLERANCES})
        result=analyze_seed_sensitivity(records,options)
        self.assertIsNone(result['recommendation'])
        self.assertEqual(result['repeats'],1)

    def test_report_saves_scope_and_warning(self):
        records,options=fixture()
        result=analyze_seed_sensitivity(records,options)
        with tempfile.TemporaryDirectory() as folder:
            write_sensitivity_report(folder,result)
            self.assertTrue((Path(folder)/'subset_losses.tsv').exists())
            text=(Path(folder)/'README.md').read_text()
            self.assertIn('固定参考版本：`ref`',text)
            self.assertIn('不是充分种子数',text)

    def test_union_preserved_but_version_differences_lost(self):
        records,options=fixture()
        for r in records:
            if (r['seed_id']=='a' and r['version']=='old') or (r['seed_id']=='b' and r['version']=='ref'):
                r.update(status='not_detected',labels=[],intervals=[])
        sync_features(records)
        options['tolerances']={m:0 for m in ALL_TOLERANCES}
        result=analyze_seed_sensitivity(records,options)
        rows=[r for r in result['comparison_subset_rows'] if r['size']==1 and r['selected_seed_ids']==['c']]
        self.assertTrue(rows)
        self.assertEqual(rows[0]['lost_difference_events'],['extra_detection','missed_detection'])
        self.assertEqual(rows[0]['comparison_event_loss'],1)
        self.assertAlmostEqual(rows[0]['comparison_error'],2/3)
        union=[r for r in result['subset_rows'] if r['size']==1 and r['selected_seed_ids']==['c']]
        self.assertTrue(all(r['detection_loss']==r['category_loss']==r['interval_loss']==0 for r in union))
        self.assertEqual(result['minimum_observed_size'],3)

    def test_no_difference_events_cannot_validate_retention(self):
        records,options=fixture();options['tolerances']={m:.1 for m in ALL_TOLERANCES}
        result=analyze_seed_sensitivity(records,options)
        self.assertIsNone(result['recommendation'])
        self.assertIn('no_baseline_version_difference_events',result['size_rows'][0]['reason'])

    def test_comparison_tolerances_must_be_explicit(self):
        records,options=fixture()
        for r in records:
            r['supported_labels']=['ecDNA','BFB']
            if r['version']=='old':r['labels']=['BFB']
        sync_features(records)
        options['tolerances']={m:.1 for m in ALL_TOLERANCES if m!='comparison_event_loss'}
        result=analyze_seed_sensitivity(records,options)
        self.assertIsNone(result['recommendation'])

    def test_spatial_class_swap_cannot_disappear_from_selected_subset(self):
        records,options=fixture()
        for r in records:
            r.update(labels=['ecDNA','BFB'],supported_labels=['ecDNA','BFB'],intervals=[['chr1',0,100],['chr2',0,100]])
            swap=r['version']=='old' and r['seed_id']=='a'
            r['features']=[{'label':'ecDNA','intervals':[['chr2' if swap else 'chr1',0,100]]},
                           {'label':'BFB','intervals':[['chr1' if swap else 'chr2',0,100]]}]
        options['tolerances']={m:0 for m in ALL_TOLERANCES}
        result=analyze_seed_sensitivity(records,options)
        selected=[r for r in result['comparison_subset_rows'] if r['size']==1 and r['selected_seed_ids']==['c']]
        self.assertTrue(selected)
        self.assertAlmostEqual(selected[0]['comparison_error'],1/3)
        self.assertEqual(selected[0]['lost_difference_events'],['spatial_class_change'])
        self.assertEqual(selected[0]['comparison_event_loss'],1)
        self.assertEqual(selected[0]['absolute_metric_deviations']['interval_distance'],0)
        self.assertEqual(selected[0]['absolute_metric_deviations']['label_discordance'],0)
        self.assertEqual(result['minimum_observed_size'],3)

    def test_one_missing_feature_baseline_blocks_recommendation_for_entire_scope(self):
        records,options=fixture()
        for r in records:
            r['supported_labels']=['ecDNA','BFB']
            if r['version']=='old':r['labels']=['BFB']
        sync_features(records);del records[0]['features']
        options['tolerances']={m:.9 for m in ALL_TOLERANCES}
        result=analyze_seed_sensitivity(records,options)
        self.assertTrue(result['ready'])
        self.assertFalse(result['comparison_scope_complete'])
        self.assertIsNone(result['recommendation'])
        self.assertEqual(result['recommendation_reason'],'class_footprint_information_missing_or_unavailable')


if __name__=='__main__':unittest.main()
