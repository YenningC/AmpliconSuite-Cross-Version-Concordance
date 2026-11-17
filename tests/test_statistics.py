import copy
import math
import tempfile
import unittest
from pathlib import Path

from ascompare.statistics import analyze, interval_jaccard, merge_intervals, subtract_intervals, build_sample_similarity, AGREEMENT_METRICS
from ascompare.reporting import write_report


def record(sample, seed, version, detected=True, labels=None, intervals=None, **extra):
    extra.setdefault('supported_labels',['ecDNA','BFB','Linear','Complex-non-cyclic'])
    result = dict(sample_id=sample,seed_id=seed,version=version,status="ok" if detected else "not_detected",
                labels=labels if labels is not None else ["ecDNA"] if detected else [],
                intervals=intervals if intervals is not None else [["chr1",10,20]] if detected else [],**extra)
    result.setdefault("features", [{"label":label,"intervals":result["intervals"]} for label in result["labels"]])
    return result


OPTIONS = dict(bootstrap=100,permutations=100,min_samples=2,seed=17)


class IntervalTests(unittest.TestCase):
    def test_merges_and_bp_jaccard(self):
        a = [["chr1",0,10],["chr1",5,15],["chr2",0,5]]
        b = [["chr1",10,20]]
        self.assertEqual(merge_intervals(a),[["chr1",0,15],["chr2",0,5]])
        self.assertEqual(subtract_intervals(a,b),[["chr1",0,10],["chr2",0,5]])
        self.assertAlmostEqual(interval_jaccard(a,b),5/25)
        self.assertIsNone(interval_jaccard([],[]))
        self.assertEqual(interval_jaccard(a,[]),0)

    def test_boundaries_and_invalid(self):
        self.assertEqual(interval_jaccard([["c",0,10]],[["c",10,20]]),0)
        with self.assertRaises(ValueError):
            merge_intervals([["c",5,4]])


class PairedTests(unittest.TestCase):
    def test_sample_weight_not_seed_weight(self):
        records = []
        for seed in range(8):
            records += [record("A",str(seed),"old",False),record("A",str(seed),"ref")]
        records += [record("B","x","old"),record("B","x","ref")]
        _,stats,_ = analyze(records,"ref",OPTIONS)
        row = next(r for r in stats if r["metric"]=="detection_difference" and r["subset"]=="all_evaluable")
        self.assertEqual(row["n_samples"],2)
        self.assertEqual(row["n_regions"],9)
        self.assertEqual(row["estimate"],-.5)

    def test_patient_weight_and_bounds_use_independent_units(self):
        records = [record(s, 'seed', version, detected=(version == 'old') == (s != 'third'))
                   for s in ('first', 'second', 'third') for version in ('old', 'ref')]
        patient_ids = {'first': 'P1', 'second': 'P1', 'third': 'P2'}
        _, stats, compatibility = analyze(records, 'ref', dict(OPTIONS, patient_ids=patient_ids, min_samples=3))
        delta = next(r for r in stats if r['metric'] == 'detection_difference' and r['subset'] == 'all_evaluable')
        self.assertEqual((delta['n_independent_units'], delta['n_samples'], delta['n_regions']), (2, 3, 3))
        self.assertEqual(delta['independence_unit'], 'patient')
        self.assertEqual(delta['estimate'], 0.)
        self.assertFalse(delta['sufficient_samples'])
        self.assertIn('too_few_independent_samples', compatibility[0]['reasons'])
        with self.assertRaisesRegex(ValueError, 'every analyzed sample'):
            analyze(records, 'ref', dict(OPTIONS, patient_ids={'first': 'P1'}))

    def test_primary_pair_selection_controls_only_inference(self):
        records = [record(str(i), 'seed', v) for v in ('old', 'new', 'ref') for i in range(30)]
        options = dict(OPTIONS, margins={m: .8 for m in AGREEMENT_METRICS})
        for specification, n in (([['old', 'new']], 1), ('reference', 2), ([], 0)):
            with self.subTest(specification=specification):
                _, stats, compatible = analyze(records, 'ref', dict(options, primary_pairs=specification))
                self.assertEqual(sum(r['primary_pair'] for r in compatible), n)
                self.assertEqual(sum(r['decision'] == 'within_margin' for r in compatible), n)
                for row in compatible:
                    if row['primary_pair']:
                        self.assertEqual(row['per_metric_alpha'], .05 / (n * len(AGREEMENT_METRICS)))
                    else:
                        self.assertEqual(row['reasons'], ['exploratory_pair_not_primary'])
                        self.assertIsNone(row['per_metric_alpha'])
                exploratory = [r for r in stats if not r['primary_pair']]
                self.assertTrue(all(r['estimate'] is not None and r['p_difference'] is None for r in exploratory))

    def test_multi_labels_and_change_coordinates(self):
        records = [record("A","x","old",labels=["ecDNA","BFB"],intervals=[["c",0,10]]),
                   record("A","x","ref",labels=["ecDNA","Linear"],intervals=[["c",5,15]])]
        rows,_,_ = analyze(records,"ref",OPTIONS)
        self.assertEqual(rows[0]["labels_added"],["BFB"])
        self.assertEqual(rows[0]["labels_removed"],["Linear"])
        self.assertEqual(rows[0]["intervals_added"],[["c",0,5]])
        self.assertEqual(rows[0]["intervals_removed"],[["c",10,15]])
        self.assertEqual(rows[0]["class_change"],1)

    def test_missing_failed_unsupported_are_not_negative(self):
        records=[record("A","x","ref"),record("B","x","ref"),record("B","x","old")]
        records[-1]["status"]="failed"
        rows,stats,compat=analyze(records,"ref",OPTIONS)
        self.assertEqual([r["status"] for r in rows],["missing","failed"])
        self.assertTrue(all(r["detection_difference"] is None for r in rows))
        self.assertEqual(compat[0]["evaluable_regions"],0)
        self.assertEqual(compat[0]["decision"],"insufficient_evidence")

    def test_all_negative_does_not_establish_compatibility(self):
        records=[record(str(i),"x",v,False) for i in range(30) for v in ("old","ref")]
        _,stats,compat=analyze(records,"ref",dict(OPTIONS,margins={m:.9 for m in AGREEMENT_METRICS}))
        self.assertEqual(compat[0]["decision"],"insufficient_evidence")
        self.assertEqual(compat[0]["signal_samples"],0)
        row=next(r for r in stats if r["subset"]=="all_evaluable" and r["metric"]=="detection_agreement")
        self.assertEqual(row["estimate"],1)

    def test_not_started_jobs_do_not_count_as_negative_agreement(self):
        records=[record('A','x','ref'),record('A','x','old')]
        records[1]['status']='not_run'
        rows,_,compat=analyze(records,'ref',OPTIONS)
        self.assertFalse(rows[0]['evaluable'])
        self.assertIsNone(rows[0]['detection_agreement'])
        self.assertEqual(compat[0]['decision'],'insufficient_evidence')

    def test_zero_difference_interval_not_degenerate_for_decision(self):
        records=[record(str(i),"x",v) for i in range(3) for v in ("old","ref")]
        _,stats,compat=analyze(records,"ref",dict(OPTIONS,margins={m:.1 for m in AGREEMENT_METRICS}))
        row=next(r for r in stats if r["subset"]=="signal_positive" and r["metric"]=="detection_discordance")
        self.assertEqual(row["bootstrap_ci_high"],0)
        self.assertGreater(row["agreement_upper"],.1)
        self.assertEqual(compat[0]["decision"],"insufficient_evidence")

    def test_enough_data_and_explicit_margin_can_pass(self):
        records=[record(str(i),"x",v) for i in range(100) for v in ("old","ref")]
        margins={m:.2 for m in AGREEMENT_METRICS}
        _,_,compat=analyze(records,"ref",dict(OPTIONS,margins=margins))
        self.assertEqual(compat[0]["decision"],"within_margin")
        _,_,compat=analyze(records,"ref",OPTIONS)
        self.assertEqual(compat[0]["decision"],"insufficient_evidence")

    def test_support_difference_prevents_full_compatibility(self):
        records=[record(str(i),"x",v,supported_labels=["ecDNA"] if v=="old" else ["ecDNA","FAN"]) for i in range(30) for v in ("old","ref")]
        rows,_,compat=analyze(records,"ref",dict(OPTIONS,margins={m:.9 for m in AGREEMENT_METRICS}))
        self.assertFalse(rows[0]["full_label_comparable"])
        self.assertEqual(rows[0]["label_scope"],["ecDNA"])
        self.assertEqual(compat[0]["decision"],"insufficient_evidence")

    def test_unknown_category_support_cannot_pass(self):
        records=[record(str(i),"x",v) for i in range(30) for v in ("old","ref")]
        for r in records:
            del r['supported_labels']
        rows,_,compat=analyze(records,"ref",dict(OPTIONS,margins={m:.9 for m in AGREEMENT_METRICS}))
        self.assertFalse(rows[0]['full_label_comparable'])
        self.assertEqual(compat[0]['decision'],'insufficient_evidence')

    def test_all_pairs_and_missing_version(self):
        records=[record("A","x",v) for v in ("old","ref")]
        rows,_,compat=analyze(records,"ref",dict(OPTIONS,versions=["old","mid","ref"]))
        self.assertEqual(len(compat),3)
        self.assertEqual(len(rows),3)
        self.assertEqual(sum(r["pair_kind"]=="historical_pair" for r in rows),1)

    def test_duplicate_rejected(self):
        item=record("A","x","ref")
        with self.assertRaises(ValueError):
            analyze([item,item],"ref",OPTIONS)

    def test_class_locations_swap_despite_identical_labels_and_overall_intervals(self):
        records=[]
        for i in range(100):
            for v in ('old','ref'):
                labels=['ecDNA','BFB'];ints=[['chr1',0,100],['chr2',0,100]]
                features=[{'label':'ecDNA','intervals':[ints[0 if v=='ref' else 1]]},
                          {'label':'BFB','intervals':[ints[1 if v=='ref' else 0]]}]
                records.append(record(str(i),'__whole_sample__',v,labels=labels,intervals=ints,features=features))
        pairs,stats,compat=analyze(records,'ref',dict(OPTIONS,analysis_scope='full_suite',margins={m:.2 for m in AGREEMENT_METRICS}))
        self.assertEqual(pairs[0]['label_discordance'],0)
        self.assertEqual(pairs[0]['interval_jaccard'],1)
        self.assertEqual(pairs[0]['class_footprint_distance'],1)
        self.assertEqual(len(pairs[0]['class_footprint_details']),2)
        self.assertTrue(all(r['distance']==1 for r in pairs[0]['class_footprint_details']))
        self.assertEqual(compat[0]['decision'],'exceeds_margin')
        self.assertIn('same_cycle_order',compat[0]['compatibility_metric_scope']['not_established'])

    def test_missing_features_never_become_perfect_spatial_agreement(self):
        records=[record(str(i),'x',v) for i in range(100) for v in ('old','ref')]
        for r in records:del r['features']
        pairs,_,compat=analyze(records,'ref',dict(OPTIONS,margins={m:.2 for m in AGREEMENT_METRICS}))
        self.assertTrue(all(r['class_footprint_distance'] is None for r in pairs))
        self.assertEqual(compat[0]['decision'],'insufficient_evidence')
        self.assertIn('class_footprint_information_missing',compat[0]['reasons'])

    def test_missing_same_class_feature_detected_by_union_validation(self):
        records=[record('A','x',v,intervals=[['c',0,10],['c',20,30]]) for v in ('old','ref')]
        records[0]['features'][0]['intervals']=[['c',0,10]]
        pairs,_,_=analyze(records,'ref',OPTIONS)
        self.assertIsNone(pairs[0]['class_footprint_distance'])
        self.assertIn('feature_union_differs',pairs[0]['class_footprint_reason'])

    def test_spatial_metric_requires_new_margin_and_uses_macro_label_mean(self):
        records=[record(str(i),'x',v) for i in range(100) for v in ('old','ref')]
        old_margins={m:.2 for m in AGREEMENT_METRICS if m!='class_footprint_distance'}
        _,_,compat=analyze(records,'ref',dict(OPTIONS,margins=old_margins))
        self.assertEqual(compat[0]['decision'],'insufficient_evidence')
        a=record('A','x','old',labels=['ecDNA','BFB'],intervals=[['c',0,10],['c',100,110]],
                 features=[{'label':'ecDNA','intervals':[['c',0,10]]},{'label':'BFB','intervals':[['c',100,110]]}])
        b=record('A','x','ref',labels=['ecDNA'],intervals=[['c',0,10]])
        pairs,_,_=analyze([a,b],'ref',OPTIONS)
        self.assertEqual(pairs[0]['class_footprint_jaccard'],.5)

    def test_no_seeds_preserved_but_not_counted_as_regions_or_negatives(self):
        records=[record('positive','a',v) for v in ('old','ref')]
        records += [dict(sample_id=f'no_seed_{i}',seed_id='__no_seeds__',version=v,status='no_seeds',labels=[],intervals=[]) for i in range(3) for v in ('old','ref')]
        rows,stats,compat=analyze(records,'ref',OPTIONS)
        self.assertEqual(len(rows),4)
        self.assertEqual(compat[0]['total_regions'],1)
        self.assertEqual(compat[0]['excluded_regions'],0)
        self.assertEqual(compat[0]['no_seed_samples'],3)
        detection=next(r for r in stats if r['metric']=='detection_agreement' and r['subset']=='all_evaluable')
        self.assertEqual(detection['n_samples'],1)
        self.assertEqual(detection['n_regions'],1)
        self.assertEqual(compat[0]['decision'],'insufficient_evidence')

    def test_full_suite_zero_seed_success_is_pipeline_negative_not_aa_agreement(self):
        records=[]
        for i in range(5):
            for v in ('old','ref'):
                positive=i<2
                records.append(record(str(i),'__whole_sample__',v,positive,
                    seed_count=1 if positive else 0,seed_intervals=[['chr1',10,20]] if positive else [],
                    detection_stage='ac' if positive else 'seed_processing',aa_executed=positive,ac_executed=positive))
        rows,stats,compat=analyze(records,'ref',dict(OPTIONS,analysis_scope='full_suite'))
        allrow=next(r for r in stats if r['metric']=='detection_agreement' and r['subset']=='all_evaluable')
        signal=next(r for r in stats if r['metric']=='detection_agreement' and r['subset']=='signal_positive')
        self.assertEqual(allrow['n_samples'],5)
        self.assertEqual(signal['n_samples'],2)
        self.assertEqual(compat[0]['no_seed_samples'],0)
        negative=next(r for r in rows if r['sample_id']=='4')
        self.assertEqual(negative['detection_stage'],'seed_processing')
        self.assertFalse(negative['aa_executed'])
        self.assertIsNone(negative['seed_interval_jaccard'])
        self.assertIn('including CNVkit',compat[0]['scope'])

    def test_full_suite_seed_stage_metrics_are_separate(self):
        records=[record('A','__whole_sample__','old',seed_count=1,seed_intervals=[['c',0,10]]),
                 record('A','__whole_sample__','ref',seed_count=2,seed_intervals=[['c',5,15],['c',15,20]])]
        rows,stats,_=analyze(records,'ref',dict(OPTIONS,analysis_scope='full_suite'))
        self.assertEqual(rows[0]['seed_count_difference'],-1)
        self.assertAlmostEqual(rows[0]['seed_interval_jaccard'],.25)
        seed_count=next(r for r in stats if r['metric']=='seed_count_difference')
        self.assertEqual(seed_count['subset'],'seed_stage_evaluable')
        self.assertIsNone(seed_count['conservative_ci_low'])

    def test_full_suite_failed_record_stays_excluded(self):
        records=[record('A','__whole_sample__','old',False),record('A','__whole_sample__','ref',False)]
        records[0]['status']='failed'
        _,stats,compat=analyze(records,'ref',dict(OPTIONS,analysis_scope='full_suite'))
        self.assertEqual(compat[0]['excluded_regions'],1)
        self.assertTrue(all(r['estimate'] is None for r in stats))

    def test_full_suite_verified_seed_outputs_can_describe_later_aa_failure(self):
        records=[record('A','__whole_sample__',v,seed_count=1,seed_intervals=[['c',0,10]]) for v in ('old','ref')]
        records[0]['status']='failed';records[0]['stage']='aa'
        rows,stats,compat=analyze(records,'ref',dict(OPTIONS,analysis_scope='full_suite'))
        self.assertFalse(rows[0]['evaluable'])
        self.assertEqual(rows[0]['seed_interval_jaccard'],1)
        self.assertEqual(next(r for r in stats if r['metric']=='seed_count_difference')['n_samples'],1)
        self.assertEqual(next(r for r in stats if r['metric']=='detection_agreement' and r['subset']=='all_evaluable')['n_samples'],0)

    def test_full_suite_report_scope_and_seed_table(self):
        records=[record('A','__whole_sample__',v,False,seed_count=0,seed_intervals=[],detection_stage='seed_processing',aa_executed=False,ac_executed=False) for v in ('old','ref')]
        pairs,stats,compat=analyze(records,'ref',dict(OPTIONS,analysis_scope='full_suite'))
        with tempfile.TemporaryDirectory() as folder:
            write_report(folder,records,pairs,stats,compat,'ref',synthetic=True)
            text=(Path(folder)/'report.md').read_text()
            self.assertIn('完整 Suite 流程',text)
            self.assertNotIn('未验证 CNVkit 版本影响',text)
            self.assertIn('不按相同 amplicon 编号对齐',text)
            self.assertTrue((Path(folder)/'seed_stage_differences.tsv').is_file())
            self.assertIn('Full Suite, including CNVkit',(Path(folder)/'version_ladder.svg').read_text())
            self.assertTrue((Path(folder)/'sample_similarity.tsv').is_file())
            self.assertTrue((Path(folder)/'sample_version_heatmap.png').is_file())
            self.assertEqual(len(list((Path(folder)/'sample_ladders').glob('*.pdf'))),1)

    def test_per_sample_table_preserves_missing_and_no_seeds(self):
        records=[record('A','x','old'),record('A','x','ref'),record('A','y','ref')]
        records += [dict(sample_id='N',seed_id='__no_seeds__',version=v,status='no_seeds',labels=[],intervals=[]) for v in ('old','ref')]
        pairs,_,_=analyze(records,'ref',OPTIONS)
        rows=build_sample_similarity(records,'ref',pairs)
        partial=next(r for r in rows if r['sample_id']=='A' and r['version']=='old')
        self.assertEqual(partial['comparison_state'],'partially_evaluable')
        self.assertEqual(partial['n_input_units'],2)
        self.assertEqual(partial['n_evaluable_units'],1)
        empty=next(r for r in rows if r['sample_id']=='N' and r['version']=='old')
        self.assertEqual(empty['comparison_state'],'no_seeds')
        self.assertIsNone(empty['detection_agreement'])

    def test_per_sample_both_negative_has_no_structural_agreement(self):
        records=[record('A','__whole_sample__',v,False,seed_count=0,seed_intervals=[]) for v in ('old','ref')]
        pairs,_,_=analyze(records,'ref',dict(OPTIONS,analysis_scope='full_suite'))
        rows=build_sample_similarity(records,'ref',pairs,'full_suite')
        self.assertTrue(all(r['comparison_state']=='both_not_detected' for r in rows))
        self.assertTrue(all(r['interval_jaccard'] is None for r in rows))
        self.assertTrue(all(r['detection_agreement']==1 for r in rows))
        self.assertTrue(all(r['seed_count']==0 for r in rows))
        self.assertTrue(all(r['uncertainty']=='no_single_sample_confidence_interval' for r in rows))

    def test_sample_plot_filename_cannot_escape_output(self):
        records=[record('../escape','x',v) for v in ('old','ref')]
        pairs,stats,compat=analyze(records,'ref',OPTIONS)
        with tempfile.TemporaryDirectory() as folder:
            write_report(folder,records,pairs,stats,compat,'ref',synthetic=True)
            files=list((Path(folder)/'sample_ladders').glob('*.svg'))
            self.assertEqual(len(files),1)
            self.assertTrue(files[0].name.startswith('escape_'))
            self.assertIn('Fixed reference: ref',files[0].read_text())

    def test_report_artifacts_and_reference(self):
        records=[record("A","x",v) for v in ("old","ref")]
        pairs,stats,compat=analyze(records,"ref",OPTIONS)
        with tempfile.TemporaryDirectory() as folder:
            write_report(folder,records,pairs,stats,compat,"ref",synthetic=True)
            out=Path(folder)
            self.assertIn("固定参考版本：`ref`",(out/"report.md").read_text())
            self.assertIn("人工合成测试数据",(out/"report.md").read_text())
            for name in ("differences.tsv","statistics.tsv","compatibility.tsv","version_ladder.png","compatibility.pdf"):
                self.assertTrue((out/name).is_file(),name)


if __name__ == "__main__":
    unittest.main()
