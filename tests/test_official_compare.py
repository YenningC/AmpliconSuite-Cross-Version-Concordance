import csv
import json
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from ascompare.official_compare import export_cycles, _validate_results_table, _output, compare_completed_records, valid_ec_metrics, summarize_official_metrics
from ascompare.full_suite import WHOLE_SAMPLE


class OfficialAdapterTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.cycles=self.root/'input_cycles.txt'
        self.cycles.write_text('Segment\t1\tchr1\t100\t199\nSegment\t2\tchr2\t1000\t1000\n'
                               'Cycle=1;Copy_count=8;Segments=1+,2-,1+\n'
                               'Cycle=2;Copy_count=3.5;Segments=0+,1-,2+,0-\n')

    def tearDown(self):
        self.temp.cleanup()

    def test_order_repeats_coordinates_and_path_closure(self):
        out=self.root/'transport.bed'
        metadata=export_cycles([self.cycles],out)
        with out.open() as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(metadata['rows'],5)
        self.assertEqual(rows[0]['end'],'200')
        self.assertEqual(rows[1]['end'],'1001')
        self.assertEqual([r['strand'] for r in rows],['+','-','+','-','+'])
        self.assertEqual([r['iscyclic'] for r in rows],['True','True','True','False','False'])
        self.assertEqual(rows[0]['circ_id'],rows[2]['circ_id'])
        self.assertEqual(rows[3]['estimated_cn'],'3.5')

    def test_amplicon_cycle_ids_remain_distinct(self):
        out=self.root/'transport.bed';metadata=export_cycles([self.cycles,self.cycles],out)
        with out.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(len({r['circ_id'] for r in rows}),4)
        self.assertFalse(metadata['inter_amplicon_edges_added'])
        self.assertEqual(metadata['amplicon_file_count'],2)
        self.assertTrue(all(len({r['circ_id'].split('_')[0] for r in rows if r['circ_id']==cycle})==1
                            for cycle in {r['circ_id'] for r in rows}))

    def test_internal_source_sink_refused(self):
        self.cycles.write_text(self.cycles.read_text().replace('1+,2-,1+','1+,0+,2-'))
        with self.assertRaises(ValueError):export_cycles([self.cycles],self.root/'out.bed')

    def test_missing_segment_refused(self):
        self.cycles.write_text(self.cycles.read_text().replace('1+,2-,1+','1+,7-,1+'))
        with self.assertRaises(ValueError):export_cycles([self.cycles],self.root/'out.bed')

    def test_missing_coordinates_not_silently_ignored(self):
        table=self.root/'result.tsv'
        table.write_text('Sample name\tFeature ID\tClassification\tLocation\nS\tS_ecDNA_1\tecDNA\t[]\n')
        with self.assertRaises(ValueError):_validate_results_table(table)

    def test_partial_malformed_native_coordinates_are_rejected(self):
        table=self.root/'result.tsv'
        for location in ("['chr1:100-200', 'invalid']","['chr1:200-100']","['chr1:100-100']"):
            with self.subTest(location=location):
                table.write_text('Sample name\tFeature ID\tClassification\tLocation\nS\tS_ecDNA_1\tecDNA\t'+location+'\n')
                with self.assertRaises(ValueError):_validate_results_table(table)
        table.write_text("Sample name\tFeature ID\tClassification\tLocation\nS\tS_NA\tNone\t[]\n")
        self.assertEqual(len(_validate_results_table(table)),1)

    def test_output_cannot_modify_source_tree(self):
        with self.assertRaises(ValueError):_output(self.root/'new',[self.root])

    def test_optional_diagnostics_keep_missing_and_failed_distinct(self):
        records=[dict(version=v,sample_id='S',seed_id='seed',status='ok',run_dir=str(self.root/v)) for v in ('ref','old')]
        options={'results':{'source':'unused','python':'unused'}}
        rows=compare_completed_records(records,self.root/'report',options,'ref')
        self.assertEqual(rows[0]['status'],'missing')
        records[1]['status']='failed'
        rows=compare_completed_records(records,self.root/'report',options,'ref')
        self.assertEqual(rows[0]['status'],'not_evaluable')

    def test_ec_multiamplicon_defaults_for_whole_sample_only(self):
        records=[]
        for version in ('ref','old'):
            run=self.root/version;(run/'aa').mkdir(parents=True)
            for n in (1,2):
                (run/'aa'/f'reconstruction_amplicon{n}_cycles.txt').write_text(self.cycles.read_text())
            records.append(dict(version=version,sample_id='S',seed_id='seed',status='ok',run_dir=str(run)))
        options={'ec':{'source':'unused','python':'unused'}}
        rows=compare_completed_records(records,self.root/'report',options,'ref')
        self.assertEqual(rows[0]['status'],'unsupported')
        for record in records:record['seed_id']=WHOLE_SAMPLE
        with patch('ascompare.official_compare.ec_comparator',return_value={'status':'complete'}) as caller:
            rows=compare_completed_records(records,self.root/'report',options,'ref')
        self.assertEqual(rows[0]['status'],'complete')
        self.assertEqual(rows[0]['seed_id'],'__whole_sample__')
        self.assertIn('/whole_sample/',str(caller.call_args.kwargs['out']))
        self.assertEqual(len(caller.call_args.args[0]),2)
        self.assertEqual(len(caller.call_args.args[1]),2)

    def test_optional_plot_failure_retained_without_changing_run_status(self):
        records=[dict(version=v,sample_id='S',seed_id=WHOLE_SAMPLE,status='failed') for v in ('ref','old')]
        diagnostic=self.root/'report/official_comparisons';diagnostic.mkdir(parents=True)
        (diagnostic/'official_plots.tsv').write_text('path\nstale_high_similarity.png\n')
        with patch('ascompare.official_plots.plot_official_metrics',side_effect=ImportError('test missing matplotlib')):
            rows=compare_completed_records(records,self.root/'report',{'ec':{'source':'unused'}},'ref')
        self.assertEqual(rows[0]['status'],'not_evaluable')
        status=json.loads((self.root/'report/official_comparisons/plots_status.json').read_text())
        self.assertEqual(status['status'],'failed')
        self.assertIn('ImportError',status['reason'])
        self.assertNotIn('stale_high_similarity',(diagnostic/'official_plots.tsv').read_text())

    def test_ec_metrics_accept_description_but_reject_missing_nan(self):
        keys=('cn_hamming_dist','cn_hamming_norm_dist','cn_cos_dist','cn_jc_dist','breakpoint_dist','fragments_dist','cycles_dist','total_cost')
        data={'distances':dict.fromkeys(keys,0.)};data['distances']['total_cost_description']={'cycles_dist':.5}
        self.assertTrue(valid_ec_metrics(data))
        data['distances']['cn_cos_dist']=float('nan');self.assertFalse(valid_ec_metrics(data))
        self.assertFalse(valid_ec_metrics({}))
        self.assertFalse(valid_ec_metrics({'distances':None}))

    def test_native_metric_units_status_and_no_probability(self):
        run=self.root/'diagnostic';(run/'official_raw').mkdir(parents=True)
        (run/'provenance.json').write_text('{"tool_version":"0.1.3"}')
        (run/'official_raw/metrics.json').write_text(json.dumps({'distances':{'cn_hamming_dist':2000,'cn_jc_dist':.3,'total_cost':1.5,'total_cost_description':{'x':1}}}))
        rows=summarize_official_metrics([dict(reference_version='ref',candidate_version='old',sample_id='S',seed_id='R',tool='ec',status='complete',directory=str(run)),
                                         dict(tool='ec',status='unsupported',reason='missing structure')],self.root)
        self.assertEqual(len(rows),4)
        self.assertEqual(rows[0]['unit'],'bp')
        self.assertEqual(rows[1]['direction'],'lower_is_closer')
        self.assertEqual(rows[2]['documented_max'],5)
        self.assertIsNone(rows[3]['value'])
        self.assertTrue(all(r['compatibility_cutoff'] is None for r in rows))


if __name__=='__main__':unittest.main()
