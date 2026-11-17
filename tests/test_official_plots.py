import importlib.util
import tempfile
import unittest
from pathlib import Path

from ascompare.common import table
from ascompare.full_suite import WHOLE_SAMPLE
from ascompare.official_plots import plot_groups, plot_official_metrics, NORMALIZED


def records():
    base=dict(sample_id='SYNTHETIC_PLOT',seed_id=WHOLE_SAMPLE,reference_version='suite-REF',candidate_version='suite-OLD',tool='ec',status='complete')
    return [{**base,'metric':metric,'value':.25} for metric,_ in NORMALIZED]+[
        {**base,'metric':'total_cost','value':1.25},
        {**base,'candidate_version':'suite-UNAVAILABLE','metric':'','value':'','status':'unsupported'}]


class OfficialPlotTests(unittest.TestCase):
    def test_native_distances_never_inverted(self):
        group=plot_groups(records())[0]
        self.assertEqual(group['points']['suite-OLD','cn_cos_dist']['value'],.25)
        self.assertEqual(group['points']['suite-OLD','total_cost']['value'],1.25)
        self.assertEqual(group['statuses']['suite-UNAVAILABLE'],'unsupported')

    def test_selected_seeds_not_pooled(self):
        rows=records();rows.append({**rows[0],'seed_id':'another_seed'})
        self.assertEqual(len(plot_groups(rows)),2)

    def test_invalid_values_not_plotted_as_zero(self):
        rows=records();rows[0]['value']='nan'
        group=plot_groups(rows)[0]
        self.assertIsNone(group['points']['suite-OLD','cn_hamming_norm_dist']['value'])

    def test_duplicate_metrics_rejected(self):
        rows=records()
        with self.assertRaises(ValueError):plot_groups(rows+[rows[0]])

    @unittest.skipUnless(importlib.util.find_spec('matplotlib'),'matplotlib optional in minimal test interpreter')
    def test_export_six_panels_and_separate_cost_all_formats(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'official_comparisons').mkdir()
            table(root/'official_comparisons/official_metrics.tsv',records())
            outputs=plot_official_metrics(root)
            self.assertEqual(len(outputs),6)
            self.assertEqual({Path(r['path']).suffix for r in outputs},{'.png','.pdf','.svg'})
            self.assertTrue(all(Path(r['path']).stat().st_size>1000 for r in outputs))
            svg=next(Path(r['path']) for r in outputs if r['figure']=='normalized_distances' and r['path'].endswith('.svg')).read_text()
            self.assertIn('suite-REF',svg)
            self.assertIn('SYNTHETIC TEST DATA',svg)
            self.assertIn('Whole-sample reconstructed path/cycle collection',svg)
            self.assertNotIn('Selected seed: __whole_sample__',svg)
            self.assertIn('unsupported',svg)
            self.assertNotIn('Confidence interval',svg)

    @unittest.skipUnless(importlib.util.find_spec('matplotlib'),'matplotlib optional in minimal test interpreter')
    def test_whole_sample_reference_failure_uses_matching_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'official_comparisons').mkdir()
            table(root/'official_comparisons/official_metrics.tsv',records())
            outputs=plot_official_metrics(root,records=[dict(sample_id='SYNTHETIC_PLOT',seed_id=WHOLE_SAMPLE,
                version='suite-REF',status='failed',reference_scope='human_only_incomplete_reference')])
            svg=next(Path(r['path']) for r in outputs if r['figure']=='normalized_distances' and r['path'].endswith('.svg')).read_text()
            self.assertIn('reference failed',svg)
            self.assertIn('viral contigs unavailable',svg)

    @unittest.skipUnless(importlib.util.find_spec('matplotlib'),'matplotlib optional in minimal test interpreter')
    def test_negative_without_reconstruction_has_no_zero_reference_anchor(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'official_comparisons').mkdir()
            row=records()[0];row.update(status='not_evaluable',metric='',value='')
            table(root/'official_comparisons/official_metrics.tsv',[row])
            outputs=plot_official_metrics(root,records=[dict(sample_id='SYNTHETIC_PLOT',seed_id=WHOLE_SAMPLE,
                version='suite-REF',status='not_detected')])
            svg=next(Path(r['path']) for r in outputs if r['figure']=='normalized_distances' and r['path'].endswith('.svg')).read_text()
            self.assertIn('no reconstructed structure',svg)
            self.assertIn('not evaluable',svg)


if __name__=='__main__':unittest.main()
