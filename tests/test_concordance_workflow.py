"""Public loader/orchestrator integration; synthetic result files, no BAM calls."""
from contextlib import ExitStack
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ascompare.concordance import load_manifest,main,run_concordance
from ascompare.concordance_demo import write_demo_manifest,SOURCE_VERSION,REFERENCE_VERSION
from ascompare.common import atomic_json,digest


class ConcordanceWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name).resolve()
        self.manifest=write_demo_manifest(self.root/'synthetic_inputs')

    def tearDown(self):self.temp.cleanup()

    def change_manifest(self,update):
        data=json.loads(self.manifest.read_text());update(data);atomic_json(self.manifest,data)

    def run_results(self,out=None,plots=False):
        config,lock,entries=load_manifest(self.manifest)
        out=Path(out or self.root/'output')
        with ExitStack() as stack:
            for target in ('ascompare.concordance._validate_bam','ascompare.concordance.run_bounded_downstream',
                           'ascompare.concordance.verify_runtime','ascompare.concordance_cnv.run_local_upstream',
                           'ascompare.concordance_import.reclassify_result'):
                stack.enter_context(patch(target,side_effect=AssertionError('No BAM/runtime operation permitted in this test')))
            if not plots:stack.enter_context(patch('ascompare.concordance_report._plots'))
            summary=run_concordance(config,lock,entries,out)
        records=json.loads((out/'downstream_records.json').read_text())
        return summary,records

    def test_three_sample_public_workflow_has_separate_endpoints_and_synthetic_outputs(self):
        summary,records=self.run_results(plots=True)
        self.assertTrue(summary['synthetic']);self.assertEqual(summary['selected_samples'],3)
        self.assertEqual(summary['reference_version'],REFERENCE_VERSION)
        aa=[r for r in records if r['analysis_kind']=='aa_ac']
        ac=[r for r in records if r['analysis_kind']=='ac_only']
        self.assertEqual(len(aa),6);self.assertTrue(all(r['status']=='ok' and r['comparison_eligible'] for r in aa))
        self.assertTrue(all(r['status']=='missing' and not r['comparison_eligible'] for r in ac if r['version']==REFERENCE_VERSION))
        self.assertTrue(all(r['source_graph_identity'] for r in ac if r['version']==SOURCE_VERSION))
        multi=next(r for r in aa if r['sample_id']=='SYNTHETIC_03' and r['version']==REFERENCE_VERSION)
        self.assertEqual(len(multi['amplicons']),2);self.assertIn('FAN',multi['labels'])
        remote=next(r for r in aa if r['sample_id']=='SYNTHETIC_02' and r['version']==REFERENCE_VERSION)
        self.assertIn(['chr7',1000000,1020000],remote['reconstruction_intervals'])
        self.assertFalse(list(self.root.rglob('*.bam')))
        self.assertTrue(all(x['decision']=='insufficient_evidence' for x in summary['compatibility_combined']))
        for path in (self.root/'output').rglob('report.md'):
            text=path.read_text()
            self.assertTrue('SYNTHETIC' in text or '人工合成测试数据' in text,str(path))
        for path in (self.root/'output').glob('*.svg'):
            self.assertIn('SYNTHETIC TEST DATA',path.read_text())

    def test_no_bam_single_source_cannot_become_an_aa_version_comparison(self):
        self.change_manifest(lambda data:[s.update(comparison_results=[]) for s in data['samples']])
        summary,records=self.run_results()
        target=[r for r in records if r['version']==REFERENCE_VERSION]
        self.assertEqual(len(target),6)
        self.assertTrue(all(r['status']=='missing' and not r['comparison_eligible'] for r in target))
        with (self.root/'output/sample_similarity.tsv').open() as stream:
            rows=list(csv.DictReader(stream,delimiter='\t'))
        comparisons=[r for r in rows if r['target_version']==REFERENCE_VERSION]
        self.assertTrue(all(r['status']=='not_evaluable' and r['interval_jaccard']=='' for r in comparisons))
        self.assertTrue(all(r['decision']=='insufficient_evidence' for r in summary['compatibility_combined']))

    def test_original_result_inputs_remain_read_only_and_unchanged(self):
        files=[p for p in self.manifest.parent.rglob('*') if p.is_file()]
        before={p:digest(p) for p in files}
        for path in files:path.chmod(0o444)
        self.run_results()
        self.assertEqual(before,{p:digest(p) for p in files})
        self.assertEqual(set(files),{p for p in self.manifest.parent.rglob('*') if p.is_file()})

    def test_output_may_not_be_nested_in_original_results(self):
        original=self.manifest.parent/'results/SYNTHETIC_01'/SOURCE_VERSION
        output=original/'nested_output'
        with self.assertRaisesRegex(ValueError,'overlaps'):self.run_results(output)
        self.assertFalse(output.exists())

    def test_wrong_components_direct_version_conflicts_and_failed_sources_are_not_comparable(self):
        def edit(data):
            data['samples'][0]['comparison_results'][0]['component_versions']['ac']='9.9.9'
            data['samples'][1]['comparison_results'][0]['aa_version']='9.9.r9'
        self.change_manifest(edit)
        stages=self.manifest.parent/'results/SYNTHETIC_03'/SOURCE_VERSION/'stages.json'
        atomic_json(stages,{'aa':{'status':'complete'},'ac':{'status':'failed','reason':'SYNTHETIC old failure'}})
        _,records=self.run_results()
        lookup={(r['sample_id'],r['version'],r['analysis_kind']):r for r in records}
        self.assertFalse(lookup['SYNTHETIC_01',REFERENCE_VERSION,'aa_ac']['comparison_eligible'])
        self.assertFalse(lookup['SYNTHETIC_02',REFERENCE_VERSION,'aa_ac']['comparison_eligible'])
        old=lookup['SYNTHETIC_03',SOURCE_VERSION,'aa_ac']
        self.assertEqual(old['status'],'failed');self.assertFalse(old['comparison_eligible'])

    def test_manifest_full_suite_is_rejected_before_execution(self):
        self.change_manifest(lambda data:data.update(analysis_scope='full_suite'))
        with self.assertRaisesRegex(ValueError,'Whole-genome'):load_manifest(self.manifest)

    def test_patient_grouping_requires_complete_identifiers(self):
        self.change_manifest(lambda data:data['samples'][0].update(patient_id='SYNTHETIC_PATIENT'))
        with self.assertRaisesRegex(ValueError,'every sample'):
            load_manifest(self.manifest)
        self.change_manifest(lambda data:[s.update(patient_id='SYNTHETIC_PATIENT') for s in data['samples']])
        config, _, _ = load_manifest(self.manifest)
        self.assertEqual({s['patient_id'] for s in config['samples']}, {'SYNTHETIC_PATIENT'})

    def test_prepare_panel_uses_installed_entry_and_keeps_reference(self):
        out = self.root / 'prepare'
        with patch('ascompare.versions.prepare_sources', side_effect=lambda entry, destination: {'id': entry['id']}) as fetch:
            self.assertEqual(main(['prepare-panel', '--out', str(out), '--versions', SOURCE_VERSION]), 0)
        self.assertEqual(fetch.call_count, 2)
        panel = json.loads((out / 'prepared_panel.json').read_text())
        self.assertEqual([row['id'] for row in panel['sources']], [REFERENCE_VERSION, SOURCE_VERSION])

    def test_check_records_unrecognized_native_header_without_running_bam(self):
        bad_dir = self.manifest.parent / 'results/SYNTHETIC_01' / REFERENCE_VERSION / 'ac'
        profile = next(bad_dir.glob('*_amplicon_classification_profiles.tsv'))
        profile.write_text(profile.read_text().replace('amplicon_number', 'unsupported_header'))
        with patch('ascompare.concordance.run_bounded_downstream', side_effect=AssertionError('Check cannot run AA')):
            self.assertEqual(main(['check', '--manifest', str(self.manifest), '--out', str(self.root/'check')]), 2)
        rows = json.loads((self.root/'check'/'import_format_check.json').read_text())
        self.assertEqual(len(rows), 6)
        target = next(r for r in rows if r['version'] == REFERENCE_VERSION and r['sample_id'] == 'SYNTHETIC_01')
        self.assertNotEqual(target['status'], 'ok')
        self.assertIn('unsupported_header', target['profile_columns'])

    def test_cli_demo_uses_real_workflow_and_leaves_reusable_manifest(self):
        out=self.root/'cli_demo'
        with patch('ascompare.concordance_report._plots'):
            self.assertEqual(main(['demo','--out',str(out)]),0)
        proof=json.loads((out/'demo_provenance.json').read_text())
        self.assertFalse(proof['bam_used']);self.assertFalse(proof['AA_executed']);self.assertFalse(proof['AC_executed'])
        config,_,_=load_manifest(proof['manifest'])
        self.assertEqual(config['mode'],'results_only');self.assertEqual(config['runtimes'],[])


if __name__=='__main__':unittest.main()
