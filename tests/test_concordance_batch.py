"""Batch import validation and public CLI workflow on invented results only."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ascompare.concordance import load_manifest, main
from ascompare.concordance_batch import import_samples
from ascompare.concordance_demo import write_demo_manifest


class ConcordanceBatchTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.tables = self.root / 'tables'
        self.settings = self.root / 'settings'
        self.tables.mkdir(); self.settings.mkdir()
        self.tsv = self.tables / 'samples.tsv'
        self.config = self.settings / 'shared.json'
        self.output = self.root / 'prepared' / 'manifest.json'
        self.base = dict(mode='results_only', target_versions=['suite-1.6.0', 'suite-1.6.2'],
                         statistics=dict(primary_pairs=[], min_samples=20, margins={}))
        self.config.write_text(json.dumps(self.base))
        self.header = ['sample_id', 'patient_id', 'source_version', 'reference_build',
                       'result_version', 'aa_dir', 'ac_dir']

    def tearDown(self):
        self.temporary.cleanup()

    def write_rows(self, rows, header=None):
        with self.tsv.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(header or self.header)
            writer.writerows(rows)

    def row(self, sample='SYNTHETIC_01', patient='SYNTHETIC_PATIENT_01',
            source='suite-1.6.0', version='', aa='source/aa', ac='source/ac'):
        return [sample, patient, source, 'hg19', version, aa, ac]

    def test_repeated_version_rows_merge_and_preserve_sample_and_patient_order(self):
        self.write_rows([
            self.row('SYNTHETIC_02', 'SYNTHETIC_PATIENT_01', version='suite-1.6.2', aa='target/aa', ac='target/ac'),
            self.row('SYNTHETIC_01', 'SYNTHETIC_PATIENT_01'),
            self.row('SYNTHETIC_02', 'SYNTHETIC_PATIENT_01'),
        ])
        result = import_samples(self.tsv, self.config, self.output)
        self.assertEqual([sample['sample_id'] for sample in result['samples']], ['SYNTHETIC_02', 'SYNTHETIC_01'])
        sample = result['samples'][0]
        self.assertEqual(sample['patient_id'], 'SYNTHETIC_PATIENT_01')
        self.assertEqual(sample['source_result']['aa_dir'], str(self.tables / 'source/aa'))
        self.assertEqual(sample['comparison_results'][0]['version'], 'suite-1.6.2')
        self.assertEqual(result['batch_import']['result_row_count'], 3)
        self.assertEqual(result['batch_import']['sample_count'], 2)
        self.assertEqual(len(result['batch_import']['samples_tsv']['sha256']), 64)
        self.assertEqual(load_manifest(self.output)[0], result)

    def test_config_paths_and_tsv_paths_have_independent_bases_without_reading_bam(self):
        self.base.update(mode='with_bam', cnv_scope='local', reference=dict(build='hg19', fasta='reference/hg19.fa',
                         aa_data_repo='reference/AA_DATA_REPO'),
                         local_cnv=dict(reference_mode='flat', background_bed='background.bed'),
                         runtimes=[dict(id='suite-1.6.2', aa_python='env/bin/python',
                                        aa_script='sources/AA.py', extra_path=['other/bin'])])
        self.config.write_text(json.dumps(self.base))
        header = ['sample_id', 'source_version', 'reference_build', 'alignment', 'index',
                  'seeds_bed', 'run_dir', 'source_cnvkit_version', 'source_seed_bed',
                  'baseline_cns', 'baseline_provenance', 'aa_version', 'ac_version']
        self.write_rows([['SYNTHETIC_01', 'suite-1.6.0', 'hg19', 'tumor.bam', 'tumor.bam.bai',
            'seeds.bed', 'results', '0.9.13', 'original_seeds.bed', 'original.cns',
            'original.json', '1.6.r0', '2.0.0']], header)
        # None of these placeholder alignment/resource files exist. Conversion
        # must not open them or launch a scientific runtime.
        with patch('ascompare.concordance._validate_bam', side_effect=AssertionError('BAM opened')):
            result = import_samples(self.tsv, self.config, self.output)
        sample = result['samples'][0]
        self.assertEqual(sample['alignment'], str(self.tables / 'tumor.bam'))
        self.assertEqual(sample['baseline_provenance'], str(self.tables / 'original.json'))
        self.assertEqual(sample['source_result']['aa_version'], '1.6.r0')
        self.assertEqual(sample['seeds_bed'], str(self.tables / 'seeds.bed'))
        self.assertEqual(result['reference']['fasta'], str(self.settings / 'reference/hg19.fa'))
        self.assertEqual(result['local_cnv']['background_bed'], str(self.settings / 'background.bed'))
        self.assertEqual(result['runtimes'][0]['extra_path'], [str(self.settings / 'other/bin')])
        self.assertEqual(result['runtimes'][0]['aa_python'], str(self.settings / 'env/bin/python'))

    def test_suite_result_directory_discovery_and_utf8_bom(self):
        root = self.tables / 'suite_result'
        (root / 'SYNTHETIC_01_AA_results').mkdir(parents=True)
        (root / 'SYNTHETIC_01_classification').mkdir()
        self.tsv.write_text('\ufeffsample_id\tsource_version\treference_build\tresult_dir\n'
            'SYNTHETIC_01\tsuite-1.6.0\thg19\tsuite_result\n', encoding='utf-8')
        result = import_samples(self.tsv, self.config, self.output)
        self.assertEqual(result['samples'][0]['source_result']['aa_dir'], str(root / 'SYNTHETIC_01_AA_results'))

    def test_unknown_duplicate_or_missing_headers_and_ragged_rows_are_rejected(self):
        cases = [
            ('sample_id\tsource_version\treference_build\tunexpected\nSYNTHETIC_01\tsuite-1.6.0\thg19\tx\n', 'Unknown'),
            ('sample_id\tsource_version\tsource_version\n', 'duplicate'),
            ('sample_id\tsource_version\n', 'needs columns'),
            ('sample_id\tsource_version\treference_build\nSYNTHETIC_01\tsuite-1.6.0\n', 'fields'),
            ('sample_id\tsource_version\treference_build\nSYNTHETIC_01\t\thg19\n', 'required'),
            ('sample_id\tsource_version\treference_build\n', 'no sample rows'),
        ]
        for text, message in cases:
            with self.subTest(message=message):
                self.tsv.write_text(text)
                with self.assertRaisesRegex(ValueError, message):
                    import_samples(self.tsv, self.config, self.output)
                self.assertFalse(self.output.exists())

    def test_sample_conflicts_duplicate_versions_and_missing_source_rows_are_rejected(self):
        cases = [
            ([self.row(), self.row()], 'duplicate result_version'),
            ([self.row(), self.row(patient='OTHER_PATIENT', version='suite-1.6.2')], 'conflicting patient_id'),
            ([self.row(), self.row(source='suite-1.6.1', version='suite-1.6.2')], 'conflicting source_version'),
            ([self.row(version='suite-1.6.2')], 'no result row'),
            ([self.row(aa='', ac='')], 'provide result_dir'),
            ([self.row(), self.row('SYNTHETIC_02', patient='')], 'every sample'),
        ]
        for rows, message in cases:
            with self.subTest(message=message):
                self.write_rows(rows)
                with self.assertRaisesRegex(ValueError, message):
                    import_samples(self.tsv, self.config, self.output)
                self.assertFalse(self.output.exists())

    def test_invalid_versions_and_global_settings_use_normal_manifest_validation(self):
        for rows, change, message in [
            ([self.row(source='suite-99.0')], {}, 'source lock'),
            ([self.row(), self.row(version='suite-99.0')], {}, 'absent from lock'),
            ([self.row()], {'reference_version': 'suite-1.6.1'}, 'fixes its reference'),
            ([self.row()], {'analysis_scope': 'full_suite'}, 'Whole-genome'),
            ([self.row()], {'samples': [{'sample_id': 'MUST_NOT_MERGE'}]}, 'omit samples'),
            ([self.row()], {'threads': 2}, 'serially'),
        ]:
            with self.subTest(message=message):
                self.write_rows(rows)
                self.config.write_text(json.dumps({**self.base, **change}))
                with self.assertRaisesRegex(ValueError, message):
                    import_samples(self.tsv, self.config, self.output)
                self.assertFalse(self.output.exists())

    def test_with_bam_requires_explicit_source_seed_and_complete_normal_pair(self):
        self.base.update(mode='with_bam', cnv_scope='local', reference=dict(build='hg19', fasta='ref.fa', aa_data_repo='resources'))
        self.config.write_text(json.dumps(self.base))
        self.write_rows([self.row()])
        with self.assertRaisesRegex(ValueError, 'alignment, index, seeds_bed'):
            import_samples(self.tsv, self.config, self.output)
        header = ['sample_id', 'source_version', 'reference_build', 'alignment', 'index',
                  'seeds_bed', 'run_dir', 'normal_alignment', 'normal_index']
        self.write_rows([['SYNTHETIC_01', 'suite-1.6.0', 'hg19', 't.bam', 't.bam.bai',
                         'seed.bed', 'result', 'n.bam', '']], header)
        with self.assertRaisesRegex(ValueError, 'supplied together'):
            import_samples(self.tsv, self.config, self.output)

    def test_output_and_original_inputs_are_never_overwritten(self):
        self.write_rows([self.row()])
        before = {path: path.read_bytes() for path in (self.tsv, self.config)}
        import_samples(self.tsv, self.config, self.output)
        saved = self.output.read_bytes()
        self.assertEqual(main(['import-samples', '--samples', str(self.tsv), '--config', str(self.config),
                               '--out', str(self.output)]), 2)
        self.assertEqual(self.output.read_bytes(), saved)
        self.assertEqual(before, {path: path.read_bytes() for path in before})
        with self.assertRaisesRegex(ValueError, 'already exists'):
            import_samples(self.tsv, self.config, self.config)

    def test_imported_three_sample_manifest_runs_public_check_and_run_without_bam(self):
        fixture_manifest = write_demo_manifest(self.tables / 'SYNTHETIC_INPUTS')
        fixture = json.loads(fixture_manifest.read_text())
        config = {key: value for key, value in fixture.items() if key != 'samples'}
        self.config.write_text(json.dumps(config))
        header = ['sample_id', 'patient_id', 'source_version', 'reference_build',
                  'result_version', 'run_dir', 'aa_version', 'ac_version']
        rows = []
        for n, sample in enumerate(fixture['samples']):
            for version, spec in [(sample['source_version'], sample['source_result']),
                                  *[(r['version'], r) for r in sample['comparison_results']]]:
                rows.append([sample['sample_id'], f'SYNTHETIC_PATIENT_{n // 2 + 1:02d}',
                    sample['source_version'], sample['reference_build'], version,
                    str(Path('SYNTHETIC_INPUTS') / spec['run_dir']),
                    spec['component_versions']['aa'], spec['component_versions']['ac']])
        self.write_rows(rows, header)
        before = self.tsv.read_bytes(), self.config.read_bytes()
        self.assertEqual(main(['import-samples', '--samples', str(self.tsv), '--config', str(self.config),
                               '--out', str(self.output)]), 0)
        with patch('ascompare.concordance._validate_bam', side_effect=AssertionError('BAM must not be opened')):
            self.assertEqual(main(['check', '--manifest', str(self.output), '--out', str(self.root / 'check')]), 0)
            with patch('ascompare.concordance_report._plots'):
                self.assertEqual(main(['run', '--manifest', str(self.output), '--out', str(self.root / 'run')]), 0)
        summary = json.loads((self.root / 'run/summary.json').read_text())
        self.assertEqual(summary['selected_samples'], 3)
        self.assertTrue(summary['synthetic'])
        records = json.loads((self.root / 'run/downstream_records.json').read_text())
        existing = [row for row in records if row['analysis_kind'] == 'aa_ac']
        self.assertEqual(len(existing), 6)
        self.assertTrue(all(row['comparison_eligible'] for row in existing))
        self.assertEqual(before, (self.tsv.read_bytes(), self.config.read_bytes()))

    def test_invalid_middle_sample_does_not_prevent_later_sample_or_hide_progress(self):
        fixture_manifest = write_demo_manifest(self.tables / 'SYNTHETIC_INPUTS')
        fixture = json.loads(fixture_manifest.read_text())
        self.config.write_text(json.dumps({key: value for key, value in fixture.items() if key != 'samples'}))
        header = ['sample_id', 'source_version', 'reference_build', 'result_version', 'run_dir']
        rows = []
        for sample in fixture['samples']:
            for version, spec in [(sample['source_version'], sample['source_result']),
                                  *[(r['version'], r) for r in sample['comparison_results']]]:
                rows.append([sample['sample_id'], sample['source_version'], sample['reference_build'],
                             version, str(Path('SYNTHETIC_INPUTS') / spec['run_dir'])])
        self.write_rows(rows, header)
        import_samples(self.tsv, self.config, self.output)
        source = self.tables / 'SYNTHETIC_INPUTS/results/SYNTHETIC_02/suite-1.6.0/aa'
        next(source.glob('*_summary.txt')).unlink()  # An intentionally invalid synthetic source result.
        out = self.root / 'run_with_invalid_middle'
        with patch('ascompare.concordance._validate_bam', side_effect=AssertionError('BAM must not be opened')):
            check = self.root / 'check_with_invalid_middle'
            self.assertEqual(main(['check', '--manifest', str(self.output), '--out', str(check)]), 2)
            progress = json.loads((check / 'sample_progress.json').read_text())
            self.assertEqual([row['state'] for row in progress], ['ready', 'input_failed', 'ready'])
            from ascompare.concordance_import import import_result
            import_order = []
            def tracked_import(spec, sample_id, *args, **kwargs):
                import_order.append(sample_id)
                return import_result(spec, sample_id, *args, **kwargs)
            with patch('ascompare.concordance_report._plots'), \
                 patch('ascompare.concordance_import.import_result', side_effect=tracked_import):
                self.assertEqual(main(['run', '--manifest', str(self.output), '--out', str(out)]), 0)
        self.assertEqual(import_order, ['SYNTHETIC_01', 'SYNTHETIC_01', 'SYNTHETIC_03', 'SYNTHETIC_03'])
        progress = json.loads((out / 'sample_progress.json').read_text())
        self.assertEqual([row['sample_id'] for row in progress], ['SYNTHETIC_01', 'SYNTHETIC_02', 'SYNTHETIC_03'])
        self.assertEqual([row['order'] for row in progress], [1, 2, 3])
        self.assertEqual([row['state'] for row in progress], ['finished', 'input_failed', 'finished'])
        with (out / 'sample_progress.tsv').open(newline='') as handle:
            progress_tsv = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual([row['sample_id'] for row in progress_tsv], [row['sample_id'] for row in progress])
        self.assertEqual([row['state'] for row in progress_tsv], [row['state'] for row in progress])
        self.assertTrue(all(row['updated_utc'] for row in progress))
        records = json.loads((out / 'downstream_records.json').read_text())
        invalid = [row for row in records if row['sample_id'] == 'SYNTHETIC_02']
        later = [row for row in records if row['sample_id'] == 'SYNTHETIC_03' and row['analysis_kind'] == 'aa_ac']
        self.assertTrue(invalid)
        self.assertTrue(all(row['status'] == 'failed' and row['stage'] == 'input_check' for row in invalid))
        self.assertEqual(len(later), 2)
        self.assertTrue(all(row['status'] == 'ok' and row['comparison_eligible'] for row in later))


if __name__ == '__main__':
    unittest.main()
