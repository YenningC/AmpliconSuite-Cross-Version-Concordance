"""Resolve the source version before choosing a local comparator runtime."""
import json
from pathlib import Path
import tempfile
import unittest

from ascompare.concordance_cnv_source import identify_source_cnvkit


class SourceVersionTests(unittest.TestCase):
    def sample(self, root, **values):
        return dict(sample_id='SYNTHETIC', source_result={'result_dir': str(root)},
                    source_version='suite-1.6.1', reference_build='hg19', **values)

    def test_original_metadata_can_differ_from_native_recipe(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'SYNTHETIC_run_metadata.json'
            p.write_text(json.dumps(dict(cnvkit_version='0.9.14', cnvkit_cmd='cnvkit.py batch SYNTHETIC.bam',
                                        ref_genome='hg19', **{'AmpliconSuite-pipeline_version': '1.6.1'})))
            sample = self.sample(tmp)
            identify_source_cnvkit(sample)
            self.assertEqual(sample['source_cnvkit_version'], '0.9.14')
            self.assertEqual(sample['source_cnvkit_evidence']['origin'], 'source_result_run_metadata')
            self.assertEqual(len(sample['source_cnvkit_evidence']['sha256']), 64)
            conflict = self.sample(tmp, source_cnvkit_version='0.9.13')
            identify_source_cnvkit(conflict)
            self.assertIsNone(conflict['source_cnvkit_version'])
            self.assertIn('conflicts', conflict['_source_cnvkit_issue'])

    def test_skipped_or_ambiguous_metadata_does_not_assign_a_version(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'SYNTHETIC_run_metadata.json'
            p.write_text(json.dumps(dict(cnvkit_version='0.9.13', cnvkit_cmd='NA')))
            sample = self.sample(tmp)
            identify_source_cnvkit(sample)
            self.assertIsNone(sample['source_cnvkit_version'])
            specified = self.sample(tmp, source_cnvkit_version='0.9.14')
            identify_source_cnvkit(specified)
            self.assertEqual(specified['source_cnvkit_evidence']['origin'], 'sample_declaration')
            p.with_name('second_run_metadata.json').write_text(p.read_text())
            sample = self.sample(tmp)
            identify_source_cnvkit(sample)
            self.assertIn('Multiple', sample['_source_cnvkit_issue'])
