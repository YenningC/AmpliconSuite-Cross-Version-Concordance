"""Interrupted synthetic batches retain a truthful, resumable sample ledger."""
from contextlib import ExitStack
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ascompare.common import atomic_json, digest, table
from ascompare.concordance import main
from ascompare.concordance_demo import write_demo_manifest


class ConcordanceProgressTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name).resolve()
        self.manifest = write_demo_manifest(self.root / 'SYNTHETIC_INPUTS')
        self.out = self.root / 'output'
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.addCleanup(self.temp.cleanup)
        for target in ('ascompare.concordance._validate_bam',
                       'ascompare.concordance.run_bounded_downstream',
                       'ascompare.concordance.verify_runtime',
                       'ascompare.concordance_cnv.run_local_upstream',
                       'ascompare.concordance_import.reclassify_result'):
            self.stack.enter_context(patch(target, side_effect=AssertionError('No BAM/runtime operation permitted')))
        self.stack.enter_context(patch('ascompare.concordance_report._plots'))
        self.stack.enter_context(patch('ascompare.concordance_cohort_plot.plot_cohort_compatibility'))

    def run_public(self):
        return main(['run', '--manifest', str(self.manifest), '--out', str(self.out)])

    def progress(self):
        progress = json.loads((self.out / 'sample_progress.json').read_text())
        with (self.out / 'sample_progress.tsv').open(newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual([r['sample_id'] for r in rows], [r['sample_id'] for r in progress])
        self.assertEqual([r['state'] for r in rows], [r['state'] for r in progress])
        return progress

    def test_keyboard_interrupt_marks_active_sample_and_same_manifest_resumes(self):
        inputs = {p: digest(p) for p in self.manifest.parent.rglob('*') if p.is_file()}
        with patch('ascompare.concordance_import.import_result', side_effect=KeyboardInterrupt):
            self.assertEqual(self.run_public(), 130)
        # All inputs are preflighted first; ready means not yet analyzed.
        progress = self.progress()
        self.assertEqual([r['state'] for r in progress], ['interrupted', 'ready', 'ready'])
        self.assertEqual([r['order'] for r in progress], [1, 2, 3])
        self.assertTrue(all(r['updated_utc'] for r in progress))
        self.assertEqual(self.run_public(), 0)
        self.assertEqual([r['state'] for r in self.progress()], ['finished'] * 3)
        records = json.loads((self.out / 'downstream_records.json').read_text())
        imported = [r for r in records if r['analysis_kind'] == 'aa_ac']
        self.assertEqual(len(imported), 6)
        self.assertTrue(all(r['status'] == 'ok' for r in imported))
        self.assertEqual(inputs, {p: digest(p) for p in inputs})

    def test_unexpected_processing_error_marks_only_active_sample_aborted(self):
        with patch('ascompare.concordance_import.import_result', side_effect=RuntimeError('SYNTHETIC processing failure')):
            self.assertEqual(self.run_public(), 2)
        self.assertEqual([r['state'] for r in self.progress()], ['aborted', 'ready', 'ready'])

    def test_interrupt_after_first_sample_preserves_checkpoint_and_old_artifacts(self):
        from ascompare.concordance_import import import_result

        def interrupt_second(spec, sample_id, *args, **kwargs):
            if sample_id == 'SYNTHETIC_02':
                raise KeyboardInterrupt
            return import_result(spec, sample_id, *args, **kwargs)

        with patch('ascompare.concordance_import.import_result', side_effect=interrupt_second):
            self.assertEqual(self.run_public(), 130)
        self.assertEqual([r['state'] for r in self.progress()], ['finished', 'interrupted', 'ready'])
        checkpoint = json.loads((self.out / 'downstream_records.json').read_text())
        self.assertEqual({r['sample_id'] for r in checkpoint}, {'SYNTHETIC_01'})
        originals = {p: digest(p) for p in (self.out / 'imports/SYNTHETIC_01').rglob('*') if p.is_file()}
        self.assertTrue(originals)
        self.assertEqual(self.run_public(), 0)
        self.assertEqual([r['state'] for r in self.progress()], ['finished'] * 3)
        self.assertEqual(originals, {p: digest(p) for p in originals})

    def test_rejected_changed_manifest_does_not_rewrite_prior_progress(self):
        with patch('ascompare.concordance_import.import_result', side_effect=KeyboardInterrupt):
            self.assertEqual(self.run_public(), 130)
        # Simulate a stale running marker left by an uncatchable process exit.
        # Rejection of a different manifest must not rewrite another run's log.
        progress = self.progress()
        progress[0]['state'] = 'running'
        atomic_json(self.out / 'sample_progress.json', progress)
        table(self.out / 'sample_progress.tsv', progress)
        paths = [self.out / name for name in ('sample_progress.json', 'sample_progress.tsv', 'provenance.json')]
        before = {p: p.read_bytes() for p in paths}
        data = json.loads(self.manifest.read_text())
        data['aa_parameters'] = {'downsample': 10}
        atomic_json(self.manifest, data)
        self.assertEqual(self.run_public(), 2)
        self.assertEqual(before, {p: p.read_bytes() for p in paths})


if __name__ == '__main__':
    unittest.main()
