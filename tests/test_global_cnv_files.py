import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ascompare import concordance_global_cnv as wgs
from ascompare.versions import load_lock


class GlobalFilesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve() / 'case'
        self.bin = self.root / 'case_cnvkit_output'
        self.bin.mkdir(parents=True)
        self.cns = self.bin / 'tumor.cns'
        self.contents = 'chromosome\tstart\tend\tgene\tlog2\tprobes\nchr1\t0\t50\t-\t1\t2\n'
        self.cns.write_text(self.contents)
        self.tumor, self.normal = (self.root / n for n in ('tumor.bam', 'normal.bam'))
        for path in (self.tumor, self.normal):
            path.write_bytes(b'fixture')
            Path(str(path) + '.bai').write_bytes(b'index')
        self.fasta = self.root / 'ref.fa'
        self.fasta.write_text('>chr1\n' + 'A' * 50 + '\n')
        Path(str(self.fasta) + '.fai').write_text('chr1\t50\t6\t50\t51\n')
        self.sample = dict(
            sample_id='SYNTHETIC',
            reference_build='hg19',
            source_version='suite-1.6.2',
            source_result={'result_dir': str(self.root)},
            alignment=str(self.tumor),
            index=str(self.tumor) + '.bai',
            normal_alignment=str(self.normal),
            normal_index=str(self.normal) + '.bai',
        )
        self.config = dict(reference=dict(build='hg19', fasta=str(self.fasta), aa_data_repo=str(self.root)))
        self.metadata = self.root / 'case_run_metadata.json'
        self.meta = {
            'AmpliconSuite-pipeline_version': '1.6.2',
            'ref_genome': 'hg19',
            'cnvkit_version': '0.9.13',
            'cnvkit_cmd': f'python cnvkit.py batch {self.tumor} -m wgs --fasta {self.bin}/ref_reduced.fa -p 4 -d {self.bin} --normal {self.normal} ; python cnvkit.py segment {self.bin}/tumor.cnr -p 4 -m cbs -o {self.cns}',
        }
        self.metadata.write_text(json.dumps(self.meta))
        self.entry = next(e for e in load_lock()['versions'] if e['id'] == 'suite-1.6.2')
        self.runtime = dict(id='suite-1.6.2', cnvkit_version='0.9.13')

    def run_global(self):
        return wgs.global_cns(self.sample, self.runtime, self.entry, self.config, self.root / 'cache')

    def test_native_archive_requires_no_python(self):
        result = self.run_global()
        self.assertEqual(result['cns'], str(self.cns))
        self.assertEqual(result['evidence_level'], 'historical_metadata_only')
        self.assertTrue(Path(result['provenance']).exists())

    def test_archive_list_beyond_source(self):
        self.sample['global_cnv_results'] = [self.sample.pop('source_result')]
        self.assertEqual(self.run_global()['reuse'], 'archived_same_version_cns')

    def test_truncated_archive_cns_starts_from_complete_cnr(self):
        self.cns.write_text(self.contents.replace('50', '25'))
        (self.bin / 'tumor.cnr').write_text(self.contents)
        archive, rejected = wgs._historical_inputs(self.sample, self.config, '0.9.13')
        self.assertEqual(archive['kind'], 'cnr')
        self.assertIn('chromosome extents', archive['cns_rejection'])

    def test_malformed_cns_falls_back_to_native_ratios(self):
        self.cns.write_text('invalid header\n')
        (self.bin / 'tumor.cnr').write_text(self.contents)
        archive, _ = wgs._historical_inputs(self.sample, self.config, '0.9.13')
        self.assertEqual(archive['kind'], 'cnr')
        self.assertIn('Invalid archived CNS', archive['cns_rejection'])

    def test_checksum_mismatch_does_not_authorize_cns_import(self):
        self.sample['source_result']['cns_sha256'] = 'wrong'
        archive, rejected = wgs._historical_inputs(self.sample, self.config, '0.9.13')
        self.assertIsNone(archive)
        self.assertTrue(rejected)

    def test_cross_directory_import_verifies_prior_and_preserves_metadata(self):
        import sys

        self.metadata.unlink()
        self.sample.pop('normal_alignment')
        self.sample.pop('normal_index')
        refdir = self.root / 'hg19'
        refdir.mkdir()
        (refdir / 'hg19_cnvkit_filtered_ref.cnn').write_text('synthetic')
        self.runtime.update(cnvkit_python=sys.executable, cnvkit_script='cnvkit.py')
        runtime_id = {'version': '0.9.13'}

        def stub(cmd, stage, out, env, timeout):
            (out / 'tumor.cnr').write_text(self.contents)
            (out / 'tumor.cns').write_text(self.contents)

        previous = self.root / 'previous'
        with (
            patch.object(wgs, '_verified_probe', return_value=(runtime_id, self.root)),
            patch.object(wgs, '_native_cmd', side_effect=stub) as native,
        ):
            record = wgs.global_cns(
                self.sample, self.runtime, self.entry, self.config, previous / 'global_cnv_cache'
            )
            self.assertEqual(native.call_count, 2)
            receipt = Path(record['provenance'])
            original = receipt.read_bytes()
            self.config['reuse_global_cnv_from'] = str(previous)
            current = self.run_global()
            self.assertEqual(native.call_count, 2)
            self.assertEqual(current['reuse'], 'prior_same_recipe_cache')
            self.assertEqual(original, receipt.read_bytes())
            Path(record['cns']).write_text(self.contents.replace('50', '25'))
            self.assertEqual(Path(current['cns']).read_text(), self.contents)
            next_result = wgs.global_cns(
                self.sample, self.runtime, self.entry, self.config, self.root / 'third'
            )
            self.assertEqual(next_result['reuse'], 'native_global_batch')
            self.assertEqual(native.call_count, 4)
            del self.config['reuse_global_cnv_from']
            resumed = wgs.global_cns(self.sample, self.runtime, self.entry, self.config, previous / 'global_cnv_cache')
            self.assertEqual(resumed['reuse'], 'same_recipe_cnr_checkpoint')
            self.assertEqual(native.call_count, 5)

    def test_local_subset_extra_options_wrong_normal_and_rescaled_are_rejected(self):
        for before, after in [
            ('-m wgs', '-m wgs --targets local.bed'),
            ('-m wgs', '-m amplicon'),
            ('-m cbs', '-m haar'),
            (str(self.normal), str(self.tumor)),
            ('tumor.cns', 'tumor.call.cns'),
        ]:
            with self.subTest(value=after):
                self.metadata.write_text(
                    json.dumps(dict(self.meta, cnvkit_cmd=self.meta['cnvkit_cmd'].replace(before, after)))
                )
                archive, rejected = wgs._historical_inputs(self.sample, self.config, '0.9.13')
                self.assertIsNone(archive)
                self.assertTrue(rejected)

    def test_archived_cnr_can_start_at_segmentation(self):
        self.cns.unlink()
        cnr = self.bin / 'tumor.cnr.gz'
        with gzip.open(cnr, 'wt') as f:
            f.write(self.contents)
        with patch.object(wgs, '_fresh', return_value={'fixture': True}) as new:
            self.run_global()
        self.assertEqual(new.call_args.args[-2]['kind'], 'cnr')
        self.assertEqual(new.call_args.args[-2]['file'], str(cnr))

    def test_source_version_label_alone_is_insufficient(self):
        self.metadata.unlink()
        with patch.object(wgs, '_fresh', return_value={}) as new:
            self.run_global()
        self.assertIsNone(new.call_args.args[-2])
        self.assertTrue(new.call_args.args[-1])

    def test_malformed_or_empty_cns_is_not_a_negative(self):
        for contents in (self.contents.splitlines()[0] + '\n', self.contents.replace('50', '500')):
            self.cns.write_text(contents)
            with patch.object(wgs, '_fresh', return_value={}) as new:
                self.run_global()
            self.assertIsNone(new.call_args.args[-2])
            self.assertTrue(new.call_args.args[-1])

    def test_no_cross_version_cnr(self):
        self.meta['cnvkit_version'] = '0.9.10'
        self.metadata.write_text(json.dumps(self.meta))
        with patch.object(wgs, '_fresh', return_value={}) as new:
            self.run_global()
        self.assertIsNone(new.call_args.args[-2])

    def test_runtime_recipe_mismatch_is_an_explicit_failure(self):
        self.runtime['cnvkit_version'] = '0.9.14'
        with self.assertRaises(wgs.LocalFailure) as err:
            self.run_global()
        self.assertEqual(err.exception.stage, 'version_check')


class DepthCacheTest(unittest.TestCase):
    def test_exact_input_cache_bypasses_only_duplicate_native_call(self):
        import pysam
        from ascompare import concordance_global_depth as depth

        with (
            tempfile.TemporaryDirectory() as tmp,
            patch.object(depth, 'engine', return_value=['synthetic']),
            patch.object(pysam, 'bedcov', return_value='chr1\t0\t10\t42\n') as native,
        ):
            root = Path(tmp)
            bed = root / 'bins.bed'
            bed.write_text('chr1\t0\t10\n')
            spec = dict(
                cache=tmp,
                events=str(root / 'events'),
                alignments={'tumor.bam': 'sample1'},
                reference=None,
                environment='flags',
            )
            depth.install(spec)
            args = (str(bed), 'tumor.bam')
            self.assertEqual(pysam.bedcov(*args, split_lines=False), native.return_value)
            self.assertEqual(pysam.bedcov(*args, split_lines=False), native.return_value)
            self.assertEqual(native.call_count, 1)
            # New counting options, changed bin definitions and changed output
            # checksum all require recomputation.
            pysam.bedcov(*args, '-Q', '1', split_lines=False)
            bed.write_text('chr1\t0\t15\n')
            pysam.bedcov(*args, split_lines=False)
            for payload in root.rglob('*.gz'):
                payload.write_bytes(b'corrupt')
            pysam.bedcov(*args, split_lines=False)
            self.assertEqual(native.call_count, 4)
            self.assertEqual(
                [json.loads(x)['hit'] for x in (root / 'events').read_text().splitlines()],
                [False, True, False, False, False],
            )
