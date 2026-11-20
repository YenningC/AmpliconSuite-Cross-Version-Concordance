"""Stage-control fixtures. Native multi-release integration is run separately."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from ascompare.common import digest
from ascompare.concordance_components import _seed_bed_text
from ascompare.concordance_cnv import LocalFailure
from ascompare import concordance_global_seeds as global_seeds
from ascompare.versions import _digest


def locked(root, entrypoint, files):
    root.mkdir()
    for name, text in files.items():
        target = root / name
        target.parent.mkdir(exist_ok=True, parents=True)
        target.write_text(text)
    hashes = {f: digest(root / f) for f in files}
    manifest = dict(commit='a' * 40, repo=root.name, files=hashes)
    (root / '.ascompare-source.json').write_text(json.dumps(manifest))
    return dict(repo=root.name, entrypoint=entrypoint, commit=manifest['commit'], manifest_sha256=_digest(hashes),
                version='synthetic')


class GlobalSeeds(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.tmp = Path(temp.name).resolve()
        ref = self.tmp / 'ref' / 'hg19'
        ref.mkdir(parents=True)
        (ref / 'hg19_noAlt.fa.fai').write_text('chr1\t900000\n')
        self.bed_text = 'chr1\t400\t500\t6.3\tinput.bed\nchr1\t200\t300\t7.0\tinput.bed\n'
        self.cns = self.tmp / 'tumor.cns'
        self.cns.write_text('chromosome\tstart\tend\tgene\tlog2\nchr1\t0\t900000\t-\t0\n')
        suite, aa = self.tmp / 'suite', self.tmp / 'AA'
        self.entry = dict(
            aa=locked(aa, 'src/AA.py', {'src/AA.py': '', 'src/amplified_intervals.py': '# synthetic'}),
            suite=locked(suite, 'suite.py', {
                'suite.py': """
def get_ref_sizes(path):
    assert open(path).read() == 'chr1\\t900000\\n'
    return {'chr1': '899999'}
def get_ref_centromeres(build):
    assert AA_REPO and build == 'hg19'
    return {}
def convert_cnvkit_cns_to_bed(path, base, cnsfile=None):
    assert open(cnsfile).readline().endswith('log2\\n')
    name = path + base + '_CNV_CALLS.bed'
    open(name, 'w').write('chr1\\t0\\t900000\\tCNVkit\\t2.0\\n')
    return name
""", 'paalib/cnv_prefilter.py': """
def prefilter_bed(bed, build, cents, sizes, gain, out):
    assert sizes == {'chr1':'899999'} and gain == 4.5 and cents == {}
    open('filtered.bed','w').write(open(bed).read())
    return 'filtered.bed'
"""}))
        self.runtime = dict(aa_script=str(aa / 'src/AA.py'), aa_python=sys.executable, ac_python=sys.executable,
                            suite_script=str(suite / 'suite.py'), id='synthetic')
        from ascompare.concordance import fast_fingerprint
        bam, index = self.tmp / 'synthetic.bam', self.tmp / 'synthetic.bam.bai'
        bam.touch()
        index.touch()
        self.sample = dict(sample_id='test', reference_build='hg19', alignment=str(bam), index=str(index),
                           _input_identities=[fast_fingerprint(bam), dict(path=str(index), sha256=digest(index))])
        policy = self.tmp / 'policy.json'
        policy.write_text('{}')
        info = dict(packages=dict(returncode=0, stdout='synthetic==1'))
        self.config = dict(reference={'aa_data_repo': str(ref.parent)}, parameter_policy='suite_defaults',
                           _reference_lengths={'chr1': 900000}, _io_policy_path=str(policy),
                           _reference_identity=[dict(test=True)], _seed_verified_runtime=dict(aa=info, ac=info))
        self.mock = self.enter_patch(patch.object(global_seeds, '_run_tool', side_effect=self.simulate_filter))

    def enter_patch(self, manager):
        obj = manager.start()
        self.addCleanup(manager.stop)
        return obj

    def simulate_filter(self, interpreter, script, args, policy, cwd, env, stage, timeout):
        self.assertEqual(script.name, 'amplified_intervals.py')
        self.assertEqual(args[args.index('--bed') + 1], 'filtered.bed')
        self.assertNotIn('--no_cstats', args)
        (cwd / global_seeds.FINAL_BED).write_text(self.bed_text)
        p = cwd / 'seed_processing.log'
        p.write_text('success')
        return p

    def run_attempt(self, n):
        return global_seeds.global_suite_seeds(self.sample, self.runtime, self.entry, self.config, self.cns,
                                               self.tmp / ('attempt_%04d' % n))

    def test_native_order_and_fields_are_identity_and_resume_checked(self):
        seeds = self.run_attempt(1)
        self.assertEqual(_seed_bed_text(seeds), self.bed_text)
        self.assertEqual(self.run_attempt(2), seeds)
        self.assertEqual(self.mock.call_count, 1)
        # Tampering with the producer filter log invalidates the whole filter.
        (self.tmp / 'attempt_0001/seed_processing.log').write_text('edited')
        self.assertEqual(self.run_attempt(3), seeds)
        self.assertEqual(self.mock.call_count, 2)
        self.assertNotEqual(_seed_bed_text(list(reversed(seeds))), _seed_bed_text(seeds))
        other = [dict(s, source_fields=['different']) for s in seeds]
        self.assertNotEqual(_seed_bed_text(other), _seed_bed_text(seeds))

    def test_empty_is_evaluable_only_after_success(self):
        self.bed_text = ''
        self.assertEqual(self.run_attempt(1), [])
        self.assertEqual(self.run_attempt(2), [])
        self.cns.write_text(self.cns.read_text().replace('900000', '800000'))
        self.mock.side_effect = LocalFailure('failed', 'seed_processing', 'filter failed')
        with self.assertRaisesRegex(LocalFailure, 'filter failed'):
            self.run_attempt(3)

    def test_mutated_bam_rejects_cached_seed_filter(self):
        self.run_attempt(1)
        Path(self.sample['alignment']).write_text('changed')
        with self.assertRaisesRegex(ValueError, 'BAM/index changed'):
            self.run_attempt(2)
        self.assertEqual(self.mock.call_count, 1)

    def test_reject_called_malformed_or_empty_cns_before_filter(self):
        self.cns = self.tmp / 'sample.call.cns'
        self.cns.touch()
        with self.assertRaisesRegex(ValueError, 'uncalled'):
            self.run_attempt(1)
        self.cns = self.tmp / 'tumor.cns'
        self.cns.write_text('chromosome\tstart\tend\tgene\tlog2\n')
        with self.assertRaisesRegex(ValueError, 'no segments'):
            self.run_attempt(2)
        self.cns.write_text('chromosome\tstart\tend\tlog2\tgene\nchr1\t0\t10\t0\t-\n')
        with self.assertRaisesRegex(ValueError, 'order'):
            self.run_attempt(3)
        self.mock.assert_not_called()

    def test_missing_bed_after_success_and_changed_source_are_detected(self):
        def no_bed(*args):
            log = args[4] / 'seed_processing.log'
            log.write_text('success')
            return log
        self.mock.side_effect = no_bed
        with self.assertRaisesRegex(LocalFailure, 'did not produce'):
            self.run_attempt(1)
        Path(self.runtime['suite_script']).write_text('changed')
        with self.assertRaisesRegex(ValueError, 'changed source'):
            self.run_attempt(2)
