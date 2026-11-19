import csv
import gzip
import io
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path

from scripts.compare_cnv_profiles import main


class CompareCnvProfilesTests(unittest.TestCase):
    def setUp(self):
        work = tempfile.TemporaryDirectory()
        self.addCleanup(work.cleanup)
        self.root = Path(work.name)
        self.header = 'chromosome\tstart\tend\tlog2\tweight\n'
        self.source = 'chr1\t0\t3\t0\t1\nchr1\t5\t10\t2\t1\nchr1\t10\t11\t8\t0\n'
        self.local = 'chr1\t1\t7\t1\t1\nchr1\t7\t11\t0\t1\n'
        self.bed = self.root/'panel.bed'
        self.bed.write_text('chr1\t0\t11\nchr2\t1\t7\n')
        self.args = ['--seeds-bed', str(self.bed), '--gain-cn', '4', '--out', str(self.root/'report.tsv')]
        for label in ['source', 'local', 'seed']:
            self.args.extend(['--'+label+'-build', 'hg19'])
        for label, stage, data in [('source','cnr',self.source), ('local','cnr',self.local),
                                   ('source','cns','chr1\t0\t11\t2\t1\n'),
                                   ('local','cns','chr1\t0\t11\t1\t1\n')]:
            gz = label == 'source' and stage == 'cnr'
            name = self.root/(label+'.'+stage+('.gz' if gz else ''))
            if gz:
                with gzip.open(name, 'wt') as handle:
                    handle.write(self.header+data)
            else:
                name.write_text(self.header+data)
            self.args += ['--'+label+'-'+stage,str(name)]

    def test_unequal_bins_gaps_and_cns_are_measured_on_the_same_paired_bases(self):
        main(self.args)
        with (self.root/'report.tsv').open() as h:
            cnr, cns, empty_cnr, empty_cns = csv.DictReader(h, delimiter='\t')
        self.assertEqual(int(cnr['paired_bp']), 7)
        self.assertEqual(int(cns['paired_bp']), 7)
        self.assertAlmostEqual(float(cnr['source_mean_log2']), 10/7)
        self.assertAlmostEqual(float(cnr['local_mean_log2']), 4/7)
        self.assertAlmostEqual(float(cnr['paired_mean_delta_log2']), -6/7)
        self.assertAlmostEqual(float(cnr['source_gain_bp_fraction']), 5/7)
        self.assertAlmostEqual(float(cnr['local_gain_bp_fraction']), 4/7)
        self.assertEqual(float(cns['paired_mean_delta_log2']), -1)
        for row in [empty_cnr, empty_cns]:
            self.assertEqual(row['status'], 'no_paired_bp')
            self.assertEqual(row['local_mean_log2'], '')

    def test_overlap_and_invalid_log_rejected_before_output(self):
        segment = self.root/'source.cns'
        for extra, reason in [('chr1\t3\t5\t1\t1','overlapping'), ('chr2\t5\t9\tnan\t1','invalid')]:
            with self.subTest(extra=extra), redirect_stderr(io.StringIO()) as stderr, self.assertRaises(SystemExit):
                segment.write_text(self.header+'chr1\t0\t11\t1\t1\n'+extra+'\n')
                main(self.args)
            self.assertIn(reason, stderr.getvalue())
            self.assertFalse((self.root/'report.tsv').exists())

    def test_build_mismatch_and_existing_output(self):
        with redirect_stderr(io.StringIO()) as stderr, self.assertRaises(SystemExit):
            main(self.args + ['--local-build','hg38'])
        self.assertIn('must match exactly', stderr.getvalue())
        out = self.root/'report.tsv'
        out.write_text('existing')
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            main(self.args)
        self.assertEqual(out.read_text(), 'existing')
