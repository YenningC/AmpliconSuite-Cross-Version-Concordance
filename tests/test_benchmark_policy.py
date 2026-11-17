import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam
from ascompare.benchmark_policy import guarded_command


class BoundedPolicyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.bam = self.root / 'SYNTHETIC.bam'
        with pysam.AlignmentFile(self.bam, 'wb', header={'HD': {'SO':'coordinate'}, 'SQ':[{'SN':'chr1','LN':1000000}]}) as f:
            for n in range(5):
                r = pysam.AlignedSegment()
                r.query_name = 'synthetic_' + str(n)
                r.query_sequence = 'A' * 50
                r.flag = 0; r.reference_id = 0; r.reference_start = 1000+n*100
                r.mapping_quality = 60; r.cigar = [(0,50)]
                r.query_qualities = pysam.qualitystring_to_array('I'*50)
                f.write(r)
        pysam.index(str(self.bam))
        self.policy = self.root / 'policy.json'
        self.policy.write_text(json.dumps(dict(alignments={str(self.bam):str(self.bam)+'.bai'},
            max_fetch_bp=10000,max_total_fetch_bp=20000,max_total_reads=100)))

    def tearDown(self):
        self.tmp.cleanup()

    def run_script(self, body):
        script = self.root / 'caller.py'
        script.write_text('import pysam\nbam=pysam.AlignmentFile('+repr(str(self.bam))+',"rb")\n'+body+'\n')
        return subprocess.run(guarded_command(sys.executable, script, [], self.policy), cwd=self.root,
                              capture_output=True, text=True)

    def test_finite_fetch_and_index_queries(self):
        r = self.run_script('assert sum(1 for _ in bam.fetch("chr1",900,2000))==5\nassert bam.mapped==5')
        self.assertEqual(r.returncode,0,r.stderr)
        audits = list(self.root.glob('bam_access.*.json'))
        report = json.loads(audits[0].read_text())
        self.assertEqual(report['queries'],1)
        self.assertEqual(report['yielded_reads'],5)

    def test_unbounded_fetch_and_iteration_and_whole_contig(self):
        for body in ('list(bam.fetch())', 'list(bam)', 'list(bam.fetch("chr1",0,1000000))',
                     'list(bam.fetch(until_eof=True))'):
            with self.subTest(body=body):
                r = self.run_script(body)
                self.assertNotEqual(r.returncode,0)
                self.assertIn('prohibited',r.stderr)

    def test_caught_forbidden_access_still_fails(self):
        r = self.run_script('try:\n list(bam.fetch())\nexcept RuntimeError:\n pass')
        self.assertNotEqual(r.returncode,0)
        self.assertIn('Caller caught',r.stderr)

    def test_budget_exhaustion(self):
        r = self.run_script('for _ in range(30):\n list(bam.fetch("chr1",900,2000))')
        self.assertNotEqual(r.returncode,0)
        self.assertIn('Cumulative',r.stderr)

    def test_no_index_creation_or_external_scan(self):
        for body in ('pysam.index('+repr(str(self.bam))+')',
                     'import subprocess\nsubprocess.run(["samtools","flagstat",'+repr(str(self.bam))+'])'):
            r = self.run_script(body)
            self.assertNotEqual(r.returncode,0)

    def test_bedcov_is_bounded(self):
        bed = self.root / 'targets.bed'; bed.write_text('chr1\t900\t2000\n')
        r = self.run_script('x=pysam.bedcov('+repr(str(bed))+','+repr(str(self.bam))+')\nassert int(x.split()[-1])==250')
        self.assertEqual(r.returncode,0,r.stderr)
        bed.write_text('chr1\t0\t1000000\n')
        r = self.run_script('pysam.bedcov('+repr(str(bed))+','+repr(str(self.bam))+')')
        self.assertNotEqual(r.returncode,0)

    def test_chunked_genome_tiling_is_forbidden(self):
        policy=json.loads(self.policy.read_text())
        policy.update(max_total_fetch_bp=2000000,max_fetch_bp=100000)
        self.policy.write_text(json.dumps(policy))
        r=self.run_script('for start in range(0,1000000,100000):\n list(bam.fetch("chr1",start,start+100000))')
        self.assertNotEqual(r.returncode,0)
        self.assertIn('genome tiling',r.stderr)


if __name__=='__main__': unittest.main()
