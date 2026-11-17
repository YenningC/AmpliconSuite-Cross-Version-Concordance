import tempfile
import unittest
from pathlib import Path
from ascompare.inputs import read_seeds,validate,validate_full
from ascompare.runner import prepare_alignment
import pysam

class InputTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.ref=self.root/'repo'/'hg19';self.ref.mkdir(parents=True)
        self.fasta=self.ref/'toy.fa';self.fasta.write_text('>chr1\n'+'ACGT'*250+'\n');pysam.faidx(str(self.fasta))
        (self.ref/'file_list.txt').write_text('fa_file toy.fa\nchrLen_file toy.fa.fai\n')
        self.bed=self.root/'seeds.bed';self.bed.write_text('chr1\t1\t100\tseed1\n')
        self.config={'reference':{'build':'hg19','fasta':str(self.fasta),'aa_data_repo':str(self.ref.parent)}}
    def bam(self,extra=False):
        p=self.root/'test.bam';sq=[{'SN':'chr1','LN':1000}]
        if extra:sq.append({'SN':'virus','LN':100})
        with pysam.AlignmentFile(str(p),'wb',header={'HD':{'SO':'coordinate'},'SQ':sq}) as f:
            r=pysam.AlignedSegment();r.query_name='test';r.query_sequence='A'*20;r.reference_id=0;r.reference_start=10;r.cigar=[(0,20)];r.mapping_quality=60;f.write(r)
        pysam.index(str(p));return p
    def row(self,bam):
        return {'sample_id':'S1','alignment':str(bam),'index':str(bam)+'.bai','seeds_bed':str(self.bed),'reference_build':'hg19'}
    def test_explicit_empty_seed_is_not_negative(self):
        r=self.row(self.bam());self.bed.write_text('')
        rows,_=validate(self.config,[r]);self.assertEqual(rows[0]['status'],'failed')
        self.config['allow_empty_seeds']=True
        rows,_=validate(self.config,[r]);self.assertEqual(rows[0]['status'],'no_seeds');self.assertEqual(rows[0]['seeds'],[])
    def test_full_suite_checks_normal_and_needs_no_seed(self):
        import shutil
        bam=self.bam();normal=self.root/'normal.bam';shutil.copy2(bam,normal);shutil.copy2(str(bam)+'.bai',str(normal)+'.bai')
        r=self.row(bam);r.pop('seeds_bed');r.update(normal_alignment=str(normal),normal_index=str(normal)+'.bai')
        rows,_=validate_full(self.config,[r]);self.assertEqual(rows[0]['status'],'ok')
        r['normal_index']=str(self.root/'missing.bai')
        rows,_=validate_full(self.config,[r]);self.assertEqual(rows[0]['status'],'failed');self.assertIn('Matched normal',rows[0]['reason'])
    def test_full_suite_rejects_same_tumor_normal_file(self):
        bam=self.bam();r=self.row(bam);r.pop('seeds_bed');r.update(normal_alignment=str(bam),normal_index=str(bam)+'.bai')
        rows,_=validate_full(self.config,[r]);self.assertEqual(rows[0]['status'],'failed');self.assertIn('same alignment',rows[0]['reason'])
    def test_good(self):
        rows,_=validate(self.config,[self.row(self.bam())]);self.assertEqual(rows[0]['status'],'ok')
    def test_build_mismatch(self):
        r=self.row(self.bam());r['reference_build']='GRCh38';rows,_=validate(self.config,[r]);self.assertEqual(rows[0]['status'],'failed')
    def test_missing_contig_strict_and_explicit_scope(self):
        r=self.row(self.bam(True));rows,_=validate(self.config,[r]);self.assertEqual(rows[0]['status'],'failed')
        self.config['reference']['allowed_missing_reference_contigs']=['virus'];rows,_=validate(self.config,[r]);self.assertEqual(rows[0]['missing_reference_contigs'],['virus'])
    def test_bad_coordinates(self):
        self.bed.write_text('chr1\t10\t1001\tx\n')
        with self.assertRaises(ValueError):read_seeds(self.bed,{'chr1':1000})
    def test_duplicate_seed(self):
        self.bed.write_text('chr1\t1\t20\tx\nchr1\t50\t60\tx\n')
        with self.assertRaises(ValueError):read_seeds(self.bed,{'chr1':1000})
    def test_cram_whole_file_conversion_is_disabled(self):
        bam=self.bam();cram=self.root/'test.cram'
        with pysam.AlignmentFile(str(bam)) as src, pysam.AlignmentFile(str(cram),'wc',header=src.header,reference_filename=str(self.fasta)) as dst:
            for r in src:dst.write(r)
        pysam.index(str(cram));r=self.row(cram);r['index']=str(cram)+'.crai'
        with self.assertRaisesRegex(ValueError,'conversion is disabled'):
            prepare_alignment(r,self.config,self.root/'stage')

if __name__=='__main__':unittest.main()
