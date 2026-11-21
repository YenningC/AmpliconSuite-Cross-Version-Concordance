"""Optional native WGS controls for exact reuse and parallel execution."""

import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from ascompare import concordance_global_depth
from ascompare.concordance_global_cnv import global_cns
from ascompare.versions import load_lock


NATIVE_VERSIONS = (
    ('suite-1.5.2', '0.9.10', 'GLOBAL_CNV_PYTHON_0910'),
    ('suite-1.6.2', '0.9.13', 'GLOBAL_CNV_PYTHON_0913'),
)
PRODUCTS = (
    'normal.targetcoverage.cnn',
    'tumor.targetcoverage.cnn',
    'reference.cnn',
    'tumor.cnr',
    'tumor.cns',
)


class NativeGlobalTest(unittest.TestCase):
    def make_inputs(self, root):
        import pysam

        chromosomes = ('chr1', 'chr2', 'chr3')
        chromosome_size = 1_000_000
        fasta = root / 'synthetic.fa'
        with fasta.open('w') as stream:
            for chrom in chromosomes:
                stream.write(f'>{chrom}\n' + 'ACGT' * (chromosome_size // 4) + '\n')
        pysam.faidx(str(fasta))
        resources = root / 'hg19'
        resources.mkdir()
        (resources / 'hg19_noAlt.fa.fai').write_bytes(Path(str(fasta) + '.fai').read_bytes())
        reference = dict(build='hg19', fasta=str(fasta), aa_data_repo=str(root))
        header = dict(
            HD={'VN': '1.6', 'SO': 'coordinate'},
            SQ=[{'SN': chrom, 'LN': chromosome_size} for chrom in chromosomes],
        )
        sample = dict(sample_id='SYNTHETIC', reference_build='hg19')
        for stem in ('normal', 'tumor'):
            bam = root / (stem + '.bam')
            with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
                for chrom_index, chrom in enumerate(chromosomes):
                    gain_step = 12 if stem == 'tumor' and chrom == 'chr1' else 75
                    for first, last, step in ((0, 200_000, 75), (200_000, 500_000, gain_step),
                                              (500_000, 990_000, 75)):
                        for start in range(first, last, step):
                            read = pysam.AlignedSegment()
                            read.query_name = f'SYNTHETIC_{chrom}_{start}'
                            read.query_sequence = 'ACGT' * 25
                            read.flag = 0
                            read.reference_id = chrom_index
                            read.reference_start = start
                            read.mapping_quality = 60
                            read.cigar = [(0, 100)]
                            out.write(read)
            pysam.index(str(bam))
            prefix = 'normal_' if stem == 'normal' else ''
            sample[prefix + 'alignment'], sample[prefix + 'index'] = str(bam), str(bam) + '.bai'
        return sample, reference

    def run_native(self, command, env, cwd):
        result = subprocess.run(
            command, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, timeout=180,
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        return result.stdout

    def native_environment(self, python, extra, root):
        return dict(
            os.environ,
            PATH=os.pathsep.join([str(Path(python).parent), *extra, os.environ['PATH']]),
            NUMEXPR_MAX_THREADS='1',
            MPLCONFIGDIR=str(root / 'mplcache'),
        )

    def native_batch(self, python, script, extra, sample, fasta, out, processes):
        out.mkdir()
        env = self.native_environment(python, extra, out.parent)
        self.run_native(
            [python, script, 'batch', sample['alignment'], '-m', 'wgs', '--normal',
             sample['normal_alignment'], '--fasta', fasta, '-p', str(processes), '-d', out],
            env, out,
        )
        self.run_native(
            [python, script, 'segment', out / 'tumor.cnr', '-m', 'cbs', '-p', str(processes),
             '-o', out / 'tumor.cns'],
            env, out,
        )

    def assert_equal_products(self, expected, actual):
        for name in PRODUCTS:
            with self.subTest(product=name):
                self.assertEqual((expected / name).read_bytes(), (actual / name).read_bytes())
        with (actual / 'tumor.cns').open() as stream:
            self.assertTrue(any(
                row['chromosome'] == 'chr1' and float(row['log2']) > 1.2
                for row in csv.DictReader(stream, delimiter='\t')
            ), 'The synthetic amplification should survive segmentation')

    @unittest.skipUnless(
        all(os.getenv(envkey) for _, _, envkey in NATIVE_VERSIONS),
        'Provide both pinned historical CNVkit interpreter paths for the native WGS control',
    )
    def test_native_cross_version_cache_matches_unmodified_batch(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            sample, reference = self.make_inputs(root)
            extra = [p for p in os.getenv('CNVKIT_SMOKE_EXTRA_PATH', '').split(os.pathsep) if p]
            serial_cache = root / 'shared-serial-cache'
            lock = {entry['id']: entry for entry in load_lock()['versions']}
            previous_engine = None
            for suite, version, envkey in NATIVE_VERSIONS:
                with self.subTest(version=version):
                    python = os.environ[envkey]
                    script = str(Path(python).parent / 'cnvkit.py')
                    runtime = dict(
                        id=suite, cnvkit_version=version, cnvkit_python=python, cnvkit_script=script,
                        cnvkit_env={'MPLCONFIGDIR': str(root / 'mplcache')}, extra_path=extra,
                    )
                    result = global_cns(sample, runtime, lock[suite], dict(reference=reference), serial_cache)
                    serial = Path(result['cns']).parent
                    engine = json.loads(self.run_native(
                        [python, '-c', 'import json, runpy, sys; '
                         'print(json.dumps(runpy.run_path(sys.argv[1])["engine"]()))',
                         concordance_global_depth.__file__],
                        self.native_environment(python, extra, root), root,
                    ))
                    reuse = engine == previous_engine
                    serial_counts = json.loads(Path(result['provenance']).read_text())['bedcov_counts']
                    self.assertEqual(serial_counts, {
                        stem: {'hits': int(reuse), 'misses': int(not reuse)} for stem in ('normal', 'tumor')
                    })
                    previous_engine = engine

                    control = root / ('plain-' + version)
                    reduced = next((serial_cache / 'reference').glob('*/*_reduced.fa'))
                    self.native_batch(python, script, extra, sample, reduced, control, 1)
                    self.assert_equal_products(control, serial)

                    parallel = global_cns(
                        sample, runtime, lock[suite], dict(reference=reference, global_cnv={'processes': 4}),
                        root / ('cold-parallel-' + version),
                    )
                    parallel_dir = Path(parallel['cns']).parent
                    self.assert_equal_products(control, parallel_dir)
                    for stage in ('batch', 'segment'):
                        with (parallel_dir / f'global_cnv_{stage}.log').open() as stream:
                            argv = json.loads(stream.readline())['argv']
                        self.assertEqual(argv[argv.index('-p') + 1], '4')
                    parallel_counts = json.loads(Path(parallel['provenance']).read_text())['bedcov_counts']
                    self.assertEqual(set(parallel_counts), {'normal', 'tumor'})
                    for counts in parallel_counts.values():
                        self.assertEqual(counts['hits'], 0)
                        self.assertGreater(counts['misses'], 0)

    @unittest.skipUnless(os.getenv('GLOBAL_CNV_PYTHON_0914'), 'Provide a CNVkit 0.9.14 interpreter')
    def test_native_0914_parallel_batch_matches_serial(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            sample, reference = self.make_inputs(root)
            extra = [p for p in os.getenv('CNVKIT_SMOKE_EXTRA_PATH', '').split(os.pathsep) if p]
            python = os.environ['GLOBAL_CNV_PYTHON_0914']
            script = str(Path(python).parent / 'cnvkit.py')
            version = self.run_native([python, script, 'version'], self.native_environment(python, extra, root), root)
            self.assertEqual(version.strip().splitlines()[-1], '0.9.14')
            serial, parallel = root / 'native-serial', root / 'native-parallel'
            self.native_batch(python, script, extra, sample, reference['fasta'], serial, 1)
            self.native_batch(python, script, extra, sample, reference['fasta'], parallel, 4)
            self.assert_equal_products(serial, parallel)
