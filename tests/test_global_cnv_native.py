"""Opt-in positive control: both native CNVkit versions, with and without reuse."""

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from ascompare.concordance_global_cnv import global_cns
from ascompare.versions import load_lock


@unittest.skipUnless(
    all(os.getenv(k) for k in ('GLOBAL_CNV_PYTHON_0910', 'GLOBAL_CNV_PYTHON_0913')),
    'Provide both historical CNVkit interpreter paths for the native WGS control',
)
class NativeGlobalTest(unittest.TestCase):
    def test_native_cross_version_cache_matches_unmodified_batch(self):
        import pysam

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            fasta = root / 'synthetic.fa'
            fasta.write_text('>chr1\n' + 'ACGT' * 250000 + '\n')
            pysam.faidx(str(fasta))
            resources = root / 'hg19'
            resources.mkdir()
            (resources / 'hg19_noAlt.fa.fai').write_bytes(Path(str(fasta) + '.fai').read_bytes())
            cfg = dict(reference=dict(build='hg19', fasta=str(fasta), aa_data_repo=str(root)))
            sample = dict(sample_id='SYNTHETIC', reference_build='hg19')
            for stem in ('normal', 'tumor'):
                bam = root / (stem + '.bam')
                header = dict(HD={'VN': '1.6', 'SO': 'coordinate'}, SQ=[{'SN': 'chr1', 'LN': 1000000}])
                with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
                    for start in range(0, 990000, 75):
                        read = pysam.AlignedSegment()
                        read.query_name = f'SYNTHETIC_{start}'
                        read.query_sequence = 'ACGT' * 25
                        read.flag = 0
                        read.reference_id = 0
                        read.reference_start = start
                        read.mapping_quality = 60
                        read.cigar = [(0, 100)]
                        out.write(read)
                pysam.index(str(bam))
                prefix = 'normal_' if stem == 'normal' else ''
                sample[prefix + 'alignment'], sample[prefix + 'index'] = str(bam), str(bam) + '.bai'
            versions = [
                ('suite-1.5.2', '0.9.10', 'GLOBAL_CNV_PYTHON_0910'),
                ('suite-1.6.2', '0.9.13', 'GLOBAL_CNV_PYTHON_0913'),
            ]
            for suite, version, envkey in versions:
                python = os.environ[envkey]
                script = str(Path(python).parent / 'cnvkit.py')
                extra = os.getenv('CNVKIT_SMOKE_EXTRA_PATH', '').split(os.pathsep)
                runtime = dict(
                    id=suite,
                    cnvkit_version=version,
                    cnvkit_python=python,
                    cnvkit_script=script,
                    extra_path=[p for p in extra if p],
                )
                entry = next(e for e in load_lock()['versions'] if e['id'] == suite)
                result = global_cns(sample, runtime, entry, cfg, root / 'global_cnv_cache')
                producer = Path(result['cns']).parent
                control = root / ('plain-' + version)
                control.mkdir()
                env = dict(
                    os.environ,
                    PATH=os.pathsep.join([str(Path(python).parent), *extra, os.environ['PATH']]),
                    NUMEXPR_MAX_THREADS='1',
                    MPLCONFIGDIR=str(root / 'mplcache'),
                )
                reduced = next((root / 'global_cnv_cache/reference').glob('*/*_reduced.fa'))
                for args in (
                    [
                        'batch',
                        sample['alignment'],
                        '-m',
                        'wgs',
                        '--normal',
                        sample['normal_alignment'],
                        '--fasta',
                        reduced,
                        '-p',
                        '1',
                        '-d',
                        control,
                    ],
                    ['segment', control / 'tumor.cnr', '-m', 'cbs', '-p', '1', '-o', control / 'tumor.cns'],
                ):
                    subprocess.run(
                        [python, script, *args],
                        check=True,
                        cwd=control,
                        env=env,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                    )
                for name in (
                    'normal.targetcoverage.cnn',
                    'tumor.targetcoverage.cnn',
                    'reference.cnn',
                    'tumor.cnr',
                    'tumor.cns',
                ):
                    self.assertEqual(
                        (producer / name).read_bytes(), (control / name).read_bytes(), version + ' ' + name
                    )
                # Instrumentation is local; different bedcov builds may legitimately miss.
                self.assertEqual(
                    set(json.loads(Path(result['provenance']).read_text())['bedcov_counts']),
                    {'tumor', 'normal'},
                )
