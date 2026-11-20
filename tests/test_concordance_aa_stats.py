import json
from pathlib import Path
import tempfile
import unittest

from ascompare.common import digest
from ascompare.concordance_aa_stats import seed_stage_measurement, stage_measurement, validated_values
from ascompare.versions import load_lock


class NativeHandoffTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name).resolve()
        self.attempt = self.root/'upstream'
        private = self.attempt/'private_AA_DATA_REPO'
        private.mkdir(parents=True)
        self.ref = self.root/'resources'
        (self.ref/'hg19').mkdir(parents=True)
        (private/'hg19').symlink_to(self.ref/'hg19')
        self.bam = self.root/'synthetic file.bam'
        self.bam.write_bytes(b'synthetic identity only')
        self.stats = private/'coverage.stats'
        self.fields = ['42','41','4','43','42','5','150','320','15','275','365','3','0.95','3',
                       str(self.bam.stat().st_size)]
        self.stats.write_text('\t'.join([str(self.bam), *self.fields])+'\n')
        self.sample = dict(alignment=str(self.bam), sample_id='SYNTHETIC', reference_build='hg19')
        self.sample['seeds'] = [dict(seed_id='SYNTHETIC_SEED')]
        self.sample['_input_identities'] = [dict(path=str(self.bam), bytes=self.bam.stat().st_size,
                                                mtime_ns=self.bam.stat().st_mtime_ns)]
        (self.root/'src').mkdir()
        self.runtime = dict(id='suite-1.3.3', aa_python='python', aa_script=str(self.root/'src/AmpliconArchitect.py'))
        entry = next(e for e in load_lock()['versions'] if e['id']==self.runtime['id'])
        self.verified = dict(aa=dict(entry['aa'], source_root=str(self.root), python=str(Path('python').absolute())), ready=True)
        self.config = dict(reference={'aa_data_repo':str(self.ref)}, _io_policy_path='policy.json')
        self.prov = self.attempt/'upstream_provenance.json'
        helper = self.root/'src/amplified_intervals.py'
        helper.write_text('# synthetic helper')
        args = ['--ref','hg19','--bam',str(self.bam),'--no_cstats']
        self.prov.write_text(json.dumps(dict(version=self.runtime['id'],io_policy='policy.json',
            amplified_intervals=dict(path=str(helper), sha256=digest(helper), args=args))))
        self.rows = [dict(status='retained',version=self.runtime['id'],sample_id='SYNTHETIC',
                          seed_id='SYNTHETIC_SEED',provenance=str(self.prov))]
        self.log = self.attempt/'seed_processing.log'
        self.log.write_text(json.dumps({'argv':['python','guard','--policy','policy.json','--script',str(helper),'--',*args]})+'\n')

    def measurement(self):
        return seed_stage_measurement(self.attempt, self.rows, self.sample, self.runtime, self.verified, self.config)

    def test_exact_native_handoff_with_whitespace_and_content_binding(self):
        evidence = self.measurement()
        self.assertEqual(evidence['mode'], 'same_version_seed_stage')
        self.assertEqual(evidence['sha256'], digest(self.stats))
        job = self.root/'downstream'
        (job/'reference').mkdir(parents=True)
        alias = Path(stage_measurement(evidence, job, self.bam))
        self.assertEqual(alias.resolve(), self.bam)
        self.assertEqual((job/'reference'/'coverage.stats').read_text().split(), [str(alias),*self.fields])
        self.stats.write_text(self.stats.read_text().replace('\t42\t','\t41\t',1))
        with self.assertRaisesRegex(ValueError, 'changed'):
            stage_measurement(evidence, job, self.bam)

    def test_reject_different_runtime_or_incomplete_stage(self):
        self.rows[0]['status'] = 'failed'
        self.assertIsNone(self.measurement())
        self.rows[0]['status'] = 'retained'
        self.verified['aa']['commit'] = 'different'
        self.assertIsNone(self.measurement())

    def test_suite_defaults_accept_verified_native_stats_without_no_cstats(self):
        self.config['parameter_policy'] = 'suite_defaults'
        self.assertEqual(self.measurement()['mode'], 'native_measurement')
        upstream = json.loads(self.prov.read_text())
        upstream['amplified_intervals']['args'].remove('--no_cstats')
        self.prov.write_text(json.dumps(upstream))
        command = json.loads(self.log.read_text())
        command['argv'].remove('--no_cstats')
        self.log.write_text(json.dumps(command) + '\n')
        evidence = self.measurement()
        self.assertEqual(evidence['mode'], 'same_version_seed_stage')
        self.assertEqual(evidence['fields'], self.fields)
        job = self.root/'downstream'
        (job/'reference').mkdir(parents=True)
        self.assertEqual(Path(stage_measurement(evidence, job, self.bam, allow_alias=False)), self.bam)
        self.assertEqual((job/'reference'/'coverage.stats').read_bytes(), self.stats.read_bytes())
        self.assertFalse((job/'coverage_input.bam').exists())

    def test_command_and_measurement_integrity(self):
        self.log.write_text(self.log.read_text().replace('policy.json','other.json'))
        self.assertEqual(self.measurement()['mode'], 'native_measurement')
        for field, bad in [(0,'nan'), (14,'1'), (11,'1'), (13,'1')]:
            with self.subTest(field=field):
                values = self.fields.copy()
                values[field] = bad
                self.stats.write_text('\t'.join([str(self.bam), *values])+'\n')
                with self.assertRaises(ValueError):
                    validated_values(self.stats,self.bam)
