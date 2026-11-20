"""Component cache integration on guarded synthetic subprocesses, never patient data.

These callers are deliberately small test programs. Native official AA/AC A/B
checks are a separate recorded experiment, not a claim made by this fixture.
"""
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import test_concordance_bounded_workflow as bounded_fixture
from ascompare.common import atomic_json, digest
from ascompare.concordance import load_manifest, run_bounded_downstream
from ascompare.concordance_demo import write_demo_manifest
from test_parsing import add_synthetic_amplicon


class ComponentCacheIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.fixture = bounded_fixture.BoundedConcordanceWorkflowTests(
            'test_actual_guarded_call_resume_and_cache_labels_are_reparsed_from_raw')
        self.fixture.setUp()
        self.root = self.fixture.root
        self.config = self.fixture.config
        self.cache = self.root / '.component_cache'
        self.config['_component_cache_dir'] = str(self.cache)
        self.config['aa_parameters']['runmode'] = 'CYCLES'
        self.fixture.aa.write_text(self.fixture.aa.read_text() +
            "\nprint('SYNTHETIC runmode', a[a.index('--runmode')+1])\n")

    def tearDown(self):
        self.fixture.tearDown()

    def run_job(self, version='suite-1.3.3', job=None, verified=None, runtime_updates=None):
        runtime = self.fixture.make_runtime(version)
        runtime.update(runtime_updates or {})
        entry = self.fixture.entries[version]
        verification = verified or self.fixture.make_verified(runtime, entry)
        return run_bounded_downstream(self.fixture.sample, self.fixture.seeds,
            runtime, entry, verification, self.config, self.root / (job or version))

    def assert_ok(self, record):
        self.assertEqual(record['status'], 'ok', record.get('reason'))
        self.assertTrue(record['comparison_eligible'])
        self.assertEqual(record['labels'], ['ecDNA'])

    def test_same_aa_in_two_suite_recipes_runs_once_and_distinct_ac_twice(self):
        first = self.run_job()
        second = self.run_job('suite-1.3.4')
        self.assert_ok(first)
        self.assert_ok(second)
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AC'])
        self.assertFalse(first['component_cache']['aa']['hit'])
        self.assertTrue(second['component_cache']['aa']['hit'])
        self.assertFalse(second['component_cache']['ac']['hit'])
        self.assertEqual(first['component_cache']['aa']['key'],
                         second['component_cache']['aa']['key'])
        self.assertNotEqual(first['component_cache']['ac']['key'],
                            second['component_cache']['ac']['key'])
        self.assertEqual(second['component_versions'], {'aa': '1.4.r2', 'ac': '1.2.1'})
        for record in (first, second):
            stages = json.loads((Path(record['run_dir']) / 'stages.json').read_text())
            self.assertEqual(stages['aa']['status'], 'complete')
            self.assertEqual(stages['ac']['status'], 'complete')
            self.assertTrue((Path(record['run_dir']) / 'ac.input').is_file())

    def test_suite_defaults_leave_native_hash_and_classification_unchanged(self):
        self.config['parameter_policy'] = 'suite_defaults'
        self.config.pop('aa_parameters')
        self.config.pop('ac_parameters')
        self.fixture.aa.write_text(self.fixture.aa.read_text() +
            "\nassert 'PYTHONHASHSEED' not in os.environ\n")
        with patch.dict(os.environ):
            os.environ.pop('PYTHONHASHSEED', None)
            first = self.run_job()
            second = self.run_job('suite-1.3.4')
            resumed = self.run_job()
        for record in (first, second):
            self.assert_ok(record)
            attempt = Path(record['run_dir'])
            stages = json.loads((attempt / 'stages.json').read_text())
            aa = stages['aa']['argv']
            self.assertEqual(aa[aa.index('--runmode') + 1], 'FULL')
            self.assertEqual(aa[aa.index('--downsample') + 1], '10')
            self.assertNotIn('--no_cstats', aa)
            self.assertNotIn('--min_size', stages['ac']['argv'])
            self.assertNotIn('--min_flow', stages['ac']['argv'])
            self.assertIsNone(json.loads((attempt / 'aa_background.json').read_text())['python_hash_seed'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AA', 'AC'])
        self.assertNotEqual(first['component_cache']['aa']['key'], second['component_cache']['aa']['key'])
        self.assertEqual(first['run_dir'], resumed['run_dir'])

    def test_suite_defaults_accept_native_background_recomputation(self):
        self.config['parameter_policy'] = 'suite_defaults'
        bam = self.fixture.bam
        stats = self.root / 'seed-stage.stats'
        # The 11 in the native field is computed coverage support, not a
        # requested minimum to substitute in the AA command.
        fields = ['20','20','3','20','20','5','150','330','25','255','405','11','0.95','3',str(bam.stat().st_size)]
        stats.write_text('\t'.join([str(bam), *fields])+'\n')
        self.fixture.aa.write_text(self.fixture.aa.read_text() + """
stats=pathlib.Path(os.environ['AA_DATA_REPO'])/'coverage.stats'
with stats.open('a') as handle:
    handle.write(stats.read_text())
print('cstats not found, generating coverage statistics')
""")
        runtime = self.fixture.make_runtime('suite-1.3.7')
        entry = self.fixture.entries['suite-1.3.7']
        verified = self.fixture.make_verified(runtime, entry)
        measurement = dict(mode='same_version_seed_stage', path=str(stats), fields=fields, sha256=digest(stats))
        result = run_bounded_downstream(self.fixture.sample, self.fixture.seeds, runtime, entry,
                                        verified, self.config, self.root/'native', measurement)
        self.assert_ok(result)
        record = json.loads((Path(result['run_dir'])/'aa_background.json').read_text())
        self.assertFalse(record['reuse_confirmed'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC'])

    def test_ac_parameter_change_reuses_aa_and_reexecutes_ac(self):
        first = self.run_job()
        self.config['ac_parameters']['min_size'] = 25
        second = self.run_job()
        self.assert_ok(second)
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AC'])
        self.assertTrue(second['component_cache']['aa']['hit'])
        self.assertFalse(second['component_cache']['ac']['hit'])
        self.assertNotEqual(first['run_dir'], second['run_dir'])

    def test_aa_parameter_change_invalidates_aa_and_ac(self):
        self.assert_ok(self.run_job())
        self.config['aa_parameters']['downsample'] = 10
        result = self.run_job()
        self.assert_ok(result)
        self.assertFalse(result['component_cache']['aa']['hit'])
        self.assertFalse(result['component_cache']['ac']['hit'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AA', 'AC'])

    def test_full_and_cycles_have_distinct_aa_cache_keys_and_recorded_argv(self):
        cycles = self.run_job()
        self.config['aa_parameters']['runmode'] = 'FULL'
        full = self.run_job()
        self.assert_ok(full)
        self.assertNotEqual(cycles['component_cache']['aa']['key'], full['component_cache']['aa']['key'])
        for mode, record in [('CYCLES', cycles), ('FULL', full)]:
            stage = json.loads((Path(record['run_dir']) / 'stages.json').read_text())['aa']
            self.assertEqual(stage['argv'][stage['argv'].index('--runmode') + 1], mode)
            self.assertIn('SYNTHETIC runmode ' + mode, (Path(record['run_dir']) / 'aa.log').read_text())
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AA', 'AC'])

    def test_verified_identity_change_cannot_reuse_aa(self):
        self.assert_ok(self.run_job())
        runtime = self.fixture.make_runtime('suite-1.3.3')
        verified = self.fixture.make_verified(runtime, self.fixture.entries['suite-1.3.3'])
        verified['aa']['manifest_sha256'] = 'SYNTHETIC_DIFFERENT_SOURCE'
        result = self.run_job(verified=verified)
        self.assert_ok(result)
        self.assertFalse(result['component_cache']['aa']['hit'])
        self.assertEqual(self.fixture.trace().count('AA'), 2)

    def test_actual_aa_script_change_invalidates_even_with_same_declared_version(self):
        self.assert_ok(self.run_job())
        self.fixture.aa.write_text(self.fixture.aa.read_text() + '\n# SYNTHETIC source changed\n')
        result = self.run_job()
        self.assert_ok(result)
        self.assertFalse(result['component_cache']['aa']['hit'])
        self.assertEqual(self.fixture.trace().count('AA'), 2)

    def test_ac_script_change_preserves_aa_cache(self):
        self.assert_ok(self.run_job())
        self.fixture.ac.write_text(self.fixture.ac.read_text() + '\n# SYNTHETIC classifier changed\n')
        result = self.run_job()
        self.assert_ok(result)
        self.assertTrue(result['component_cache']['aa']['hit'])
        self.assertFalse(result['component_cache']['ac']['hit'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AC'])

    def test_changed_seed_coordinates_cannot_borrow_previous_aa(self):
        self.assert_ok(self.run_job())
        self.fixture.seeds[0] = dict(self.fixture.seeds[0], start=111)
        result = self.run_job()
        self.assert_ok(result)
        self.assertFalse(result['component_cache']['aa']['hit'])
        self.assertEqual(self.fixture.trace().count('AA'), 2)

    def test_effective_hash_seed_is_part_of_aa_identity(self):
        self.assert_ok(self.run_job())
        env = dict(self.fixture.make_runtime('suite-1.3.3')['env'], PYTHONHASHSEED='1234')
        result = self.run_job(runtime_updates={'env': env})
        self.assert_ok(result)
        self.assertFalse(result['component_cache']['aa']['hit'])
        self.assertEqual(self.fixture.trace().count('AA'), 2)

    def test_policy_file_content_change_invalidates_aa(self):
        self.assert_ok(self.run_job())
        policy = json.loads(self.fixture.policy.read_text())
        policy['max_total_reads'] += 1
        atomic_json(self.fixture.policy, policy)
        result = self.run_job()
        self.assert_ok(result)
        self.assertFalse(result['component_cache']['aa']['hit'])
        self.assertEqual(self.fixture.trace().count('AA'), 2)

    def test_unverified_runtime_cannot_borrow_successful_component(self):
        self.assert_ok(self.run_job())
        runtime = self.fixture.make_runtime('suite-1.3.4')
        verified = self.fixture.make_verified(runtime, self.fixture.entries['suite-1.3.4'])
        verified.update(ready=False, diagnostics=['SYNTHETIC verification failure'])
        result = self.run_job('suite-1.3.4', verified=verified)
        self.assertEqual(result['status'], 'failed')
        self.assertFalse(result['comparison_eligible'])
        self.assertIn('verification failure', result['reason'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC'])

    def test_failed_ac_keeps_successful_aa_for_retry(self):
        marker = self.root / 'SYNTHETIC_ac_failure'
        marker.touch()
        self.fixture.ac.write_text(self.fixture.ac.read_text() +
            '\nif pathlib.Path(' + repr(str(marker)) + ').exists(): sys.exit(7)\n')
        failed = self.run_job()
        self.assertEqual((failed['status'], failed['stage']), ('failed', 'ac'))
        self.assertFalse(failed['comparison_eligible'])
        marker.unlink()
        recovered = self.run_job()
        self.assert_ok(recovered)
        self.assertTrue(recovered['component_cache']['aa']['hit'])
        self.assertFalse(recovered['component_cache']['ac']['hit'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AC'])

    def test_failed_aa_is_never_published_as_success(self):
        marker = self.root / 'SYNTHETIC_aa_failure'
        marker.touch()
        self.fixture.aa.write_text(self.fixture.aa.read_text() +
            '\nif pathlib.Path(' + repr(str(marker)) + ').exists(): sys.exit(9)\n')
        failed = self.run_job()
        self.assertEqual((failed['status'], failed['stage']), ('failed', 'aa'))
        self.assertFalse(failed['comparison_eligible'])
        marker.unlink()
        recovered = self.run_job()
        self.assert_ok(recovered)
        self.assertFalse(recovered['component_cache']['aa']['hit'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AA', 'AC'])

    def test_corrupt_shared_aa_payload_forces_native_aa_retry(self):
        self.assert_ok(self.run_job())
        graph = next((self.cache / 'aa').rglob('*_graph.txt'))
        graph.write_text('SYNTHETIC truncated graph\n')
        result = self.run_job('suite-1.3.4')
        self.assert_ok(result)
        self.assertFalse(result['component_cache']['aa']['hit'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AA', 'AC'])

    def test_malformed_native_aa_output_is_not_a_reusable_success(self):
        marker = self.root / 'SYNTHETIC_malformed_aa'
        marker.touch()
        self.fixture.aa.write_text(self.fixture.aa.read_text() + '\nif pathlib.Path(' + repr(str(marker)) +
            ").exists():\n (prefix.parent/'reconstruction_amplicon1_graph.txt').write_text('SYNTHETIC malformed graph\\n')\n")
        failed = self.run_job()
        self.assertNotIn(failed['status'], ('ok', 'not_detected'))
        self.assertFalse(failed['comparison_eligible'])
        marker.unlink()
        recovered = self.run_job()
        self.assert_ok(recovered)
        self.assertFalse(recovered['component_cache']['aa']['hit'])
        self.assertEqual(self.fixture.trace().count('AA'), 2)

    def test_bad_amplicon_preserves_other_valid_amplicon_and_is_not_cached(self):
        add_synthetic_amplicon(self.fixture.template, kind='Cyclic')
        bad = self.fixture.template / 'aa/reconstruction_amplicon2_graph.txt'
        bad.write_text('SYNTHETIC malformed second amplicon\n')
        feature = self.fixture.template / 'ac/classification_classification_bed_files/reconstruction_amplicon2_ecDNA_1_intervals.bed'
        feature.write_text('chr1\t100\t200\nchr2\t1000\t1100\n')
        self.fixture.ac.write_text(self.fixture.ac.read_text().replace('assert len(rows)==1', 'assert len(rows)==2'))
        result = self.run_job()
        self.assertEqual(result['status'], 'invalid')
        self.assertEqual(result['coverage_status'], 'partial')
        self.assertEqual(result['evaluable_amplicon_count'], 1)
        self.assertEqual([(a['amplicon_id'], a['status']) for a in result['amplicons']],
                         [('amplicon1', 'ok'), ('amplicon2', 'invalid')])
        self.assertEqual(result['features'][0]['label'], 'ecDNA')
        self.assertFalse(result['comparison_eligible'])
        self.assertFalse(list((self.cache / 'aa').glob('*/manifest.json')))
        self.assertFalse(list((self.cache / 'ac').glob('*/manifest.json')))

    def test_corrupt_ac_payload_reuses_aa_but_reexecutes_ac(self):
        self.assert_ok(self.run_job())
        bed = next((self.cache / 'ac').rglob('*intervals.bed'))
        bed.write_text('SYNTHETIC invalid BED\n')
        result = self.run_job(job='same_components_new_job')
        self.assert_ok(result)
        self.assertTrue(result['component_cache']['aa']['hit'])
        self.assertFalse(result['component_cache']['ac']['hit'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AC'])

    def test_unreadable_component_manifest_is_a_miss_and_retries_aa(self):
        self.assert_ok(self.run_job())
        manifest = next((self.cache / 'aa').rglob('manifest.json'))
        manifest.write_text('{SYNTHETIC truncated JSON')
        recovered = self.run_job('suite-1.3.4')
        self.assert_ok(recovered)
        self.assertFalse(recovered['component_cache']['aa']['hit'])
        self.assertEqual(self.fixture.trace().count('AA'), 2)

    def test_full_job_resume_reports_cache_hits_without_reexecuting_components(self):
        first = self.run_job()
        resumed = self.run_job()
        self.assert_ok(resumed)
        self.assertEqual(first['run_dir'], resumed['run_dir'])
        self.assertTrue(resumed['downstream_cache_hit'])
        self.assertTrue(resumed['component_cache']['aa']['hit'])
        self.assertTrue(resumed['component_cache']['ac']['hit'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC'])

    def test_zero_amplicons_can_reuse_aa_and_never_claim_ac_execution(self):
        self.fixture.aa.write_text(self.fixture.aa.read_text() +
            "\nfor p in prefix.parent.glob('*_amplicon*'): p.unlink()\n" +
            "pathlib.Path(str(prefix)+'_summary.txt').write_text('#Amplicons = 0\\n')\n")
        first = self.run_job()
        second = self.run_job('suite-1.3.4')
        for record in (first, second):
            self.assertEqual(record['status'], 'not_detected')
            self.assertFalse(record['detected'])
            stages = json.loads((Path(record['run_dir']) / 'stages.json').read_text())
            self.assertEqual(stages['ac']['status'], 'skipped_zero_amplicons')
        self.assertTrue(second['component_cache']['aa']['hit'])
        self.assertNotIn('ac', second['component_cache'])
        self.assertEqual(self.fixture.trace(), ['AA'])

    def test_zero_amplicons_with_orphan_graph_is_not_negative_or_cached(self):
        self.fixture.aa.write_text(self.fixture.aa.read_text() +
            "\npathlib.Path(str(prefix)+'_summary.txt').write_text('#Amplicons = 0\\n')\n")
        result = self.run_job()
        self.assertNotIn(result['status'], ('ok', 'not_detected'))
        self.assertFalse(result['comparison_eligible'])
        self.assertFalse(list((self.cache / 'aa').glob('*/manifest.json')))

    def test_changed_bam_after_input_validation_cannot_reuse_cached_aa(self):
        self.assert_ok(self.run_job())
        self.fixture.bam.write_bytes(self.fixture.bam.read_bytes() + b'SYNTHETIC changed BAM')
        result = self.run_job()
        self.assertEqual(result['status'], 'failed')
        self.assertFalse(result['comparison_eligible'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC'])

    def test_changed_index_after_input_validation_cannot_reuse_cached_aa(self):
        self.assert_ok(self.run_job())
        index = Path(str(self.fixture.bam) + '.bai')
        index.write_bytes(index.read_bytes() + b'SYNTHETIC changed index')
        result = self.run_job()
        self.assertEqual(result['status'], 'failed')
        self.assertFalse(result['comparison_eligible'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC'])

    def test_cache_manifest_never_persists_arbitrary_environment_secrets(self):
        dummy_secret = 'SYNTHETIC_FAKE_TOKEN_78266432'
        env = dict(self.fixture.make_runtime('suite-1.3.3')['env'], PRIVATE_SERVICE_TOKEN=dummy_secret)
        self.assert_ok(self.run_job(runtime_updates={'env': env}))
        manifests = list(self.cache.rglob('manifest.json'))
        self.assertTrue(manifests)
        self.assertTrue(all(dummy_secret not in p.read_text() for p in manifests))

    def test_missing_feature_after_exit_zero_stays_missing_and_is_not_ac_hit(self):
        marker = self.root / 'SYNTHETIC_missing_feature'
        marker.touch()
        self.fixture.ac.write_text(self.fixture.ac.read_text() + '\nif pathlib.Path(' + repr(str(marker)) +
            ").exists():\n for f in prefix.parent.rglob('*intervals.bed'): f.unlink()\n")
        failed = self.run_job()
        self.assertEqual(failed['status'], 'missing')
        self.assertFalse(failed['comparison_eligible'])
        marker.unlink()
        recovered = self.run_job()
        self.assert_ok(recovered)
        self.assertTrue(recovered['component_cache']['aa']['hit'])
        self.assertFalse(recovered['component_cache']['ac']['hit'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AC'])

    def test_explicit_external_cache_reuse_leaves_producer_unchanged(self):
        producer = self.root / 'producer'
        self.config['_component_cache_dir'] = str(producer / '.component_cache')
        first = self.run_job()
        self.assert_ok(first)
        before = {str(p): digest(p) for p in producer.rglob('*') if p.is_file()}
        self.config['_component_cache_dir'] = str(self.root / 'consumer/.component_cache')
        self.config['reuse_aa_from'] = str(producer)
        result = self.run_job('suite-1.3.4')
        self.assert_ok(result)
        self.assertTrue(result['component_cache']['aa']['hit'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AC'])
        self.assertEqual(before, {str(p): digest(p) for p in producer.rglob('*') if p.is_file()})

    def test_missing_aa_package_probe_prevents_job_and_recipe_reuse(self):
        previous = None
        for version in ('suite-1.3.3', 'suite-1.3.3', 'suite-1.3.4'):
            runtime = self.fixture.make_runtime(version)
            verified = self.fixture.make_verified(runtime, self.fixture.entries[version])
            verified['aa'].pop('packages')
            result = self.run_job(version, verified=verified)
            self.assert_ok(result)
            self.assertFalse(result['component_cache']['aa']['hit'])
            self.assertFalse(result['component_cache']['ac']['hit'])
            if previous:
                self.assertNotEqual(previous, result['component_cache']['aa']['key'])
            previous = result['component_cache']['aa']['key']
        self.assertEqual(self.fixture.trace(), ['AA', 'AC'] * 3)

    def test_failed_ac_package_probe_reuses_only_verified_aa(self):
        runtime = self.fixture.make_runtime('suite-1.3.3')
        verified = self.fixture.make_verified(runtime, self.fixture.entries['suite-1.3.3'])
        verified['ac']['packages'] = {'returncode': 1, 'stdout': '', 'stderr': 'SYNTHETIC inventory failure'}
        first = self.run_job(verified=verified)
        second = self.run_job(verified=verified)
        self.assert_ok(first)
        self.assert_ok(second)
        self.assertTrue(second['component_cache']['aa']['hit'])
        self.assertFalse(second['component_cache']['ac']['hit'])
        self.assertNotEqual(first['component_cache']['ac']['key'], second['component_cache']['ac']['key'])
        self.assertEqual(self.fixture.trace(), ['AA', 'AC', 'AC'])


class ComponentPackageIdentityTests(unittest.TestCase):
    def test_successful_inventory_retains_its_existing_stable_identity(self):
        import hashlib
        from ascompare.concordance_component_identity import _packages
        info = {'packages': {'returncode': 0, 'stdout': 'numpy==1.0\nscipy==2.0\n'}}
        expected = dict(returncode=0, count=2,
            freeze_sha256=hashlib.sha256(json.dumps(['numpy==1.0', 'scipy==2.0']).encode()).hexdigest())
        self.assertEqual(_packages(info, 'aa'), expected)
        self.assertEqual(_packages(info, 'aa'), _packages(info, 'aa'))

    def test_unknown_or_malformed_inventory_always_gets_a_fresh_identity(self):
        from ascompare.concordance_component_identity import _packages
        for probe in (None, [], {'returncode': 3, 'stdout': 'numpy==1.0'},
                      {'returncode': 0}, {'returncode': 0, 'stdout': None}):
            with self.subTest(probe=probe):
                first = _packages({'packages': probe}, 'aa')
                second = _packages({'packages': probe}, 'aa')
                self.assertEqual(first['status'], 'unverified')
                self.assertNotEqual(first, second)


class ComponentModeManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.manifest = write_demo_manifest(Path(self.temp.name) / 'SYNTHETIC_inputs')

    def tearDown(self):
        self.temp.cleanup()

    def test_public_loader_accepts_only_full_or_cycles_reconstruction(self):
        data = json.loads(self.manifest.read_text())
        for mode in ('FULL', 'CYCLES'):
            with self.subTest(mode=mode):
                data['aa_parameters'] = {'runmode': mode}
                atomic_json(self.manifest, data)
                config, _, _ = load_manifest(self.manifest)
                self.assertEqual(config['aa_parameters']['runmode'], mode)

    def test_public_loader_rejects_modes_without_required_reconstruction_outputs(self):
        data = json.loads(self.manifest.read_text())
        for mode in ('BPGRAPH', 'SVVIEW', 'cycles', 'SYNTHETIC_UNKNOWN'):
            with self.subTest(mode=mode):
                data['aa_parameters'] = {'runmode': mode}
                atomic_json(self.manifest, data)
                with self.assertRaisesRegex(ValueError, 'runmode'):
                    load_manifest(self.manifest)

    def test_explicit_external_reuse_requires_an_existing_cache_output(self):
        data = json.loads(self.manifest.read_text())
        data['reuse_aa_from'] = str(self.manifest.parent / 'SYNTHETIC_absent_output')
        atomic_json(self.manifest, data)
        with self.assertRaisesRegex(ValueError, 'reuse_aa_from'):
            load_manifest(self.manifest)

    def test_suite_defaults_require_removing_algorithm_overrides(self):
        data = json.loads(self.manifest.read_text())
        data['parameter_policy'] = 'suite_defaults'
        overrides = (
            ('aa_parameters', {'runmode': 'CYCLES'}),
            ('ac_parameters', {'min_size': 50000}),
            ('local_cnv', {'cngain': 6}),
        )
        for key, value in overrides:
            with self.subTest(key=key):
                atomic_json(self.manifest, {**data, key: value})
                with self.assertRaisesRegex(ValueError, 'suite_defaults'):
                    load_manifest(self.manifest)
        atomic_json(self.manifest, data)
        self.assertEqual(load_manifest(self.manifest)[0]['parameter_policy'], 'suite_defaults')


if __name__ == '__main__':
    unittest.main()
