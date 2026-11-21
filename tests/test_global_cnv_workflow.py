"""Synthetic public workflow and source reproduction regressions."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pysam

from ascompare.common import atomic_json, digest
from ascompare.concordance import load_manifest, run_concordance
from ascompare.concordance_global_workflow import _cnvkit_runtime, _source_evidence


class GlobalWorkflowTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        self.genome = self.root / 'refs' / 'hg19'
        self.genome.mkdir(parents=True)
        self.fasta = self.genome / 'toy.fa'
        self.fasta.write_text('>chr1\n' + 'A' * 1000 + '\n')
        pysam.faidx(str(self.fasta))
        (self.genome / 'hg19_noAlt.fa.fai').write_text('chr1\t1000\n')
        (self.genome / 'hg19_cnvkit_filtered_ref.cnn').touch()
        self.bam = self.root / 'toy.bam'
        with pysam.AlignmentFile(str(self.bam), 'wb', header={
                'HD':{'SO':'coordinate'}, 'SQ':[{'SN':'chr1','LN':1000}]}) as stream:
            pass
        pysam.index(str(self.bam))
        self.source = 'suite-1.6.1'
        self.target = 'suite-1.6.2'
        self.archive = self.root / 'source'
        self.archive.mkdir()
        self.config = dict(mode='with_bam', synthetic=True, target_versions=[self.target],
                           reference=dict(build='hg19', fasta=str(self.fasta), aa_data_repo=str(self.genome.parent)),
                           samples=[dict(sample_id='SYNTHETIC', source_version=self.source, reference_build='hg19',
                                         alignment=str(self.bam), index=str(self.bam)+'.bai',
                                         source_result={'result_dir':str(self.archive)})],
                           sampling={'max_seeds':'all'}, statistics={'primary_pairs':[]},
                           runtimes=[{'id':v} for v in (self.source, self.target)])
        self.manifest = self.root / 'manifest.json'

    def load(self):
        atomic_json(self.manifest,self.config)
        return load_manifest(self.manifest)

    def test_unseeded_input_still_discovers_and_routes_only_native_target_seeds(self):
        self.config['cnvkit_runtimes'] = [
            dict(version='0.9.13', python='/exact/bin/python', script='/exact/bin/cnvkit.py',
                 extra_path=['/tools/R'])
        ]
        self.config['runtimes'][0].update(
            cnvkit_python='/exact/bin/python', cnvkit_script='/exact/bin/cnvkit.py',
            env={'DYLD_LIBRARY_PATH': '/aa/lib'},
        )
        config, lock, entries = self.load()
        self.assertEqual(config['cnv_scope'],'global')
        self.assertEqual(config['parameter_policy'],'suite_defaults')
        seeds = [dict(chrom='chr1',start=100,end=225,seed_id='native_1',source_fields=['8.5'])]
        observed = []
        def native(sample,runtime,entry,config,out):
            return dict(cnvkit_version='0.9.13',cns='toy.cns',reuse='same_recipe_cache',
                        evidence_level='historical_metadata_only',provenance='toy-metadata')
        def reconstruct(sample,selected,runtime,entry,verified,config,out,measurement):
            observed.append((runtime['id'],selected))
            return dict(sample_id=sample['sample_id'],version=runtime['id'],status='not_detected',
                        stage='aa',detected=False,labels=[],amplicons=[],features=[],intervals=[],
                        supported_labels=[],comparison_eligible=True)
        out = self.root / 'run'
        # Only routing is exercised here. Native scientific outputs have separate controls.
        with (patch('ascompare.concordance_global_cnv.global_cns',side_effect=native) as cnv,
              patch('ascompare.concordance_global_seeds.global_suite_seeds',side_effect=[[],seeds]) as seed_filter,
              patch('ascompare.concordance_official_solver.official_solver',return_value=('clarabel','synthetic')),
              patch('ascompare.concordance.verify_runtime',return_value={'ready':True}),
              patch('ascompare.concordance_global_seeds.global_seed_measurement',return_value=None),
              patch('ascompare.concordance.run_bounded_downstream',side_effect=reconstruct),
              patch('ascompare.concordance_report._plots',return_value=[]),
              patch('ascompare.concordance_cohort_plot.plot_cohort_compatibility')):
            run_concordance(config,lock,entries,out)
        self.assertEqual(cnv.call_count,2)
        for call in cnv.call_args_list:
            runtime = call.args[1]
            self.assertEqual(runtime['cnvkit_python'], '/exact/bin/python')
            self.assertEqual(runtime['cnvkit_script'], '/exact/bin/cnvkit.py')
            self.assertEqual(runtime['cnvkit_env'], {})
            self.assertEqual(runtime['cnvkit_extra_path'], ['/tools/R'])
        source_seed_runtime = seed_filter.call_args_list[0].args[1]
        self.assertEqual(source_seed_runtime['env'], {'DYLD_LIBRARY_PATH': '/aa/lib'})
        self.assertNotIn('cnvkit_env', source_seed_runtime)
        self.assertEqual(observed,[(self.target,seeds)])
        records=json.loads((out/'upstream_records.json').read_text())
        self.assertEqual({r['version']:r['status'] for r in records},
                         {self.source:'not_retained',self.target:'retained'})
        sets=json.loads((out/'global_upstream/SYNTHETIC/discovery.json').read_text())
        self.assertTrue(sets[self.target]['historical_provenance_limited'])
        self.assertTrue((out/'result.tsv').is_file())

    def test_global_cnv_registry_preserves_explicit_overrides_and_alternate_installations(self):
        tool = dict(python='/registry/python', script='/registry/cnvkit.py',
                    env={'R_LIBS': '/registry/R'}, extra_path=['/registry/bin'])
        registry = {'0.9.13': tool}
        runtime = dict(cnvkit_script=tool['script'], cnvkit_env={}, cnvkit_extra_path=[])
        updated = _cnvkit_runtime(runtime, '0.9.13', registry)
        self.assertEqual(updated['cnvkit_python'], tool['python'])
        self.assertEqual(updated['cnvkit_env'], {})
        self.assertEqual(updated['cnvkit_extra_path'], [])
        for alternate in (dict(cnvkit_python='/alternate/python'), dict(cnvkit_version='0.9.14')):
            self.assertIs(_cnvkit_runtime(alternate, '0.9.13', registry), alternate)
        self.assertIs(_cnvkit_runtime(runtime, '0.9.10', registry), runtime)

    def test_aa_recovery_keeps_verified_jobs_and_retries_the_unique_limit_failure(self):
        from ascompare.concordance_global_seeds import FINAL_BED, _native_seeds

        self.config['io_policy'] = {'max_unique_bp': 100, 'max_total_fetch_bp': 900}
        config, lock, entries = self.load()
        sid = 'SYNTHETIC'
        prior, resumed = self.root / 'original_run', self.root / 'recovery'
        audit_name = 'bam_access.1.json'
        reason = f'Unique genomic coverage budget exceeded; genome tiling is prohibited; see {audit_name}'
        calls = []

        def native(sample, runtime, entry, config, out):
            path = Path(out) / runtime['id']
            path.mkdir(parents=True)
            cns = path / 'toy.cns'
            cns.write_text('SYNTHETIC')
            receipt = path / 'complete.json'
            atomic_json(receipt, dict(cns=str(cns), sha256=digest(cns)))
            return dict(cnvkit_version='0.9.13', cns=str(cns), reuse='native_global_batch',
                        evidence_level='native_verified', provenance=str(receipt))

        def seed_filter(sample, runtime, entry, config, cns, out):
            bed = Path(out) / FINAL_BED
            bed.write_text('chr1\t100\t225\t8.5\n')
            atomic_json(Path(out) / 'upstream_provenance.json', dict(version=runtime['id'], cns=cns,
                        io_policy=config['_io_policy_path'], seed_bed_sha256=digest(bed)))
            return _native_seeds(bed, config['_reference_lengths'])

        def reconstruct(sample, selected, runtime, entry, verified, config, out, measurement):
            version = runtime['id']
            policy = json.loads(Path(config['_io_policy_path']).read_text())
            limit = policy['max_unique_bp']
            calls.append((version, limit))
            job = Path(out)
            attempt = job / 'attempt_0001'
            attempt.mkdir(parents=True)
            row = dict(sample_id=sid, version=version, run_dir=str(attempt), labels=[],
                       analysis_kind='aa_ac', source_version=self.source, features=[], intervals=[],
                       amplicons=[], supported_labels=[], reconstruction_intervals=[])
            stages = attempt / 'stages.json'
            bed = attempt / 'selected_seeds.bed'
            bed.write_text('chr1\t100\t225\t8.5\n')
            if version == self.target and limit == 100:
                atomic_json(stages, {'aa': dict(status='failed', returncode=1, reason=reason)})
                atomic_json(attempt / audit_name, dict(blocked=[reason.split('; see')[0]], unique_requested_bp=100,
                    effective_max_unique_bp=100, first_rejected_unique_request=dict(limit_bp=100, projected_unique_bp=120)))
                return dict(row, status='failed', stage='aa', reason=reason, comparison_eligible=False)
            atomic_json(stages, {'aa': dict(status='complete', returncode=0),
                                 'ac': dict(status='skipped_zero_amplicons')})
            summary = attempt / 'reconstruction_summary.txt'
            summary.write_text('#Amplicons = 0\n')
            joint = dict(row, status='not_detected', stage='aa', detected=False, comparison_eligible=True)
            atomic_json(job / 'cache.json', dict(signature='fake', record=joint,
                files={str(path): digest(path) for path in (stages, bed, summary)}))
            return joint

        with (patch('ascompare.concordance_global_cnv.global_cns', side_effect=native) as cnv,
              patch('ascompare.concordance_global_seeds.global_suite_seeds', side_effect=seed_filter) as seed,
              patch('ascompare.concordance_official_solver.official_solver', return_value=('clarabel', 'synthetic')),
              patch('ascompare.concordance.verify_runtime', return_value={'ready': True}),
              patch('ascompare.concordance_global_seeds.global_seed_measurement', return_value=None),
              patch('ascompare.concordance.run_bounded_downstream', side_effect=reconstruct),
              patch('ascompare.concordance_report._plots', return_value=[]),
              patch('ascompare.concordance_cohort_plot.plot_cohort_compatibility')):
            run_concordance(config, lock, entries, prior)
            files_before = {p: digest(p) for p in prior.rglob('*') if p.is_file()}
            self.config.update(recover_aa_from=str(prior), aa_rescue_unique_bp=140)
            config, lock, entries = self.load()
            run_concordance(config, lock, entries, resumed)
            self.assertEqual((cnv.call_count, seed.call_count), (2, 2))
            self.assertEqual(calls, [(self.source, 100), (self.target, 100), (self.target, 140)])
            self.assertEqual(files_before, {p: digest(p) for p in prior.rglob('*') if p.is_file()})
            result = json.loads((resumed / 'downstream_records.json').read_text())
            self.assertEqual([r['status'] for r in result], ['not_detected', 'not_detected'])
            audit = json.loads((resumed / 'aa_recovery' / (sid + '.json')).read_text())
            self.assertEqual(audit['completed_versions_reused'], [self.source])
            self.assertEqual(audit['retries'][self.target]['status'], 'not_detected')
            (prior / 'downstream_records.json').write_text('[]')
            with self.assertRaisesRegex(ValueError, 'Previous report input changed'):
                run_concordance(config, lock, entries, self.root / 'rejected')
            self.assertEqual((cnv.call_count, seed.call_count, len(calls)), (2, 2, 3))

    def test_missing_required_discovery_contig_is_a_preflight_failure(self):
        (self.genome / 'hg19_noAlt.fa.fai').write_text('chr2\t500\n')
        config,lock,entries=self.load()
        rows=run_concordance(config,lock,entries,self.root/'check',check_only=True)['input_check']
        self.assertEqual(rows[0]['status'],'failed')
        self.assertIn('chr2',rows[0]['reason'])

    def test_original_bed_and_runtime_proof_required_to_call_replay(self):
        seedfile=self.archive/'original.bed'
        seedfile.write_text('chr1\t10\t20\n')
        path=self.archive/'original_run_metadata.json'
        atomic_json(path,{'amplified_intervals_cmd':f'aa --ref hg19 --bam {self.bam} --out {self.archive}/original'})
        sample=dict(self.config['samples'][0],seeds_bed=str(seedfile),seeds=[],
                    source_cnvkit_version='0.9.13',source_cnvkit_evidence=dict(
                        origin='source_result_run_metadata',path=str(path)))
        native=dict(chrom='chr1',start=10,end=20,seed_id='x')
        sample['seeds']=[native]
        sets={self.source:dict(status='ok',cnvkit_version='0.9.13',seeds=[native])}
        joint=dict(status='failed',version=self.source,comparison_eligible=False,sample_id='SYNTHETIC')
        with patch('ascompare.concordance_global_workflow._historical_seed_command_matches',return_value=True):
            _source_evidence(sample,joint,sets,{self.source:joint},[native],self.root/'results',self.config)
        self.assertEqual(sample['global_source_seed_calibration'],'passed')
        self.assertEqual(sample['source_bridge_status'],'unverified')
        sets[self.source]['cnvkit_version']='0.9.10'
        _source_evidence(sample,joint,sets,{},[],self.root/'results',self.config)
        self.assertEqual(sample['global_source_seed_calibration'],'unverified')
