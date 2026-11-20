"""Synthetic public workflow and source reproduction regressions."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pysam

from ascompare.common import atomic_json
from ascompare.concordance import load_manifest, run_concordance
from ascompare.concordance_global_workflow import _source_evidence


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
              patch('ascompare.concordance_global_seeds.global_suite_seeds',side_effect=[[],seeds]),
              patch('ascompare.concordance_official_solver.official_solver',return_value=('clarabel','synthetic')),
              patch('ascompare.concordance.verify_runtime',return_value={'ready':True}),
              patch('ascompare.concordance_global_seeds.global_seed_measurement',return_value=None),
              patch('ascompare.concordance.run_bounded_downstream',side_effect=reconstruct),
              patch('ascompare.concordance_report._plots',return_value=[]),
              patch('ascompare.concordance_cohort_plot.plot_cohort_compatibility')):
            run_concordance(config,lock,entries,out)
        self.assertEqual(cnv.call_count,2)
        self.assertEqual(observed,[(self.target,seeds)])
        records=json.loads((out/'upstream_records.json').read_text())
        self.assertEqual({r['version']:r['status'] for r in records},
                         {self.source:'not_retained',self.target:'retained'})
        sets=json.loads((out/'global_upstream/SYNTHETIC/discovery.json').read_text())
        self.assertTrue(sets[self.target]['historical_provenance_limited'])
        self.assertTrue((out/'result.tsv').is_file())

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
