"""Synthetic execution contracts for the fixed-source-AA CNVkit comparisons."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ascompare.common import atomic_json, digest
from ascompare.concordance import _resolve
from ascompare.local_cnv_comparator import comparison_panel, run_local_cnv_comparators

SOURCE = 'suite-1.6.1'
OTHER = 'suite-1.3.3'
ALSO = 'suite-1.2.2'
REF = 'suite-1.6.2'


class LocalComparatorTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name).resolve()
        self.helper = self.root / 'amplified_intervals.py'
        self.helper.write_text('# SYNTHETIC\n')
        self.source = dict(id=SOURCE, aa_python='/source/aa/python', aa_script=str(self.root/'AA.py'),
                           cnvkit_version='0.9.13', cnvkit_python='/v13/python', cnvkit_script='/v13/cnvkit',
                           env={'PROFILE': 'source_AA'})
        self.older = dict(self.source, id=OTHER, aa_python='/old/aa/python',
                          cnvkit_version='0.9.10', cnvkit_python='/v10/python', cnvkit_script='/v10/cnvkit',
                          env={'PROFILE': 'older_CNV'})
        self.runtimes = {r['id']: r for r in (self.source, self.older,
                                             dict(self.older, id=ALSO), dict(self.source, id=REF))}
        self.config = dict(runtimes=list(self.runtimes.values()), target_versions=[OTHER, ALSO, REF],
                           cnvkit_runtimes=[dict(version='0.9.14', python='/v14/python', script='/v14/cnvkit',
                                                 env={'PROFILE':'v14_CNV'})])
        self.sample = dict(sample_id='SYNTHETIC_1', source_version=SOURCE, source_cnvkit_version='0.9.14',
                           seeds=[dict(seed_id='a', chrom='chr1', start=0, end=10)])
        self.entry = {SOURCE: dict(aa=dict(version='1.6.r1', entrypoint='AA.py', manifest_sha256='m', commit='aa'))}
        self.probe = {SOURCE: dict(aa=dict(source_root=str(self.root), sha256='entrydigest',
                                           python=self.source['aa_python'], **self.entry[SOURCE]['aa']))}
        self.called = []

    def fake_measure(self, sample, rt, entry, config, attempt):
        attempt.mkdir(parents=True, exist_ok=True)
        version = rt['cnvkit_version']; sid = sample['sample_id']
        self.called.append((sid, version, rt['aa_python'], rt['env'], rt.get('cnvkit_env')))
        key = version + sid
        prov = attempt / 'upstream_provenance.json'
        ident = {'version':version}
        atomic_json(prov, dict(version=SOURCE, cnvkit=ident, amplified_intervals=dict(path=str(self.helper),
            sha256=digest(self.helper), args=['--gain','4.5','--cnsize_min','50000','--ref','hg19']),
            cache=dict(key=key, signature=dict(sample_id=sid, inputs={'alignment':sid}, io_contract={'max':10},
                                               parameters={'windows':[[1,2,3]]}, cnvkit=ident))))
        return [dict(sample_id=sid, source_version=SOURCE, version=SOURCE, cnvkit_version=version, seed_id='a',
                     seed=sample['seeds'][0], retention_fraction=.5, retained_intervals=[['chr1',0,5]],
                     status='not_retained', min_seed_overlap=1., comparison_eligible=True, stage='seed_processing',
                     cnv_cache_key=key, provenance=str(prov))]

    def run_case(self, native=None):
        with patch('ascompare.local_cnv_comparator.run_local_upstream', self.fake_measure):
            return run_local_cnv_comparators(self.sample,self.config,self.entry,self.runtimes,self.probe,
                                             native or [],self.root/'measures')

    def test_control_and_targets_use_fixed_filter_and_no_duplicate_cnv_arm(self):
        native = self.fake_measure(self.sample,self.source,{},self.config,self.root/'native')
        self.called.clear()
        rows = self.run_case(native)
        self.assertEqual([r['cnvkit_version'] for r in rows], ['0.9.14','0.9.13','0.9.10'])
        self.assertEqual([r['measurement_origin'] for r in rows],
                         ['controlled_local_job','native_source_job','controlled_local_job'])
        self.assertEqual({r['source_filter_aa_version'] for r in rows},{'1.6.r1'})
        self.assertEqual([(v,aa,e['PROFILE'],ce['PROFILE']) for _,v,aa,e,ce in self.called],
                         [('0.9.14','/source/aa/python','source_AA','v14_CNV'),
                          ('0.9.10','/source/aa/python','source_AA','older_CNV')])
        self.assertEqual(rows[-1]['target_suite_versions'],[OTHER,ALSO])
        self.assertEqual(len({r['local_design_id'] for r in rows}),1)
        self.assertTrue(all(r['local_design_id'] and r['actual_cnvkit_version']==r['cnvkit_version'] for r in rows))
        self.assertEqual(rows[1]['calibration_role'],'candidate_observation')

    def test_each_patient_uses_own_standard_and_data(self):
        a=self.run_case()
        self.sample.update(sample_id='SYNTHETIC_2',source_cnvkit_version='0.9.10')
        b=self.run_case()
        self.assertEqual([r['cnvkit_version'] for r in b if r['measurement_role']=='source_local_standard'],['0.9.10'])
        self.assertFalse({r['local_design_id'] for r in a}&{r['local_design_id'] for r in b})

    def test_unmatched_explicit_runtime_does_not_inherit_source_AA_loader(self):
        entry=self.config['cnvkit_runtimes'][0]
        entry.pop('env')
        entry['extra_path']=['/tested/R/bin']
        from ascompare.local_cnv_comparator import _choose_tool
        tool, reason=_choose_tool('0.9.14',self.config,self.runtimes,self.source)
        self.assertIsNone(reason)
        self.assertEqual(tool['env'],{})
        self.assertEqual(tool['extra_path'],['/tested/R/bin'])
        self.run_case()
        self.assertEqual(self.called[0][-1],{})

    def test_missing_source_environment_is_not_replaced_and_distinct_target_jobs_skip(self):
        self.config['cnvkit_runtimes']=[]
        rows=self.run_case()
        self.assertFalse(self.called)
        self.assertEqual(rows[0]['status'],'missing')
        self.assertEqual(rows[0]['source_cnvkit_version'],'0.9.14')
        self.assertTrue(all(not r['comparison_eligible'] and not r['actual_cnvkit_version'] for r in rows))
        self.assertEqual(rows[-1]['reason'],'source_local_standard_unavailable')

    def test_conflicting_environment_fails_that_arm_only(self):
        self.runtimes[ALSO]=dict(self.older,id=ALSO,cnvkit_python='/alternate/python')
        self.config['runtimes']=list(self.runtimes.values())
        result=self.run_case()
        self.assertTrue(result[0]['comparison_eligible'])
        self.assertEqual(result[-1]['status'],'unsupported')
        self.assertEqual(result[-1]['reason'],'multiple_configured_CNVkit_environments')

    def test_missing_aa_and_source_conflict_are_unavailable(self):
        self.probe[SOURCE]['aa']['version']='wrong'
        self.assertEqual(self.run_case()[0]['reason'],'source_AA_filter_not_verified')
        self.sample.update(_source_cnvkit_issue='record conflicts')
        self.assertTrue(all(r['reason']=='record conflicts' for r in self.run_case()))
        self.assertFalse(self.called)

    def test_tampered_or_missing_design_never_satisfies_standard(self):
        self.sample['source_cnvkit_version']='0.9.13'
        rows=self.fake_measure(self.sample,self.source,{},self.config,self.root/'native')
        prov=Path(rows[0]['provenance'])
        data=json.loads(prov.read_text());data['cnvkit']['version']='0.9.14';atomic_json(prov,data)
        self.called.clear()
        result=self.run_case(rows)
        self.assertIsNone(result[0]['actual_cnvkit_version'])
        self.assertEqual(result[0]['status'],'unsupported')
        self.assertFalse(self.called)

    def test_panel_retains_suite_mapping_without_claiming_installed_runtime(self):
        self.config['runtimes']=[]
        self.assertEqual(comparison_panel(self.sample,self.config),
                         {'0.9.14':[], '0.9.13':[SOURCE,REF], '0.9.10':[OTHER,ALSO]})

    def test_optional_versions_are_included_once_independent_of_runtime_catalog(self):
        self.config['local_cnv'] = {'comparison_versions':['0.9.13','0.9.12','0.9.12']}
        self.assertEqual(comparison_panel(self.sample,self.config)['0.9.12'],[])
        self.assertEqual(list(comparison_panel(self.sample,self.config)).count('0.9.13'),1)
        result=self.run_case()
        extra = next(r for r in result if r['cnvkit_version']=='0.9.12')
        self.assertEqual((extra['status'],extra['target_suite_versions']),('missing',[]))

    def test_registry_resolves_venv_python_without_following_symlink(self):
        (self.root/'base').touch()
        (self.root/'venv-python').symlink_to(self.root/'base')
        entry=dict(version='0.9.14',python='venv-python',script='bin/cnvkit',extra_path=['bin'])
        config=_resolve(dict(cnvkit_runtimes=[entry]),self.root)
        self.assertEqual(config['cnvkit_runtimes'][0]['python'],str(self.root/'venv-python'))
        self.assertEqual(config['cnvkit_runtimes'][0]['script'],str(self.root/'bin/cnvkit'))
        self.assertEqual(config['cnvkit_runtimes'][0]['extra_path'],[str(self.root/'bin')])
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            _resolve(dict(cnvkit_runtimes=[entry,entry]),self.root)
