import copy
import unittest
from ascompare.versions import load_lock
from scripts.build_panel_config import build_config


class PanelConfigTests(unittest.TestCase):
    def fixture(self):
        lock=load_lock()
        sources=[{'id':e['id'],'suite_script':'/src/'+e['id']+'/Suite.py','aa_script':'/src/'+e['id']+'/AA.py','ac_script':'/src/'+e['id']+'/AC.py'} for e in lock['versions']]
        templates=[{'id':vid,'suite_python':py,'aa_python':py,'ac_python':py,'cnvkit_script':py.replace('python','cnvkit.py'),'cnvkit_version':'0.9.14','extra_path':['/base/bin']} for vid,py in [('suite-1.6.2','/reference-env/bin/python'),('suite-1.6.1','/ac200-env/bin/python')]]
        return lock,{'reference_version':lock['reference_version'],'sources':sources},templates

    def test_official_panel_uses_correct_sources_environments_and_solvers(self):
        lock,prepared,templates=self.fixture()
        config=build_config(prepared,templates,lock,{'build':'hg19','fasta':'/ref/x.fa','aa_data_repo':'/ref'},mosek_license='/secure/mosek.lic')
        self.assertEqual(len(config['versions']),11)
        self.assertEqual(config['versions'][0]['id'],'suite-1.6.2')
        for v in config['versions']:
            self.assertTrue(v['id'].startswith('suite-'))
            self.assertIn(v['id'],v['aa_script'])
            expected='/ac200-env/bin/python' if v['id'] in ('suite-1.6.0','suite-1.6.1') else '/reference-env/bin/python'
            self.assertEqual(v['ac_python'],expected)
            self.assertEqual(v['solver'],'clarabel' if v['id'].startswith('suite-1.6.') else 'mosek')
            self.assertEqual(v['env']['MOSEKLM_LICENSE_FILE'],'/secure/mosek.lic')
        self.assertEqual(config['statistics']['margins'],{})

    def test_missing_sources_or_reference_mismatch_cannot_silently_substitute(self):
        lock,prepared,templates=self.fixture()
        missing=copy.deepcopy(prepared);missing['sources']=[x for x in missing['sources'] if x['id']!='suite-1.3.4']
        with self.assertRaisesRegex(ValueError,'suite-1.3.4'):build_config(missing,templates,lock,{})
        prepared['reference_version']='suite-1.6.1'
        with self.assertRaisesRegex(ValueError,'different fixed reference'):build_config(prepared,templates,lock,{})
