"""Portable panel construction and report integrity checks; no patient inputs."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ascompare.benchmark import main
from ascompare.versions import load_lock
from scripts.build_benchmark_manifest import build


class BenchmarkCliTests(unittest.TestCase):
    def fixture(self):
        lock = load_lock()
        prepared = dict(reference_version=lock['reference_version'], sources=[
            dict(id=e['id'], aa_script='/sources/'+e['id']+'/AA.py',
                 ac_script='/sources/'+e['id']+'/AC.py') for e in lock['versions'] if e.get('suite')])
        templates = [dict(id=vid, aa_python=py, ac_python=py, env={'MOSEKLM_LICENSE_FILE':'/secure/license'})
                     for vid, py in [('suite-1.6.2','/current/python'),('suite-1.6.1','/bfb101/python')]]
        cnv = [dict(version=v, python='/cnv/'+v+'/python', script='/cnv/'+v+'/cnvkit.py')
               for v in ('0.9.10','0.9.13')]
        manifest = dict(mode='with_bam', samples=[dict(source_version='suite-1.3.4',
                            source_result={'aa_dir':'relative/aa','ac_dir':'relative/ac'})])
        return manifest, prepared, templates, cnv, lock

    def test_exact_panel_mapping_and_inputs_are_not_mutated(self):
        args = self.fixture(); before = copy.deepcopy(args)
        result = build(*args)
        self.assertEqual(args, before)
        self.assertEqual(len(result['runtimes']), 10)
        self.assertEqual(result['reference_version'], 'suite-1.6.2')
        for runtime in result['runtimes']:
            modern = runtime['id'].startswith('suite-1.6.')
            self.assertEqual(runtime['cnvkit_version'], '0.9.13' if modern else '0.9.10')
            self.assertEqual(runtime['solver'], 'clarabel' if modern else 'mosek')
            self.assertEqual(runtime['ac_python'], '/bfb101/python' if runtime['id'] in
                             ('suite-1.6.0','suite-1.6.1') else '/current/python')
            self.assertIn(runtime['id'], runtime['aa_script'])
        self.assertEqual(result['samples'][0]['source_result']['aa_dir'],'relative/aa')

    def test_missing_historical_cnvkit_is_not_replaced_by_newest(self):
        manifest, prepared, templates, cnv, lock = self.fixture()
        with self.assertRaisesRegex(ValueError, '0.9.10'):
            build(manifest, prepared, templates, cnv[1:], lock)
        manifest['mode']='results_only'
        self.assertTrue(all('cnvkit_script' not in r for r in
                           build(manifest,prepared,templates,[],lock)['runtimes']))

    def test_report_refuses_modified_saved_results(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary)/'demo'
            with patch('ascompare.benchmark_report._plots'):
                self.assertEqual(main(['demo','--out',str(out)]),0)
                self.assertEqual(main(['report','--out',str(out)]),0)
                path=out/'downstream_records.json'
                records=json.loads(path.read_text()); records[0]['labels']=['fabricated']
                path.write_text(json.dumps(records))
                self.assertEqual(main(['report','--out',str(out)]),2)


if __name__ == '__main__':
    unittest.main()
