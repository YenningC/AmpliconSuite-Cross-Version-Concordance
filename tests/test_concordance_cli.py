"""Portable panel construction and report integrity checks; no patient inputs."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ascompare.concordance import main
from ascompare.versions import load_lock
from scripts.build_concordance_manifest import build


class ConcordanceCliTests(unittest.TestCase):
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
        self.assertEqual(len(result['runtimes']), 11)
        self.assertEqual(result['reference_version'], 'suite-1.6.2')
        for runtime in result['runtimes']:
            modern = runtime['id'].startswith('suite-1.6.')
            uses_cnvkit_0913 = modern or runtime['id'] == 'suite-1.5.3'
            self.assertEqual(runtime['cnvkit_version'], '0.9.13' if uses_cnvkit_0913 else '0.9.10')
            self.assertEqual(runtime['solver'], 'clarabel' if modern else 'mosek')
            self.assertEqual(runtime['ac_python'], '/bfb101/python' if runtime['id'] in
                             ('suite-1.6.0','suite-1.6.1') else '/current/python')
            self.assertIn(runtime['id'], runtime['aa_script'])
        self.assertEqual(result['samples'][0]['source_result']['aa_dir'],'relative/aa')

    def test_build_keeps_original_version_even_when_no_target_recipe_uses_it(self):
        manifest, prepared, templates, cnv, lock = self.fixture()
        original = {'version':'0.9.14', 'python':'/extra/bin/python', 'script':'/extra/bin/cnvkit.py'}
        manifest['samples'][0]['source_cnvkit_version'] = '0.9.14'
        cnv[0].update(env={}, extra_path=['/tools/R'])
        result = build(manifest, prepared, templates, cnv+[original], lock)
        native = next(r for r in result['runtimes'] if r['cnvkit_version'] == cnv[0]['version'])
        self.assertEqual(native['cnvkit_env'], {})
        self.assertEqual(native['cnvkit_extra_path'], ['/tools/R'])
        self.assertIn('MOSEKLM_LICENSE_FILE', native['env'])
        self.assertIn(original, result['cnvkit_runtimes'])
        self.assertNotIn('0.9.14', [r['cnvkit_version'] for r in result['runtimes']])
        manifest['cnvkit_runtimes'] = [dict(original, python='/other/bin/python')]
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            build(manifest, prepared, templates, cnv+[original], lock)

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
            with patch('ascompare.concordance_report._plots'):
                self.assertEqual(main(['demo','--out',str(out)]),0)
                self.assertTrue(json.loads((out/'provenance.json').read_text())['synthetic'])
                self.assertEqual(main(['report','--out',str(out)]),0)
                detached = Path(temporary)/'report-copy'
                original = {p.name: p.read_bytes() for p in out.glob('*records.json')}
                self.assertEqual(main(['report', '--from-run', str(out), '--out', str(detached)]), 0)
                self.assertEqual(original, {p.name: p.read_bytes() for p in out.glob('*records.json')})
                self.assertTrue((detached/'result.tsv').is_file())
                self.assertFalse(json.loads((detached/'report_provenance.json').read_text())['algorithms_rerun'])
                self.assertEqual(main(['report', '--from-run', str(out), '--out', str(detached)]), 2)
                local_record = out/'local_cnv_comparator_records.json'
                self.assertIn(local_record.name, original)
                local_record.write_text('[{"extra": true}]')
                self.assertEqual(main(['report', '--out', str(out)]), 2)
                local_record.write_bytes(original[local_record.name])
                path=out/'downstream_records.json'
                records=json.loads(path.read_text()); records[0]['labels']=['fabricated']
                path.write_text(json.dumps(records))
                self.assertEqual(main(['report','--out',str(out)]),2)


if __name__ == '__main__':
    unittest.main()
