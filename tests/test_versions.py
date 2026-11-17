import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from ascompare.versions import _check_source, _digest, load_lock, verify_runtime


class VersionTests(unittest.TestCase):
    def test_fixed_reference_and_component_fingerprints(self):
        lock = load_lock()
        self.assertEqual(lock['reference_version'], 'suite-1.6.2')
        self.assertGreaterEqual(sum(bool(v.get('suite')) for v in lock['versions']), 11)
        transition = next(v for v in lock['versions'] if v['id'] == 'suite-1.5.3')
        self.assertEqual((transition['aa']['version'], transition['ac']['version']), ('1.5.r7', '1.5.2'))
        self.assertFalse(transition.get('bfb'))
        for version in lock['versions']:
            for name in ('aa', 'ac', 'bfb', 'suite'):
                if version.get(name):
                    c = version[name]
                    self.assertEqual(len(c['commit']), 40)
                    self.assertEqual(len(c['archive_sha256']), 64)
                    self.assertEqual(len(c['manifest_sha256']), 64)

    def test_changed_source_cannot_keep_identity(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            script = root / 'tool.py'
            script.write_text('print("historical")\n')
            comp = {'commit': 'a' * 40, 'repo': 'https://github.com/AmpliconSuite/AmpliconClassifier'}
            files = {'tool.py': hashlib.sha256(script.read_bytes()).hexdigest()}
            comp['manifest_sha256'] = _digest(files)
            (root / '.ascompare-source.json').write_text(json.dumps(dict(comp, files=files)))
            self.assertEqual(_check_source(root, comp)['files_checked'], 1)
            script.write_text('print("replacement")\n')
            with self.assertRaisesRegex(ValueError, 'changed source'):
                _check_source(root, comp)

    def test_extra_python_cannot_shadow_code(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            comp = {'commit': 'a' * 40, 'repo': 'https://github.com/AmpliconSuite/AmpliconClassifier',
                    'manifest_sha256': _digest({})}
            (root / '.ascompare-source.json').write_text(json.dumps(dict(comp, files={})))
            (root / 'replacement.py').write_text('pass')
            with self.assertRaisesRegex(ValueError, 'Unexpected Python'):
                _check_source(root, comp)

    def test_identity_mismatch_never_falls_back(self):
        result = verify_runtime({'id': 'old'}, {'id': 'new'})
        self.assertFalse(result['ready'])
        self.assertEqual(result['diagnostics'], ['Runtime/lock ID mismatch'])
