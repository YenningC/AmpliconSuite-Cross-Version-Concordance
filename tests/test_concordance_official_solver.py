from pathlib import Path
import json
import os
import sys
import tempfile
import unittest

from ascompare.concordance_official_solver import official_solver
from ascompare.versions import _digest, sha256


class OfficialSolverTests(unittest.TestCase):
    def _fixture(self, directory, version, validator_code):
        commit = 'b' * 40
        suite = directory / ('AmpliconSuite-pipeline-' + commit[:12])
        (suite / 'paalib').mkdir(parents=True)
        (suite / 'paalib/__init__.py').touch()
        (suite / 'paalib/config_validator.py').write_text(validator_code)
        (suite / 'AmpliconSuite-pipeline.py').touch()
        files = {
            path.relative_to(suite).as_posix(): sha256(path)
            for path in suite.rglob('*.py')
        }
        component = dict(version=version, commit=commit, entrypoint='AmpliconSuite-pipeline.py',
                         repo='https://github.com/AmpliconSuite/AmpliconSuite-pipeline',
                         manifest_sha256=_digest(files))
        (suite / '.ascompare-source.json').write_text(json.dumps(dict(component, files=files)))
        entry = dict(suite=component, aa=dict(entrypoint='src/AmpliconArchitect.py'))
        runtime = dict(suite_script=str(suite / 'AmpliconSuite-pipeline.py'), aa_python=sys.executable,
                       aa_script=str(directory / 'AmpliconArchitect-aaaa/src/AmpliconArchitect.py'))
        return entry, runtime, suite

    def test_earlier_suite_uses_its_aa_default(self):
        self.assertEqual(official_solver({}, {'suite': {'version': '1.5.3'}}, {})[0], 'mosek')

    def test_160_and_161_call_no_argument_validator(self):
        for version in ('1.6.0', '1.6.1'):
            with self.subTest(version=version), tempfile.TemporaryDirectory() as tmp:
                entry, runtime, _ = self._fixture(
                    Path(tmp), version, 'def _validate_mosek_license():\n    return True\n'
                )
                solver, reason = official_solver(runtime, entry, os.environ)
                self.assertEqual(solver, 'mosek')
                self.assertIn('accepted', reason)

    def test_162_uses_checked_aa_interpreter_and_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            entry, runtime, _ = self._fixture(
                Path(tmp), '1.6.2', '\n'.join([
                    'import os',
                    'def _validate_mosek_license(python):',
                    '    assert python == os.environ["EXPECTED_AA_PYTHON"]',
                    '    return False',
                ]) + '\n'
            )
            # Exercise a prepared panel whose runtime lists only the AA/AC scripts.
            del runtime['suite_script']
            runtime['env'] = {'EXPECTED_AA_PYTHON': sys.executable}
            solver, reason = official_solver(runtime, entry, os.environ)
            self.assertEqual(solver, 'clarabel')
            self.assertIn('official Clarabel fallback', reason)

    def test_validator_exception_is_not_treated_as_license_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            entry, runtime, _ = self._fixture(
                Path(tmp), '1.6.0', 'def _validate_mosek_license():\n    raise RuntimeError("private license path")\n'
            )
            with self.assertRaisesRegex(ValueError, 'Could not run the locked Suite') as error:
                official_solver(runtime, entry, os.environ)
            self.assertNotIn('private license path', str(error.exception))

    def test_changed_source_is_rejected_before_running_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            entry, runtime, suite = self._fixture(
                Path(tmp), '1.6.0', 'def _validate_mosek_license():\n    return True\n'
            )
            (suite / 'paalib/config_validator.py').write_text('raise RuntimeError("modified source")\n')
            with self.assertRaisesRegex(ValueError, 'Missing or changed source file'):
                official_solver(runtime, entry, os.environ)
