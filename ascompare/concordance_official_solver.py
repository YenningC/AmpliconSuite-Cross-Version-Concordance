"""Resolve the copy-number solver with the locked Suite's license check."""

from pathlib import Path
import subprocess

from .versions import _check_source


_MODERN_SOLVER_PREFLIGHT = {
    '1.6.0': False,
    '1.6.1': False,
    '1.6.2': True,
}

_CHECK = """
import logging
import sys

logging.disable(logging.CRITICAL)
sys.path.insert(0, sys.argv[1])
from paalib.config_validator import _validate_mosek_license

args = [sys.argv[3]] if sys.argv[2] == 'with_interpreter' else []
print('licensed' if _validate_mosek_license(*args) else 'unlicensed')
"""


def _suite_root(runtime, entry):
    suite = entry['suite']
    if runtime.get('suite_script'):
        script = Path(runtime['suite_script']).resolve()
        root = script
        for _ in Path(suite['entrypoint']).parts:
            root = root.parent
        if script != root / suite['entrypoint']:
            raise ValueError('Suite script does not match the locked entrypoint')
        return root
    if runtime.get('suite_source'):
        return Path(runtime['suite_source']).resolve()

    # prepare-panel installs the verified component snapshots side by side.
    aa_root = Path(runtime['aa_script']).resolve()
    for _ in Path(entry['aa']['entrypoint']).parts:
        aa_root = aa_root.parent
    return aa_root.parent / ('AmpliconSuite-pipeline-' + suite['commit'][:12])


def official_solver(runtime, entry, environment):
    """Return the official solver to request and a provenance reason.

    Newer Suite releases check the license locally before launching AA. AA's
    own ability to retry a failed MOSEK optimization then remains in effect.
    The subprocess executes only the locked license check. Its diagnostic text
    is intentionally discarded because it can contain local license paths.
    """
    suite_version = (entry.get('suite') or {}).get('version')
    if suite_version not in _MODERN_SOLVER_PREFLIGHT:
        return 'mosek', 'Locked AA default (no Suite solver selection)'

    root = _suite_root(runtime, entry)
    _check_source(root, entry['suite'])
    aa_python = str(Path(runtime['aa_python']).absolute())
    argument_mode = 'with_interpreter' if _MODERN_SOLVER_PREFLIGHT[suite_version] else 'without_interpreter'
    try:
        result = subprocess.run(
            [aa_python, '-I', '-B', '-c', _CHECK, str(root), argument_mode, aa_python],
            env=dict(environment, **runtime.get('env', {})),
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True, timeout=30, check=True,
        )
    except (OSError, subprocess.SubprocessError):
        raise ValueError(f'Could not run the locked Suite {suite_version} MOSEK license check with the configured AA interpreter') from None

    if result.stdout.strip() == 'licensed':
        return 'mosek', f'Suite {suite_version} native MOSEK license check accepted'
    if result.stdout.strip() == 'unlicensed':
        return 'clarabel', f'Suite {suite_version} native MOSEK license check selected the official Clarabel fallback'
    raise ValueError(f'Locked Suite {suite_version} license check returned an unexpected result')
