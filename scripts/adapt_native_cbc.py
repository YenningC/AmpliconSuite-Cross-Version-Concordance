#!/usr/bin/env python3
"""Explicit isolated-PuLP CBC architecture adaptation; never edits BFB source."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--native-cbc', type=Path, required=True)
p.add_argument('--library-dir', type=Path)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
import pulp
if pulp.__version__ != '3.3.0':
    p.error('This documented adaptation requires the isolated PuLP 3.3.0 dependency')
target = Path(pulp.PULP_CBC_CMD().path).resolve()
prefix = Path(sys.prefix).resolve()
if prefix not in target.parents:
    p.error('PuLP belongs to another environment; first install a private PuLP 3.3.0 copy in this venv')
source = a.native_cbc.resolve()
if not source.is_file() or source == target:
    p.error('Choose a distinct existing native CBC executable')
h = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
probe = subprocess.run([str(source), '-stop'], capture_output=True, text=True)
if probe.returncode:
    p.error('Native CBC cannot execute: ' + probe.stderr)
backup = target.with_name(target.name + '.upstream_backup')
if backup.exists():
    p.error('An adaptation already exists; inspect its saved manifest instead of overwriting')
original = h(target)
shutil.copy2(target, backup)
shutil.copy2(source, target)
record = {'kind': 'native_cbc_binary', 'description': 'Explicit native CBC executable substitution in isolated PuLP3.3; BFB source unchanged; dependency distribution modified', 'path': str(target), 'sha256': h(target), 'original_sha256': original, 'source_path': str(source), 'source_sha256': h(source), 'version_probe': probe.stdout.strip()}
if a.library_dir:
    record['required_library_path'] = str(a.library_dir.resolve())
a.out.parent.mkdir(parents=True, exist_ok=True)
a.out.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record, indent=2))
