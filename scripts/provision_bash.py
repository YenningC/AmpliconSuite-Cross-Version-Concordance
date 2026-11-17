#!/usr/bin/env python3
"""Provision an isolated, hash-locked ARM64 Bash tool without changing an environment.

Uses package extraction only (no conda environment creation or global registration).
Applies each package's declared prefix relocation and repairs local Mach-O
signatures after relocation. Only the Bash executable is exposed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import urllib.request

PACKAGES = [
    ('bash-5.3-hea35128_0.conda', '569d0e58c60cb9282182dc648a779d245d41dd5e944c55292803450edba3a833'),
    ('readline-8.3-h8b90a29_1.conda', 'd6782b430e6dce4fe8c7ac118cd988d5a193eb4a7f60741edfaede58016b6110'),
    ('ncurses-6.6-he64c551_1.conda', '7024a48c8c0d0114ed4ab53c76bf9275d50e91ba7cea367a9aead638d3c29c68'),
]
BASE = 'https://conda.anaconda.org/conda-forge/osx-arm64/'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def binary_prefix_replace(data, old, new):
    if len(new)>len(old):
        raise ValueError('Installation prefix exceeds package placeholder length')
    pattern=re.compile(re.escape(old)+b'[^\x00]*\x00')
    def replacement(match):
        original=match.group(0)
        changed=original[:-1].replace(old,new)
        return changed+b'\x00'*(len(original)-len(changed))
    return pattern.sub(replacement,data)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--micromamba', required=True, type=Path)
    ap.add_argument('--cache-dir', type=Path)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--download', action='store_true', help='Allow public package download when absent from cache')
    args = ap.parse_args()
    if platform.system() != 'Darwin' or platform.machine() != 'arm64':
        ap.error('This exact package lock is for macOS ARM64 only')
    out = args.out.resolve()
    if out.exists():
        ap.error('Output must be a new directory; existing tools/environments are never overwritten')
    out.mkdir(parents=True)
    (out/'archives').mkdir(); (out/'packages').mkdir()
    prefix = out/'runtime'; (prefix/'bin').mkdir(parents=True); (prefix/'lib').mkdir()
    records = []; relocations=[]; signing=[]
    for filename, expected in PACKAGES:
        archive = out/'archives'/filename
        cached = args.cache_dir/filename if args.cache_dir else None
        if cached and cached.is_file():
            shutil.copy2(cached, archive)
            acquisition = {'method': 'local_cache_copy', 'cache': str(cached.resolve())}
        elif args.download:
            with urllib.request.urlopen(BASE+filename, timeout=60) as response, archive.open('wb') as dest:
                shutil.copyfileobj(response, dest)
            acquisition = {'method': 'public_download'}
        else:
            raise RuntimeError('Missing cached package; --download is required to fetch: '+filename)
        if sha(archive) != expected:
            raise RuntimeError('Package checksum mismatch: '+filename)
        extracted = out/'packages'/filename.removesuffix('.conda')
        command = [str(args.micromamba.resolve()), 'package', 'extract', '--no-rc', '--offline', str(archive), str(extracted)]
        proc = subprocess.run(command, text=True, capture_output=True)
        (out/(filename+'.extract.log')).write_text(proc.stdout+proc.stderr)
        proc.check_returncode()
        index = json.loads((extracted/'info/index.json').read_text())
        records.append({'source_url':BASE+filename, 'archive_sha256':expected, 'index':index, 'acquisition':acquisition, 'extract_argv':command})
        if (extracted/'lib').exists():
            shutil.copytree(extracted/'lib', prefix/'lib', dirs_exist_ok=True, symlinks=True)
        if index['name'] == 'bash':
            shutil.copy2(extracted/'bin/bash', prefix/'bin/bash')
        for entry in json.loads((extracted/'info/paths.json').read_text())['paths']:
            dest=prefix/entry['_path']
            placeholder=entry.get('prefix_placeholder')
            if not placeholder or not dest.is_file() or dest.is_symlink():
                continue
            before=dest.read_bytes()
            if entry.get('file_mode')=='binary':
                after=binary_prefix_replace(before,placeholder.encode(),str(prefix).encode())
            else:
                after=before.replace(placeholder.encode(),str(prefix).encode())
            if after!=before:
                dest.write_bytes(after)
                relocations.append({'file':str(dest),'package':filename,'before_sha256':hashlib.sha256(before).hexdigest(),'after_sha256':sha(dest),'placeholder':placeholder,'prefix':str(prefix)})
    for dest in [prefix/'bin/bash']+[x for x in (prefix/'lib').rglob('*') if x.is_file() and not x.is_symlink()]:
        kind=subprocess.run(['/usr/bin/file','-b',str(dest)],text=True,capture_output=True,check=True).stdout
        if 'Mach-O' not in kind:
            continue
        verified=subprocess.run(['/usr/bin/codesign','--verify',str(dest)],text=True,capture_output=True)
        if verified.returncode:
            before=sha(dest);command=['/usr/bin/codesign','--force','--sign','-',str(dest)]
            signed=subprocess.run(command,text=True,capture_output=True,check=True)
            signing.append({'file':str(dest),'before_sha256':before,'after_sha256':sha(dest),'argv':command,'stderr':signed.stderr})
    bash = prefix/'bin/bash'
    env = dict(os.environ)
    for name in ('DYLD_LIBRARY_PATH','DYLD_FALLBACK_LIBRARY_PATH','LD_LIBRARY_PATH','BASH_ENV','ENV'):
        env.pop(name, None)
    probe_script = 'declare -A probe || exit 1; probe[key]=ok; test "${probe[key]}" = ok; printf "associative_arrays=ok\\n"'
    checks = {}
    (out/'link_probe').mkdir(); linked = out/'link_probe/bash'; linked.symlink_to(bash)
    for name, command in [('version',[str(bash),'--version']), ('associative_arrays',[str(bash),'-c',probe_script]), ('cross_directory_symlink',[str(linked),'-c',probe_script])]:
        result = subprocess.run(command, text=True, capture_output=True, env=env)
        checks[name] = {'argv':command,'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr}
        (out/'runtime_probes.json').write_text(json.dumps(checks,indent=2)+'\n')
        result.check_returncode()
    inventory = {str(p.relative_to(prefix)):sha(p) for p in prefix.rglob('*') if p.is_file() and not p.is_symlink()}
    provenance = {'platform':'osx-arm64','exposed_bin':str(prefix/'bin'),'packages':records,'bash_sha256':sha(bash),'binary_inventory':inventory,'checks':checks,'relocations':relocations,'local_signing':signing,'method':'Extract verified packages into isolated tool tree, apply declared prefix relocation and local Mach-O signing; no global environment registration or existing environment mutation'}
    (out/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(json.dumps({'bin':str(prefix/'bin'),'bash_sha256':sha(bash),'provenance':str(out/'provenance.json'),'version':checks['version']['stdout'].splitlines()[0]}))


if __name__ == '__main__':
    main()
