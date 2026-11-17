"""Official release locks, verified source acquisition, and runtime checks.

No image-name inference, automatic version substitution, or remote patient access.
"""
from __future__ import annotations
import base64
import hashlib
import io
import json
import os
import posixpath
from pathlib import Path
import re
import subprocess
import tarfile
import urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
REPOS = ('AmpliconSuite-pipeline', 'AmpliconArchitect', 'AmpliconClassifier', 'BFBArchitect')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def _digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def load_lock(path=None):
    default = ROOT / 'versions.lock.json'
    if not default.is_file():
        default = Path(__file__).parent / 'versions.lock.json'
    data = json.loads(Path(path or default).read_text())
    ids = [v['id'] for v in data['versions']]
    if len(ids) != len(set(ids)) or data['reference_version'] not in ids:
        raise ValueError('Invalid version lock identifiers/reference')
    return data


def _get(url):
    # All URLs are public software repositories; no sample information is sent.
    if not url.startswith(('https://api.github.com/repos/AmpliconSuite/',
                           'https://codeload.github.com/AmpliconSuite/')):
        raise ValueError('Only official AmpliconSuite source URLs are accepted')
    req = urllib.request.Request(url, headers={'User-Agent': 'ampliconsuite-compare/0.1'})
    error = None
    for _ in range(3):
        try:
            with urllib.request.urlopen(req, timeout=45) as r:
                return r.read()
        except OSError as exc:
            error = exc
    raise error


def refresh_releases(destination):
    """Save public release history without modifying the study's fixed lock."""
    dest = Path(destination)
    dest.mkdir(parents=True, exist_ok=True)
    result = {}
    for repo in REPOS:
        data = json.loads(_get(f'https://api.github.com/repos/AmpliconSuite/{repo}/releases?per_page=100'))
        (dest / (repo + '.json')).write_text(json.dumps(data, indent=2) + '\n')
        stable = [r for r in data if not r['draft'] and not r['prerelease']]
        result[repo] = {'latest_stable': max(stable, key=lambda r: r['published_at'])['tag_name'] if stable else None,
                        'release_count': len(data)}
    result['checked_at'] = datetime.now(timezone.utc).isoformat()
    (dest / 'discovery.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def _check_source(root, comp):
    root = Path(root)
    mf = root / '.ascompare-source.json'
    if not mf.is_file():
        raise ValueError(f'Missing verified source manifest: {mf}; run prepare first')
    manifest = json.loads(mf.read_text())
    if manifest['commit'] != comp['commit'] or manifest['repo'] != comp['repo']:
        raise ValueError('Source manifest commit/repository mismatch')
    files = manifest['files']
    if comp.get('manifest_sha256') != _digest(files):
        raise ValueError('Source inventory does not match the fixed release lock')
    for rel, digest in files.items():
        fp = root / rel
        if not fp.is_file() or sha256(fp) != digest:
            raise ValueError(f'Missing or changed source file: {fp}')
    # Generated caches are harmless; unexpected code can shadow locked imports.
    extra = {str(f.relative_to(root)) for f in root.rglob('*.py')} - set(files)
    if extra:
        raise ValueError(f'Unexpected Python source: {sorted(extra)[:5]}')
    return {'source_root': str(root), 'commit': comp['commit'],
            'manifest_sha256': _digest(files), 'files_checked': len(files)}


def prepare_component(comp, destination, allow_new_digest=False):
    """Acquire one immutable commit. Normal calls require recorded archive digest."""
    if not re.fullmatch(r'[0-9a-f]{40}', comp['commit']):
        raise ValueError('A full resolved commit is required')
    dest = Path(destination).resolve()
    repo = comp['repo'].rsplit('/', 1)[-1]
    if repo not in REPOS or comp['repo'] != 'https://github.com/AmpliconSuite/' + repo:
        raise ValueError('Unexpected repository')
    target = dest / (repo + '-' + comp['commit'][:12])
    archive = dest / (repo + '-' + comp['commit'] + '.tar.gz')
    if target.exists() and comp.get('manifest_sha256'):
        _check_source(target, comp)
        return target
    dest.mkdir(parents=True, exist_ok=True)
    if not archive.exists():
        content = _get(f'https://codeload.github.com/AmpliconSuite/{repo}/tar.gz/{comp["commit"]}')
        archive.write_bytes(content)
    actual = sha256(archive)
    if comp.get('archive_sha256') != actual and not allow_new_digest:
        raise ValueError(f'Archive fingerprint missing/mismatched for {comp["tag"]}')
    files = {}
    with tarfile.open(archive) as tf:
        members = tf.getmembers()
        # Extract only ordinary files; reject links, traversal and special entries.
        for member in members:
            parts = Path(member.name).parts[1:]
            if not parts:
                continue
            rel = Path(*parts)
            if rel.is_absolute() or '..' in rel.parts:
                raise ValueError('Unsafe source archive path')
            if member.isdir():
                continue
            source_member = member
            if member.issym() or member.islnk():
                linked = posixpath.normpath(posixpath.join(posixpath.dirname(member.name), member.linkname) if member.issym() else member.linkname)
                if linked.split('/')[0] != member.name.split('/')[0] or linked.startswith('/'):
                    raise ValueError('Unsafe source archive link')
                source_member = tf.getmember(linked)
            if not source_member.isfile():
                raise ValueError('Source archive contains special entries or chained links')
            data = tf.extractfile(source_member).read()
            fp = target / rel
            fp.parent.mkdir(parents=True, exist_ok=True)
            fp.write_bytes(data)
            fp.chmod(source_member.mode & 0o755)
            files[str(rel)] = hashlib.sha256(data).hexdigest()
    if not allow_new_digest and _digest(files) != comp.get('manifest_sha256'):
        raise ValueError('Extracted source fingerprint mismatch')
    manifest = {'repo': comp['repo'], 'commit': comp['commit'], 'archive_sha256': actual, 'files': files}
    (target / '.ascompare-source.json').write_text(json.dumps(manifest, indent=2) + '\n')
    if allow_new_digest:
        comp.update(archive_sha256=actual, manifest_sha256=_digest(files))
    _check_source(target, comp)
    return target


def prepare_sources(lock_entry, destination):
    """Return exact AA/AC source paths. Does not install or change environments."""
    result = {'id': lock_entry['id']}
    for name in ('aa', 'ac', 'bfb', 'suite'):
        comp = lock_entry.get(name)
        if comp:
            root = prepare_component(comp, destination)
            result[name + '_source'] = str(root)
            script = root / comp['entrypoint']
            if script.is_file():
                result[name + '_script'] = str(script)
            elif name != 'bfb':
                raise ValueError('Locked source entrypoint is missing: ' + str(script))
    return result


def prepare_panel(destination, versions=None, lock_path=None):
    """Prepare the selected locked Suite sources and write their runtime paths."""
    from .common import atomic_json

    lock = load_lock(lock_path)
    entries = {entry['id']: entry for entry in lock['versions'] if entry.get('suite')}
    requested = list(dict.fromkeys([lock['reference_version'], *(versions or entries)]))
    unknown = [version for version in requested if version not in entries]
    if unknown:
        raise ValueError('Suite version absent from lock: ' + ', '.join(unknown))
    root = Path(destination).resolve()
    rows = []
    for version in requested:
        print('Preparing ' + version, flush=True)
        rows.append(prepare_sources(entries[version], root))
    output = root / 'prepared_panel.json'
    atomic_json(output, {'reference_version': lock['reference_version'], 'sources': rows})
    return output


def _run(args, env):
    try:
        p = subprocess.run(args, env=env, capture_output=True, text=True, timeout=90)
        return {'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {'returncode': -1, 'stdout': '', 'stderr': str(exc)}


def verify_runtime(version, lock_entry, env=None):
    """Check locked source, executable version, help, deps, and BFB integrity.

    ready means launch checks passed, NOT a completed BAM/solver validation.
    Runtime input: id, aa_python, aa_script, ac_python, ac_script.
    """
    runtime_env = dict(os.environ if env is None else env)
    runtime_env.update(version.get('env', {}))
    extra_path = version.get('extra_path', '')
    if isinstance(extra_path, list):
        extra_path = os.pathsep.join(extra_path)
    if extra_path:
        runtime_env['PATH'] = extra_path + os.pathsep + runtime_env.get('PATH', '')
    result = {'id': version['id'], 'ready': True, 'diagnostics': [], 'validation_level': 'launch_only'}
    if version['id'] != lock_entry['id']:
        return {'id': version['id'], 'ready': False, 'diagnostics': ['Runtime/lock ID mismatch']}
    result['runtime_variants'] = version.get('dependency_overrides', [])
    for override in result['runtime_variants']:
        try:
            if sha256(override['path']) != override['sha256']:
                raise ValueError('Declared dependency binary fingerprint changed')
        except (KeyError, OSError, ValueError) as exc:
            result['ready'] = False
            result['diagnostics'].append(str(exc))
    for name in ('aa', 'ac'):
        info = {}
        try:
            comp = lock_entry[name]
            script = Path(version[name + '_script']).resolve()
            if not script.is_file():
                raise ValueError(f'No script: {script}')
            # Entrypoint parents are resolved from the known locked relative path.
            root = script
            for _ in Path(comp['entrypoint']).parts:
                root = root.parent
            if script != root / comp['entrypoint']:
                raise ValueError('Script does not match locked entrypoint')
            info.update(_check_source(root, comp))
            info['sha256'] = sha256(script)
            py = str(Path(version[name + '_python']).absolute())
            info['python'] = py
            component_env = dict(runtime_env)
            component_env['PATH'] = str(Path(py).parent) + os.pathsep + runtime_env.get('PATH', '')
            for arg, key in [('--version', 'version_probe'), ('--help', 'help_probe')]:
                info[key] = _run([py, str(script), arg], component_env)
                if info[key]['returncode']:
                    raise ValueError(f'{name.upper()} {arg} failed: {info[key]["stderr"][-1500:]}')
            output = info['version_probe']['stdout'] + info['version_probe']['stderr']
            if not re.search(r'(?<![\w.])v?' + re.escape(comp['version']) + r'(?![\w.])', output):
                raise ValueError(f'Expected {comp["version"]}, got version output {output.strip()}')
            info['version'] = comp['version']
            info['help'] = info['help_probe']['stdout'] + info['help_probe']['stderr']
            info['packages'] = _run([py, '-m', 'pip', 'freeze', '--all'], runtime_env)
        except (KeyError, ValueError, OSError) as exc:
            result['ready'] = False
            result['diagnostics'].append(f'{name}: {exc}')
        result[name] = info
    if result.get('aa', {}).get('version') and ('--solver' not in result['aa'].get('help', '') or version.get('solver')=='mosek'):
        code = '''import mosek
with mosek.Env() as e:
 with e.Task(0,0) as t:
  t.appendvars(1);t.putvarbound(0,mosek.boundkey.lo,0,0);t.putcj(0,1);t.optimize()
print('mosek optimization passed')'''
        probe = _run([version['aa_python'], '-c', code], runtime_env)
        result['aa_solver_probe'] = probe
        if probe['returncode']:
            result['ready'] = False
            result['diagnostics'].append('Historical AA requires working Mosek: ' + probe['stderr'][-1500:])
    if lock_entry.get('bfb'):
        # RECORD detects local edits to pip-distributed Python code (not a solver test).
        code = '''import importlib.metadata as m,json,hashlib,base64,pathlib
r={'version':None,'changed':[]}
try:
 d=m.distribution('bfbarchitect');r['version']=d.version
 import bfbarchitect
 from bfbarchitect._version import __version__
 r['module_version']=__version__;r['import_origin']=str(pathlib.Path(bfbarchitect.__file__).resolve())
 expected=pathlib.Path(d.locate_file('bfbarchitect/__init__.py')).resolve()
 if pathlib.Path(r['import_origin'])!=expected or __version__!=d.version:r['error']='BFB import origin/version differs from verified distribution'
 for f in d.files or []:
  if str(f).endswith('.py') and f.hash and f.hash.mode=='sha256':
   p=d.locate_file(f)
   actual=base64.urlsafe_b64encode(hashlib.sha256(p.read_bytes()).digest()).decode().rstrip('=')
   if actual!=f.hash.value:r['changed'].append(str(f))
except Exception as e:r['error']=str(e)
print(json.dumps(r))'''
        probe = _run([version.get('ac_python', ''), '-c', code], runtime_env)
        result['bfb_probe'] = probe
        try:
            bfb = json.loads(probe['stdout'])
            if bfb['version'] != lock_entry['bfb']['version'] or bfb.get('changed') or bfb.get('error'):
                raise ValueError(f'BFB dependency mismatch or edited install: {bfb}')
        except (ValueError, KeyError) as exc:
            result['ready'] = False
            result['diagnostics'].append(str(exc))
    if lock_entry.get('bfb') and result.get('ac', {}).get('version'):
        code = 'from bfbarchitect.BFBSolver import reconstruct_BFB_cbc; print(reconstruct_BFB_cbc([2],[0],[1],1,max_time=15,max_threads=1))'
        probe = _run([version['ac_python'], '-c', code], runtime_env)
        result['bfb_solver_probe'] = probe
        if probe['returncode']:
            result['ready'] = False
            result['diagnostics'].append('BFB CBC fallback cannot execute: ' + probe['stderr'][-1500:])
    return result
