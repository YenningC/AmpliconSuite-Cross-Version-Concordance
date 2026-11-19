"""Component-specific identities for whole-stage AA/AC output reuse.

Suite recipe identifiers and source checkout locations are provenance, not
algorithm identities.  Unknown environment/dependency differences conservatively
cause misses.  Native AA partial caches are never passed across commits.
"""
import os
import hashlib
import json
from pathlib import Path
import platform
import re
import uuid

from .common import digest

ENGINE = 'bounded-component-stages-v1'
# None of the locked AA entrypoints/modules import these distributions.  Their
# transitive dependencies are deliberately NOT removed: AA may use those too.
AA_UNUSED_DISTRIBUTIONS = {'bfbarchitect', 'pulp', 'ampliconclassifier', 'cnvkit',
                           'ampliconsuite-pipeline'}


def _packages(info, component):
    probe = info.get('packages')
    if (not isinstance(probe, dict) or probe.get('returncode') != 0
            or not isinstance(probe.get('stdout'), str)):
        # Launch checks can succeed even when pip cannot inventory the runtime.
        # Keep ordinary execution available, but an unknown inventory must never
        # stand in for a verified shared dependency environment.  A fresh nonce
        # prevents both whole-job resume and component reuse for this invocation.
        return dict(status='unverified', reuse_disabled_reason='Package inventory unavailable',
                    execution_nonce=uuid.uuid4().hex)
    text = probe['stdout']
    rows = []
    for line in text.splitlines():
        name = re.split(r'==|===|\s*@\s*|\[', line.strip(), maxsplit=1)[0]
        name = re.sub(r'[-_.]+', '-', name).lower()
        if component == 'aa' and name in AA_UNUSED_DISTRIBUTIONS:
            continue
        if line.strip():
            rows.append(line.strip())
    # pip freeze may contain a private URL with credentials.  Persist only its
    # one-way digest, never the original requirements or arbitrary environment.
    return dict(returncode=probe.get('returncode'), count=len(rows),
                freeze_sha256=hashlib.sha256(json.dumps(sorted(rows)).encode()).hexdigest())


def _environment_identity(environment, runtime):
    material = {'PATH', 'HOME', 'LANG', 'TZ', 'VIRTUAL_ENV', 'CONDA_PREFIX',
                'AA_DATA_REPO', 'TMPDIR', 'TEMP', 'TMP', 'MPLBACKEND'}
    prefixes = ('PYTHON', 'MOSEK', 'GRB_', 'OMP_', 'OPENBLAS_', 'MKL_',
                'NUMEXPR_', 'VECLIB_', 'LD_', 'DYLD_', 'LC_')
    keys = material | set(runtime.get('env', {}))
    return {k: hashlib.sha256(v.encode()).hexdigest() for k, v in sorted(environment.items())
            if k != 'MPLCONFIGDIR' and (k in keys or k.startswith(prefixes))}


def _file_identity(path):
    p = Path(path).expanduser()
    return dict(path=str(p.absolute()), resolved=str(p.resolve()),
                sha256=digest(p) if p.is_file() else None)


def component_identity(component, runtime, entry, verified, environment):
    info = verified.get(component, {})
    comp = entry[component]
    # The actual script is hashed again even if runtime verification was earlier.
    # Source inventory validation is performed by verify_runtime before this call.
    result = dict(engine=ENGINE, component=component,
                  source={k: comp.get(k) for k in ('repo', 'commit', 'version', 'manifest_sha256')},
                  verified={k: info.get(k) for k in ('version', 'commit', 'manifest_sha256', 'sha256')},
                  entrypoint_sha256=digest(runtime[component + '_script']),
                  interpreter=_file_identity(runtime[component + '_python']),
                  packages=_packages(info, component),
                  platform=dict(system=platform.system(), machine=platform.machine()),
                  # Private output locations cannot influence scientific identity.
                  environment_sha256=_environment_identity(environment, runtime),
                  dependency_overrides=[])
    for item in verified.get('runtime_variants', runtime.get('dependency_overrides', [])):
        if item.get('component') and item['component'].lower() != component:
            continue
        result['dependency_overrides'].append({**{k: item[k] for k in
            ('component', 'path', 'sha256') if k in item},
            'actual_sha256': digest(item['path']) if item.get('path') and Path(item['path']).is_file() else None})
    if component == 'aa':
        result['solver'] = runtime.get('solver', 'clarabel' if comp['version'].startswith('1.6.') else 'mosek')
        paths = [Path.home() / 'mosek' / 'mosek.lic']
        for name in ('MOSEKLM_LICENSE_FILE', 'GRB_LICENSE_FILE'):
            value = environment.get(name)
            if value and '@' not in value:
                paths.append(Path(value).expanduser())
        result['license_file_metadata_only'] = [dict(path=str(p), exists=p.is_file(),
            bytes=p.stat().st_size if p.is_file() else None,
            mtime_ns=p.stat().st_mtime_ns if p.is_file() else None) for p in paths]
    else:
        result['bfb'] = verified.get('bfb_probe', {}).get('stdout')
        result['bfb_lock'] = entry.get('bfb')
    return result


def effective_environment(component, runtime, reference):
    env = dict(os.environ, **runtime.get('env', {}))
    env.update(AA_DATA_REPO=str(Path(reference['aa_data_repo']).resolve()),
               PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    # Preserve 0.4.2's actual binary lookup order.  A genuinely different AC bin
    # in AA's PATH conservatively causes a miss, rather than changing execution.
    order = ('aa', 'ac') if component == 'aa' else ('ac', 'aa')
    bins = [str(Path(runtime[name + '_python']).absolute().parent) for name in order]
    inherited_path = env.get('PATH', '') if component == 'aa' else os.environ.get('PATH', '')
    env['PATH'] = os.pathsep.join(bins + runtime.get('extra_path', []) + [inherited_path])
    if component == 'aa':
        env['PYTHONHASHSEED'] = env.get('PYTHONHASHSEED', '0')
    return env
