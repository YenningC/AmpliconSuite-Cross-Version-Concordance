"""Bounded, per-sample CNVkit version controls using one fixed source AA seed filter."""
import hashlib
import json
import re
from pathlib import Path

from .concordance_cnv import load_builtin_map, run_local_upstream
from .common import digest
from .inputs import safe_id


def is_exact_version(value):
    return (isinstance(value, str) and len(value) <= 100
            and re.fullmatch(r'[0-9]+\.[0-9]+\.[0-9]+[A-Za-z0-9.+_-]*', value) is not None)


def resolve_registry(items, parent):
    """Resolve CNV interpreter paths without following a virtualenv's Python symlink."""
    if not isinstance(items, list):
        raise ValueError('cnvkit_runtimes must be an array')
    resolved, seen = [], set()
    for item in items:
        if not isinstance(item, dict) or set(item) - {'version', 'python', 'script', 'env', 'extra_path'}:
            raise ValueError('Invalid cnvkit_runtimes entry or unknown keys')
        if not all(isinstance(item.get(k), str) and item[k].strip() for k in ('version', 'python', 'script')):
            raise ValueError('Each CNVkit runtime needs version, python and script')
        version = item['version']
        if not is_exact_version(version):
            raise ValueError('cnvkit_runtimes requires exact CNVkit version IDs')
        if version in seen:
            raise ValueError(f'Duplicate CNVkit runtime for {version}')
        seen.add(version)
        entry = dict(item)
        for key in ('python', 'script'):
            path = parent / Path(item[key]).expanduser()
            entry[key] = str(path.absolute() if key == 'python' else path.resolve())
        if 'extra_path' in item:
            paths = item['extra_path']
            if not isinstance(paths, list) or not all(isinstance(v, str) and v for v in paths):
                raise ValueError('CNVkit extra_path must be a list of paths')
            entry['extra_path'] = [str((parent / Path(v).expanduser()).resolve()) for v in paths]
        if 'env' in item and (not isinstance(item['env'], dict)
                              or not all(isinstance(v, str) for v in item['env'].values())):
            raise ValueError('CNVkit env must be a mapping of names to strings')
        resolved.append(entry)
    return resolved


def comparison_panel(sample, config):
    """Map desired CNVkit versions to Suite rows; recipes alone establish only intent."""
    source_cnv = sample.get('source_cnvkit_version')
    panel = {source_cnv: []} if source_cnv else {}
    recipes = load_builtin_map()['versions']
    installed = {r['id']: r.get('cnvkit_version') for r in config.get('runtimes', [])}
    for suite in dict.fromkeys([sample['source_version'], *config.get('target_versions', [])]):
        version = installed.get(suite) or recipes.get(suite, {}).get('cnvkit_version')
        if version:
            panel.setdefault(version, []).append(suite)
    for version in config.get('local_cnv', {}).get('comparison_versions', []):
        panel.setdefault(version, [])
    return panel


def _cnv_environment(runtime):
    return dict(env=runtime.get('cnvkit_env', runtime.get('env', {})),
                extra_path=runtime.get('cnvkit_extra_path', runtime.get('extra_path', [])))


def _configured_tool(runtime):
    if not all(runtime.get(k) for k in ('cnvkit_python', 'cnvkit_script')):
        return None
    return dict(python=runtime['cnvkit_python'], script=runtime['cnvkit_script'],
                **_cnv_environment(runtime))


def _choose_tool(version, config, runtimes, source):
    native_source = _configured_tool(source) if source.get('cnvkit_version') == version else None
    choices = []
    for rt in runtimes.values():
        tool = _configured_tool(rt) if rt.get('cnvkit_version') == version else None
        if tool and tool not in choices:
            choices.append(tool)
    explicit = next((r for r in config.get('cnvkit_runtimes', []) if r['version'] == version), None)
    if explicit:
        matching = [t for t in choices if (t['python'], t['script']) == (explicit['python'], explicit['script'])]
        defaults = native_source if native_source in matching else (
            matching[0] if len(matching) == 1 else {'env': {}, 'extra_path': []})
        return {**defaults, **{k: v for k, v in explicit.items() if k != 'version'}}, None
    if native_source:
        return native_source, None
    if len(choices) == 1:
        return choices[0], None
    return None, 'multiple_configured_CNVkit_environments' if choices else 'exact_CNVkit_runtime_unavailable'


def _digest_value(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _locked_aa(probe, entry, runtime):
    aa = entry['aa']
    return (runtime and probe.get('version') == aa['version'] and probe.get('manifest_sha256')
            and probe['manifest_sha256'] == aa.get('manifest_sha256') and probe.get('commit') == aa['commit']
            and probe.get('sha256') and probe.get('source_root') and probe.get('python') == runtime.get('aa_python')
            and (Path(probe['source_root']) / aa['entrypoint']).resolve() == Path(runtime['aa_script']).resolve())


def _design(rows, version, source_runtime, aa_probe, sample, config):
    """Digest measured inputs, windows and filter; exclude CNV software being compared."""
    sample_id, suite = sample['sample_id'], sample['source_version']
    paths = {r.get('provenance') for r in rows}
    if len(paths) != 1 or not all(paths):
        raise ValueError('No unique measured CNV/filter provenance')
    provenance = json.loads(Path(paths.pop()).read_text())
    cache = provenance['cache']
    signature = cache['signature']
    actual = provenance['cnvkit']['version']
    if (actual != version or signature['cnvkit'] != provenance['cnvkit']
            or signature['sample_id'] != sample_id or provenance['version'] != suite
            or any(r.get('sample_id') != sample_id or r.get('version') != suite or r.get('cnvkit_version') != version
                   or r.get('cnv_cache_key') != cache['key'] for r in rows)
            or {r['seed_id'] for r in rows} != {s['seed_id'] for s in sample['seeds']}):
        raise ValueError('Inconsistent CNV/filter measurement identity')
    helper = provenance['amplified_intervals']
    expected_helper = Path(source_runtime['aa_script']).resolve().with_name('amplified_intervals.py')
    if Path(helper['path']).resolve() != expected_helper or helper['sha256'] != digest(expected_helper):
        raise ValueError('Source AA seed filter identity differs')
    args = helper['args']
    filter_identity = dict(suite=suite, aa_manifest=aa_probe['manifest_sha256'],
                           aa_python=aa_probe['python'], env=source_runtime.get('env', {}),
                           extra_path=source_runtime.get('extra_path', []), helper_sha256=helper['sha256'],
                           gain=args[args.index('--gain')+1], min_size=args[args.index('--cnsize_min')+1],
                           build=args[args.index('--ref')+1], no_cstats='--no_cstats' in args)
    return _digest_value(dict(sample=sample_id, inputs=signature['inputs'], io_contract=signature['io_contract'],
                              parameters=signature['parameters'], reference=config.get('_reference_identity'),
                              filter=filter_identity, matching_threshold=rows[0]['min_seed_overlap']))


def run_local_cnv_comparators(sample, config, entries, runtimes, verified, native_source_rows, out):
    """Return one local observation per CNVkit version and candidate, outside formal Suite decisions."""
    from .concordance import _new_attempt
    source = sample['source_version']
    sid, source_cnv = sample['sample_id'], sample.get('source_cnvkit_version')
    source_runtime = runtimes.get(source, {})
    aa_probe = verified.get(source, {}).get('aa', {})
    source_available = False
    panel = comparison_panel(sample, config)
    output = []
    for version, suites in panel.items():
        standard = version == source_cnv
        meta = dict(source_cnvkit_version=source_cnv, cnvkit_version=version, sample_id=sid, source_version=source,
                    source_filter_suite_version=source, source_filter_aa_version=entries[source]['aa']['version'],
                    comparison_design='CNVkit_version_with_fixed_source_AA_filter',
                    measurement_role='source_local_standard' if standard else 'comparison', target_suite_versions=suites,
                    measurement_origin='unavailable', local_design_id=None, actual_cnvkit_version=None)
        problem = sample.get('_input_failed') or sample.get('_source_cnvkit_issue')
        status, stage = ('failed', 'input_check') if sample.get('_input_failed') else ('unsupported', 'version_check')
        if not is_exact_version(version):
            problem = problem or 'invalid_target_CNVkit_version'
        if not source_cnv:
            problem = problem or 'source_CNVkit_version_unknown'
        if not problem and not _locked_aa(aa_probe, entries[source], source_runtime):
            problem = 'source_AA_filter_not_verified'
        tool, tool_issue = _choose_tool(version, config, runtimes, source_runtime)
        if not problem and tool_issue:
            problem = tool_issue
            status = 'missing' if not tool else status
            if tool_issue == 'multiple_configured_CNVkit_environments':
                status = 'unsupported'
        reuse = (native_source_rows and tool and version == source_runtime.get('cnvkit_version')
                 and tool == _configured_tool(source_runtime))
        if not problem and not standard and not source_available and not reuse:
            problem, stage = 'source_local_standard_unavailable', 'source_control'
        seeds = sample.get('seeds', [])
        if not problem and not seeds:
            problem, status, stage = 'No original candidate seeds', 'no_seeds', 'input_check'
        if problem:
            for seed in seeds or [{}]:
                output.append(dict(meta, seed=seed, seed_id=seed.get('seed_id', '__no_seeds__'),
                                   status=status, stage=stage, reason=str(problem), comparison_eligible=False))
            continue
        if reuse:
            rows = native_source_rows
            meta['measurement_origin'] = 'native_source_job'
        else:
            # Run the CNV tool in its own software environment. AA always receives the source environment.
            rt = dict(source_runtime, cnvkit_python=tool['python'], cnvkit_script=tool['script'],
                      cnvkit_version=version, cnvkit_env=tool['env'], cnvkit_extra_path=tool['extra_path'])
            job = _new_attempt(Path(out) / 'local_cnv_comparators' / safe_id(sid) / version)
            print(f'  {sid}: local CNVkit {version}, source AA filter {source}', flush=True)
            local_config = dict(config, _controlled_cnv_role=meta['measurement_role'],
                                _fixed_source_aa_identity={k: aa_probe.get(k) for k in
                                    ('version', 'commit', 'manifest_sha256', 'sha256', 'python', 'packages')})
            rows = run_local_upstream(sample, rt, entries[source], local_config, job)
            meta['measurement_origin'] = 'controlled_local_job'
        valid = [r for r in rows if r.get('comparison_eligible') is True]
        if valid:
            try:
                meta['local_design_id'] = _design(rows, version, source_runtime, aa_probe, sample, config)
                meta['actual_cnvkit_version'] = version
                if standard:
                    source_available = True
            except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
                rows = [dict(r, status='unsupported', stage='version_check', comparison_eligible=False,
                             reason=f'Local CNVkit design evidence missing or inconsistent: {exc}') for r in rows]
        for r in rows:
            observed = dict(r, **meta)
            if not standard:
                observed.update(calibration_role='candidate_observation', calibration_status='unverified',
                                calibration_reason='Measured on fixed source AA for a between-CNVkit comparison')
            output.append(observed)
    return output
