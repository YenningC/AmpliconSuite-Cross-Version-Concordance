#!/usr/bin/env python3
"""Attach locked source/runtime paths to an existing-results manifest. No BAM I/O."""
import argparse
import copy
import json
from pathlib import Path
from ascompare.versions import load_lock
from ascompare.concordance_cnv import load_builtin_map


def build(manifest, prepared, templates, cnv_runtimes, lock):
    result = copy.deepcopy(manifest)
    entries = {e['id']: e for e in lock['versions'] if e.get('suite')}
    sources = {e['id']: e for e in prepared['sources']}
    template_map = {t['id']: t for t in templates}
    reference = template_map.get(lock['reference_version'])
    old_bfb = template_map.get('suite-1.6.1') or template_map.get('suite-1.6.0')
    cnv = {}
    for record in cnv_runtimes:
        version = record['version']
        if version in cnv:
            raise ValueError('Duplicate CNVkit runtime for ' + version)
        cnv[version] = record
    if result['mode'] == 'with_bam':
        # Keep exact input-version environments even if no Suite recipe uses them.
        registry = {row['version']: copy.deepcopy(row) for row in result.get('cnvkit_runtimes', [])}
        if len(registry) != len(result.get('cnvkit_runtimes', [])):
            raise ValueError('Duplicate CNVkit runtime in input manifest')
        for version, row in cnv.items():
            if version in registry and registry[version] != row:
                raise ValueError('Conflicting CNVkit runtime for ' + version)
            registry[version] = copy.deepcopy(row)
        result['cnvkit_runtimes'] = list(registry.values())
    recipes = load_builtin_map()['versions']
    if prepared['reference_version'] != lock['reference_version']:
        raise ValueError('Prepared source reference differs from the fixed lock')
    targets = result.get('target_versions', list(entries))
    requested = list(dict.fromkeys([lock['reference_version'], *targets,
                    *(s['source_version'] for s in result['samples'])]))
    runtimes = []
    for vid in requested:
        entry, source = entries[vid], sources[vid]
        template = template_map.get(vid) or (old_bfb if entry['ac']['version'] == '2.0.0' else reference)
        if not template:
            raise ValueError('Missing AA/AC runtime template for ' + vid)
        keys = ('aa_python', 'ac_python', 'extra_path', 'env', 'dependency_overrides')
        runtime = {k: copy.deepcopy(template[k]) for k in keys if k in template}
        runtime.update({k: source[k] for k in ('id', 'aa_script', 'ac_script')})
        runtime['solver'] = 'clarabel' if entry['aa']['version'].startswith('1.6.') else 'mosek'
        if result['mode'] == 'with_bam':
            ck = recipes[vid]['cnvkit_version']
            if ck not in registry:
                raise ValueError('Missing exact CNVkit runtime ' + ck + '; no latest-version fallback')
            cnv_runtime = registry[ck]
            runtime.update(cnvkit_python=cnv_runtime['python'], cnvkit_script=cnv_runtime['script'], cnvkit_version=ck)
            for key in ('env', 'extra_path'):
                if key in cnv_runtime:
                    runtime['cnvkit_' + key] = copy.deepcopy(cnv_runtime[key])
        runtimes.append(runtime)
    result.update(reference_version=lock['reference_version'], target_versions=requested, runtimes=runtimes)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True, help='Base manifest with samples and reference')
    p.add_argument('--prepared', type=Path, required=True)
    p.add_argument('--runtime-templates', type=Path, required=True)
    p.add_argument('--cnvkit-runtimes', type=Path, help='JSON list of exact version/python/script triples')
    p.add_argument('--lock', type=Path)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if a.out.exists():
        p.error('Refusing to overwrite existing manifest')
    if a.out.resolve().parent != a.input.resolve().parent:
        p.error('Keep output beside the input manifest so relative sample paths retain their meaning')
    try:
        result = build(json.loads(a.input.read_text()), json.loads(a.prepared.read_text()),
                       json.loads(a.runtime_templates.read_text()),
                       json.loads(a.cnvkit_runtimes.read_text()) if a.cnvkit_runtimes else [], load_lock(a.lock))
    except (KeyError, ValueError) as exc:
        p.error(str(exc))
    a.out.write_text(json.dumps(result, indent=2) + '\n')
    print('Prepared ' + str(len(result['runtimes'])) + ' locked runtime entries. No samples were read or executed.')


if __name__ == '__main__':
    main()
