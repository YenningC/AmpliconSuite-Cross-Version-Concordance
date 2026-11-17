#!/usr/bin/env python3
"""Compose a ten-Suite server configuration; no download or tool/sample execution."""
import argparse
import copy
import json
from pathlib import Path
from ascompare.versions import load_lock

RUNTIME_KEYS=('aa_python','ac_python','suite_python','cnvkit_script','cnvkit_version','extra_path','env','dependency_overrides')


def build_config(prepared,templates,lock,reference,*,mosek_license=None,tool_bins=(),threads=2,fingerprint='full'):
    if prepared.get('reference_version')!=lock['reference_version']:
        raise ValueError('Prepared sources and lock have different fixed reference versions')
    if threads<1:
        raise ValueError('threads must be positive')
    if mosek_license and ('\n' in mosek_license or 'START_LICENSE' in mosek_license.upper()):
        raise ValueError('Provide only a license path/server address, never license contents')
    sources={x['id']:x for x in prepared['sources']}
    envs={x['id']:x for x in templates}
    refenv=envs.get(lock['reference_version'])
    ac200=envs.get('suite-1.6.1') or envs.get('suite-1.6.0')
    if not refenv or not ac200:
        raise ValueError('Need reference and AC2.0.0 runtime templates')
    versions=[]
    official=[e for e in lock['versions'] if e.get('suite')]
    official.sort(key=lambda e:(e['id']!=lock['reference_version'],tuple(-int(x) for x in e['suite']['version'].split('.'))))
    for entry in official:
        vid=entry['id'];source=sources.get(vid)
        if not source or not all(source.get(k) for k in ('suite_script','aa_script','ac_script')):
            raise ValueError('Missing prepared official Suite source: '+vid)
        template=ac200 if entry['ac']['version']=='2.0.0' else refenv
        runtime={k:copy.deepcopy(template[k]) for k in RUNTIME_KEYS if k in template}
        runtime.update(copy.deepcopy(source))
        for key in ('aa_python','ac_python','suite_python','cnvkit_script','cnvkit_version'):
            if not runtime.get(key):
                raise ValueError('Runtime template is missing '+key)
        if runtime['suite_python']!=runtime['ac_python']:
            raise ValueError('Suite and AC interpreters must match')
        if runtime['cnvkit_version']!='0.9.14':
            raise ValueError('This controlled panel requires explicit CNVkit0.9.14')
        if not isinstance(runtime.get('extra_path',[]),list):
            raise ValueError('Template extra_path must be a list')
        runtime['extra_path']=list(dict.fromkeys(list(tool_bins)+runtime.get('extra_path',[])))
        # Current locked historical AAs use MOSEK; modern releases expose Clarabel.
        runtime['solver']='clarabel' if entry['aa']['version'].startswith('1.6.') else 'mosek'
        if mosek_license:
            runtime.setdefault('env',{})['MOSEKLM_LICENSE_FILE']=mosek_license
        versions.append(runtime)
    return {'analysis_scope':'full_suite','reference_version':lock['reference_version'],
            'reference':reference,'versions':versions,'solver':'clarabel','downsample':-1,
            'threads':threads,'alignment_fingerprint':fingerprint,
            'full_suite':{'cngain':4.5,'cnsize_min':50000,'segmentation':'cbs'},
            'statistics':{'min_samples':20,'margins':{},'bootstrap':2000,'permutations':10000,'seed':20261117}}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--runtime-templates',type=Path)
    p.add_argument('--reference-python',type=Path,help='Alternative to templates: base runtime interpreter for reference and historical source scripts')
    p.add_argument('--ac200-python',type=Path,help='Alternative to templates: BFB1.0.1/AC2.0.0 runtime interpreter')
    p.add_argument('--cnvkit-script',type=Path,help='Optional common explicit CNVkit script with explicit interpreter paths')
    p.add_argument('--tool-bin',type=Path,action='append',default=[])
    p.add_argument('--mosek-license',help='Path or server address only; applied to all entries, also affecting BFB auto-detection')
    p.add_argument('--reference-build',required=True)
    p.add_argument('--fasta',type=Path,required=True)
    p.add_argument('--aa-data-repo',type=Path,required=True)
    p.add_argument('--threads',type=int,default=2)
    p.add_argument('--alignment-fingerprint',choices=['full','sampled'],default='full')
    p.add_argument('--lock',type=Path)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();lock=load_lock(a.lock)
    if a.out.exists():p.error('Refusing to overwrite an existing configuration')
    if a.runtime_templates:
        if a.reference_python or a.ac200_python or a.cnvkit_script:p.error('Use templates or explicit interpreters, not both')
        templates=json.loads(a.runtime_templates.read_text())
    else:
        if not a.reference_python or not a.ac200_python:p.error('Provide --runtime-templates or both interpreter paths')
        templates=[]
        for vid,py in [(lock['reference_version'],a.reference_python),('suite-1.6.1',a.ac200_python)]:
            # Preserve venv identity; resolving its python symlink would bypass it.
            py=py.absolute();ck=a.cnvkit_script.absolute() if a.cnvkit_script else py.parent/'cnvkit.py'
            templates.append({'id':vid,'suite_python':str(py),'aa_python':str(py),'ac_python':str(py),
                              'cnvkit_script':str(ck),'cnvkit_version':'0.9.14','extra_path':[]})
    reference={'build':a.reference_build,'fasta':str(a.fasta.absolute()),'aa_data_repo':str(a.aa_data_repo.absolute())}
    try:
        config=build_config(json.loads(a.prepared.read_text()),templates,lock,reference,mosek_license=a.mosek_license,
                            tool_bins=[str(x.absolute()) for x in a.tool_bin],threads=a.threads,fingerprint=a.alignment_fingerprint)
    except (KeyError,TypeError,ValueError) as exc:p.error(str(exc))
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(config,indent=2)+'\n')
    print(f'Wrote {len(config["versions"])} official Suite entries; fixed reference {config["reference_version"]}. No tools or patient data were executed. Run ascompare check next.')


if __name__=='__main__':main()
