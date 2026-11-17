#!/usr/bin/env python3
"""Create isolated modern runtime environments from an existing Conda/base Python.

Requires common science dependencies, samtools, R and native CBC in the base
installation. Downloads only locked public source and public Python packages.
Never changes that base installation and never reads patient files.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from ascompare.versions import load_lock,prepare_sources
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--base-python',type=Path,required=True)
p.add_argument('--work',type=Path,required=True)
p.add_argument('--lock',type=Path)
a=p.parse_args();a.work=a.work.resolve();a.work.mkdir(parents=True,exist_ok=True)
base=str(a.base_python.absolute())
probe=json.loads(subprocess.check_output([base,'-c',"import sys,json; print(json.dumps({'prefix':sys.prefix,'base':sys.base_prefix}))"],text=True))
if probe['prefix']!=probe['base']:
    p.error('Use a Conda/base Python with scientific dependencies; a parent venv does not pass its local packages via --system-site-packages')
lock=load_lock(a.lock);entries={x['id']:x for x in lock['versions']};records=[]
for version,name,pulp_version in [('suite-1.6.2','reference-env','4.2.0'),('suite-1.6.1','ac200-env','3.3.0')]:
    target=a.work/name
    if target.exists():p.error('Refusing to overwrite existing runtime environment: '+str(target))
    source=prepare_sources(entries[version],a.work/'sources')
    subprocess.run([base,'-m','venv','--system-site-packages',str(target)],check=True)
    py=str(target/'bin/python');build=a.work/(name+'-build');build.mkdir()
    packages=[]
    for component in ('bfb','ac'):
        dest=build/component;shutil.copytree(source[component+'_source'],dest);packages.append(str(dest))
    command=[py,'-m','pip','install','--ignore-installed','pulp=='+pulp_version,'mosek==11.2.8','gurobipy==14.0.0','CNVkit==0.9.14']+packages
    with (a.work/(name+'.install.log')).open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    frozen=subprocess.check_output([py,'-m','pip','freeze','--all'],text=True)
    (a.work/(name+'.requirements.frozen.txt')).write_text(frozen)
    records.append(dict(source,suite_python=py,aa_python=py,ac_python=py,cnvkit_script=str(target/'bin/cnvkit.py'),cnvkit_version='0.9.14',extra_path=[str(a.base_python.absolute().parent)]))
(a.work/'runtime_templates.json').write_text(json.dumps(records,indent=2)+'\n')
print('Prepared runtime_templates.json; run the documented native CBC adaptation for historical ARM environments if required, then ascompare check.')
