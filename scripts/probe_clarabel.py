#!/usr/bin/env python3
"""Actual locked AA solver smoke test without patient data; expected CN=2."""
import argparse
import importlib.util
import json
from pathlib import Path
from ascompare.versions import load_lock, _check_source
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--sources',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
p.add_argument('--lock',type=Path)
a=p.parse_args();results=[]
for entry in load_lock(a.lock)['versions']:
    if entry['id'] not in ('suite-1.6.0','suite-1.6.1','suite-1.6.2'):continue
    source=a.sources/('AmpliconArchitect-'+entry['aa']['commit'][:12]);_check_source(source,entry['aa'])
    spec=importlib.util.spec_from_file_location('clarabel_solver',source/'src/clarabel_solver.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    cn=[float(x) for x in module.call_clarabel(1,1,[[0,1],[0,1]],[[1,-1],[1,-1]],[1.,1.],[-2.,-2.])]
    results.append({'version':entry['id'],'CN':cn,'analytical_optimum':2,'passed':max(abs(x-2) for x in cn)<1e-3,'synthetic':True})
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
if not results or not all(r['passed'] for r in results):raise SystemExit(2)
