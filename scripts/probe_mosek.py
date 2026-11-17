#!/usr/bin/env python3
"""Exercise each historical locked AA MOSEK interface using a synthetic flow.

No alignment data or license contents are read by this script. MOSEK itself
uses its configured license. Errors are recorded by type/code, never license text.
"""
import argparse
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
from ascompare.versions import load_lock, _check_source

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--sources',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
p.add_argument('--lock',type=Path)
a=p.parse_args()
a.sources=a.sources.resolve();a.out=a.out.resolve();a.out.parent.mkdir(parents=True,exist_ok=True)
os.chdir(a.out.parent);sys.dont_write_bytecode=True
import mosek
generic={}
try:
    with mosek.Env() as env:
        with env.Task(0,0) as task:
            task.appendvars(1);task.putvarbound(0,mosek.boundkey.lo,1,0);task.putcj(0,1)
            task.putobjsense(mosek.objsense.minimize);task.optimize()
            solution=[float(x) for x in task.getxx(mosek.soltype.itr)]
            generic={'passed':abs(solution[0]-1)<1e-6,'solution':solution}
except Exception as exc:
    generic={'passed':False,'error_type':type(exc).__name__,'error_code':str(getattr(exc,'errno','unknown'))}
results=[];by_commit={}
for entry in load_lock(a.lock)['versions']:
    if not entry.get('suite') or entry['id'] in ('suite-1.6.0','suite-1.6.1','suite-1.6.2'):
        continue
    commit=entry['aa']['commit'];record={'version':entry['id'],'AA':entry['aa']['version'],'commit':commit,'synthetic':True}
    if commit in by_commit:
        record.update(by_commit[commit]);record['shared_source_probe']=True;results.append(record);continue
    source=a.sources/('AmpliconArchitect-'+commit[:12])
    try:
        _check_source(source,entry['aa'])
        spec=importlib.util.spec_from_file_location('locked_mosek_'+commit[:12],source/'src/mosek_solver.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        # Upstream errors may print license diagnostics; do not expose those.
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            cn=[float(x) for x in module.call_mosek(1,1,[[0,1],[0,1]],[[1.,-1.],[1.,-1.]],[1.,1.],[-2.,-2.],[1.,1.],[0.,0.])]
        outcome={'CN':cn,'analytical_optimum':2,'tolerance':1e-3,'passed':max(abs(x-2) for x in cn)<1e-3,'api':'locked AA call_mosek'}
    except Exception as exc:
        outcome={'passed':False,'error_type':type(exc).__name__,'error_code':str(getattr(exc,'errno','unknown'))}
    by_commit[commit]=outcome;record.update(outcome);results.append(record)
result={'mosek_version':'.'.join(map(str,mosek.Env.getversion())),'generic_optimization':generic,'AA_interface_probes':results,'unique_AA_sources':len(by_commit)}
a.out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
if not generic.get('passed') or not all(x['passed'] for x in results):raise SystemExit(2)
