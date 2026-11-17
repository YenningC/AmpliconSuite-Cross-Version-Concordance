#!/usr/bin/env python3
"""Fetch only locked public sources. No patient files are used."""
import argparse
import json
from pathlib import Path
from ascompare.versions import load_lock, prepare_sources
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--out',type=Path,required=True)
p.add_argument('--lock',type=Path)
a=p.parse_args();lock=load_lock(a.lock);results=[]
for entry in lock['versions']:
    print('Preparing '+entry['id'],flush=True)
    results.append(prepare_sources(entry,a.out))
(a.out/'prepared_panel.json').write_text(json.dumps({'reference_version':lock['reference_version'],'sources':results},indent=2)+'\n')
