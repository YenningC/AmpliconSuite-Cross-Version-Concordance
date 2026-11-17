#!/usr/bin/env python3
"""Audit the active CNVkit package's source vs installed RECORD; no patient data."""
import argparse
import base64
import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import cnvlib

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
dist=metadata.distribution('CNVkit');root=Path(cnvlib.__file__).resolve().parent
inventory={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in root.rglob('*.py')}
records={str(Path(dist.locate_file(f)).resolve()):f for f in (dist.files or [])}
changed=[];missing=[];checked=0
for rel,sha in inventory.items():
    f=records.get(str(root/rel))
    if f is None or not f.hash or f.hash.mode!='sha256':
        missing.append(rel);continue
    expected=f.hash.value;actual=base64.urlsafe_b64encode(bytes.fromhex(sha)).decode().rstrip('=')
    checked+=1
    if expected!=actual:changed.append(rel)
status='mismatch' if changed else ('verified_against_distribution_RECORD' if not missing else 'partially_verified_RECORD_missing')
cbs=root/'segmentation/cbs.py'
result={'CNVkit':dist.version,'import_origin':str(root),'source_inventory_sha256':hashlib.sha256(json.dumps(inventory,sort_keys=True,separators=(',',':')).encode()).hexdigest(),'source_sha256':inventory,'record_status':status,'record_files_checked':checked,'record_unverifiable_files':missing,'changed_files':changed,'CBS_seed_evidence':[{'line':i,'text':line.strip()} for i,line in enumerate(cbs.read_text().splitlines(),1) if 'set.seed(' in line]}
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))
