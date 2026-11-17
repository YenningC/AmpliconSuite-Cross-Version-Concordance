#!/usr/bin/env python3
"""Package this GitHub source tree, including public docs, hidden config and previews.

Run from an installed comparison environment. Output must be outside the source
repository. Uses a code/docs/examples allowlist; local data is never collected.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import zipfile


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--wheel',action='store_true',help='Build wheel with installed setuptools; no dependencies bundled')
    args=parser.parse_args()
    source=args.source.resolve();out=args.out.resolve()
    if source==out or source in out.parents:
        parser.error('Output must be outside the repository')
    required=('README.md','README.zh-CN.md','GITHUB_UPLOAD.md','LICENSE_STATUS.md','.gitignore',
              'pyproject.toml','requirements.txt','versions.lock.json','official_tools.lock.json','validation_summary.json')
    files=[]
    for name in required:
        p=source/name
        if not p.is_file() or p.is_symlink():parser.error('Missing regular public file: '+name)
        files.append(p)
    allowed={'.py','.json','.md','.tsv','.bed','.txt','.sh','.yml','.yaml'}
    for folder in ('ascompare','tests','examples','docs','scripts','.github'):
        for p in (source/folder).rglob('*'):
            if '__pycache__' in p.parts:continue
            if p.is_symlink():parser.error('Refusing symlink: '+str(p.relative_to(source)))
            if p.is_file() and (p.suffix in allowed or p.suffix=='.png' and p.parent==source/'docs/assets'):
                rel=p.relative_to(source).as_posix()
                if any(s in rel.lower() for s in ('runtimes.local','runtime_probes','archive_v0.1','memory.md')):
                    parser.error('Internal/local file in public allowlist: '+rel)
                files.append(p)
    for name in ('LICENSE','repository_validation.json'):
        p=source/name
        if p.is_file():
            if p.is_symlink():parser.error('Refusing symlink: '+name)
            files.append(p)
    version=re.search(r'__version__\s*=\s*[\'\"]([^\'\"]+)',(source/'ascompare/__init__.py').read_text()).group(1)
    prefix='ampliconsuite-benchmarks-'+version+'-github'
    archive_root='AmpliconSuite-Benchmarks'
    out.mkdir(parents=True,exist_ok=True)
    names=[prefix+'.zip',prefix+'.tar.gz','release_manifest.json']
    if any((out/name).exists() for name in names):parser.error('Use a new output directory; release files are immutable')
    files=sorted(set(files));manifest={}
    with zipfile.ZipFile(out/names[0],'w',compression=zipfile.ZIP_DEFLATED) as z,tarfile.open(out/names[1],'w:gz') as t:
        for p in files:
            rel=p.relative_to(source).as_posix();data=p.read_bytes()
            manifest[rel]=hashlib.sha256(data).hexdigest()
            z.write(p,archive_root+'/'+rel);t.add(p,arcname=archive_root+'/'+rel,recursive=False)
    if args.wheel:
        subprocess.run([sys.executable,'-m','pip','wheel','--no-deps','--no-build-isolation',
                        '--wheel-dir',str(out),str(source)],check=True)
    artifacts={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir())
               if p.is_file() and p.suffix in ('.whl','.zip','.gz')}
    result=dict(program='AmpliconSuite Benchmarks',version=version,profile='public GitHub source',
                archive_root=archive_root,files=manifest,artifacts=artifacts,
                scope='code, public documentation and synthetic examples; no patient data, private configuration or environments')
    (out/'release_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'SHA256SUMS').write_text(''.join(value+'  '+name+'\n' for name,value in artifacts.items()))
    print(json.dumps(dict(files=len(files),artifacts=artifacts),indent=2))


if __name__=='__main__':
    main()
