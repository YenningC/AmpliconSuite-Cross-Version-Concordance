import argparse
import hashlib
import json
import os
import sys
import atexit
import shutil
from pathlib import Path
from .common import atomic_json, now, table
from .inputs import read_samples, safe_id, validate, validate_full


def main():
    parser=argparse.ArgumentParser(description='Local full-pipeline and fixed-seed AmpliconSuite version comparison')
    sub=parser.add_subparsers(dest='command',required=True)
    for command in ('check','run'):
        p=sub.add_parser(command); p.add_argument('--config',required=True); p.add_argument('--samples',required=True); p.add_argument('--out',required=True)
    p=sub.add_parser('demo'); p.add_argument('--out',required=True)
    p=sub.add_parser('report'); p.add_argument('--out',required=True)
    p=sub.add_parser('releases'); p.add_argument('--out',required=True)
    p=sub.add_parser('prepare'); p.add_argument('--version',required=True); p.add_argument('--out',required=True); p.add_argument('--lock')
    args=parser.parse_args()
    if args.command in ('check','run'):
        parser.error('Legacy BAM workflows are disabled. Use asbench check/run --manifest with the bounded 0.2 schema; whole-genome operations are prohibited.')
    from .versions import load_lock, verify_runtime, refresh_releases, prepare_sources
    if args.command=='releases':
        result=refresh_releases(args.out); print(json.dumps(result,ensure_ascii=False,default=str)); return 0
    if args.command=='prepare':
        lock=load_lock(args.lock)
        entry=next(x for x in lock['versions'] if x['id']==args.version)
        print(json.dumps(prepare_sources(entry,args.out),indent=2,default=str)); return 0
    if args.command=='demo':
        from .demo import demo
        demo(args.out); return 0
    out=Path(args.out).resolve(); out.mkdir(parents=True,exist_ok=True)
    import fcntl
    lock_handle=(out/'.run.lock').open('a+')
    try:
        fcntl.flock(lock_handle.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        raise ValueError('Another process is using this output directory')
    atexit.register(lock_handle.close)
    if args.command=='report':
        from .reporting import write_reports
        config=json.loads((out/'config.resolved.json').read_text())
        write_reports(out,json.loads((out/'records.json').read_text()),config['reference_version'],options={**config.get('statistics',{}),'versions':[x['id'] for x in config['versions']],'analysis_scope':config.get('analysis_scope','fixed_seed')},synthetic=config.get('synthetic',False)); return 0
    config_path=Path(args.config).resolve()
    config=json.loads(config_path.read_text())
    scope=config.get('analysis_scope','fixed_seed')
    if scope not in ('fixed_seed','full_suite'):
        raise ValueError('analysis_scope must be fixed_seed or full_suite')
    if config.get('synthetic'):
        raise ValueError('Synthetic fixtures run only through the separate demo command')
    for key in ('fasta','aa_data_repo'):
        config['reference'][key]=str((config_path.parent/Path(config['reference'][key]).expanduser()).resolve())
    lock=load_lock(str((config_path.parent/Path(config['version_lock'])).resolve()) if config.get('version_lock') else None)
    entries={x['id']:x for x in lock['versions']}
    ids=[safe_id(x['id']) for x in config['versions']]
    if len(ids)!=len(set(ids)) or config['reference_version'] not in ids:
        raise ValueError('Versions must be unique and include reference_version')
    if config['reference_version']!=lock['reference_version']:
        raise ValueError('reference_version must match the fixed release lock; use a separately reviewed lock to change it')
    for runtime in config['versions']:
        for key in (('aa_python','aa_script','ac_python','ac_script','suite_python','suite_script','cnvkit_script') if scope=='full_suite' else ('aa_python','aa_script','ac_python','ac_script')):
            runtime[key]=str((config_path.parent/Path(runtime[key]).expanduser()).absolute()) if key.endswith('_python') else str((config_path.parent/Path(runtime[key]).expanduser()).resolve())
    samples=read_samples(args.samples,full_mode=scope=='full_suite')
    # Prevent any output tree from replacing or containing originals/reference/source installation.
    input_paths=[Path(row[k]).resolve() for row in samples for k in (('alignment','index','normal_alignment','normal_index') if scope=='full_suite' else ('alignment','index','seeds_bed'))]
    input_paths += [Path(config['reference']['aa_data_repo']),Path(config['reference']['fasta'])]
    for p in input_paths:
        if out==p or out in p.parents or (p.is_dir() and p in out.parents):
            raise ValueError('Output directory must be separate from source inputs/resources')
    snapshot=out/'program_source'
    if not snapshot.exists():
        shutil.copytree(Path(__file__).parent,snapshot,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    free_bytes=shutil.disk_usage(out).free
    required_free=int(config.get('min_free_disk_gb',10)*1024**3)
    cram_bytes=sum(Path(row[k]).stat().st_size for row in samples for k in ('alignment','normal_alignment') if row.get(k,'').lower().endswith('.cram') and Path(row[k]).is_file())
    required_free+=int(config.get('cram_staging_space_factor',3)*cram_bytes)
    if free_bytes<required_free:
        raise ValueError('Insufficient free output disk space for configured reserve and CRAM staging estimate')
    atomic_json(out/'storage_check.json',{'free_bytes':free_bytes,'required_reserve_bytes':required_free,'note':'CRAM size factor is a planning estimate, not a guarantee of peak output space'})
    print('Checking inputs and provenance (large-file hashing policy is recorded)...',flush=True)
    checked,sources=(validate_full if scope=='full_suite' else validate)(config,samples)
    metadata={'created_utc':now(),'config':config,'lock':lock,'inputs':sources,'manifest':str(Path(args.samples).resolve()),'scope':scope}
    canonical={k:v for k,v in metadata.items() if k!='created_utc'}
    run_hash=hashlib.sha256(json.dumps(canonical,sort_keys=True).encode()).hexdigest()
    old=out/'provenance.json'
    if old.exists() and json.loads(old.read_text()).get('run_hash')!=run_hash:
        raise ValueError('Inputs/configuration changed; use a new output directory (no stale resume)')
    atomic_json(old,{**metadata,'run_hash':run_hash}); atomic_json(out/'config.resolved.json',config)
    atomic_json(out/'input_check.json',checked); table(out/'input_check.tsv',checked)
    atomic_json(out/'versions.lock.json',lock)
    env=os.environ.copy(); env['AA_DATA_REPO']=config['reference']['aa_data_repo']; env['PYTHONDONTWRITEBYTECODE']='1'
    runtimes={}
    for runtime in config['versions']:
        entry=entries[runtime['id']]
        verified=verify_runtime(runtime,entry,env)
        if scope=='full_suite':
            from .full_suite import verify_full_runtime
            verified=verify_full_runtime(runtime,entry,env,verified)
        verified['supported_labels']=entry.get('supported_labels',entry.get('supported_classes',['ecDNA','BFB','Linear','Complex-non-cyclic']))
        runtimes[runtime['id']]=verified
    runtime_identity={}
    for ident,r in runtimes.items():
        runtime_identity[ident]={'components':{k:{f:r.get(k,{}).get(f) for f in ('version','manifest_sha256','sha256')} for k in ('aa','ac','suite','cnvkit')},
          'packages':{k:r.get(k,{}).get('packages',{}).get('stdout') for k in ('aa','ac')},
          'variants':r.get('runtime_variants',[])}
    identity_path=out/'runtime_identity.json'
    if identity_path.exists():
        previous_identity=json.loads(identity_path.read_text())
        previous_checks=json.loads((out/'runtime_check.json').read_text()) if (out/'runtime_check.json').exists() else {}
        for ident,identity in previous_identity.items():
            # Previously blocked versions have no reusable scientific result; allow them to become runnable.
            if previous_checks.get(ident,{}).get('ready',True) and runtime_identity.get(ident)!=identity:
                raise ValueError('Runtime dependency/source identity changed; use a new output directory')
    atomic_json(identity_path,runtime_identity)
    atomic_json(out/'runtime_check.json',runtimes)
    table(out/'run_metadata.tsv',[{'date_utc':metadata['created_utc'],'program_version':'0.1.0','reference_version':config['reference_version'],'build':config['reference']['build'],'coordinate_system':'BED0 half-open','run_hash':run_hash,'source':x['path'],'sha256':x.get('sha256',''),'fingerprint_method':x.get('method','full SHA256'),'parameters':config} for x in sources])
    if args.command=='check':
        print('Preflight saved to '+str(out)); return 0 if all(s['status'] in ('ok','no_seeds') for s in checked) and all(r['ready'] for r in runtimes.values()) else 2
    from .runner import run_jobs
    from .reporting import write_reports
    if scope=='full_suite':
        from .full_suite import run_full_suite
        records=run_full_suite(config,checked,runtimes,out)
    else:
        records=run_jobs(config,checked,runtimes,out)
    atomic_json(out/'records.json',records)
    if config.get('official_comparisons'):
        from .official_compare import compare_completed_records
        compare_completed_records(records,out,config['official_comparisons'],config['reference_version'])
    write_reports(out,records,config['reference_version'],options={**config.get('statistics',{}),'versions':ids,'analysis_scope':scope},synthetic=False)
    print('Results: '+str(out),flush=True)
    return 0 if all(s['status'] in ('ok','no_seeds') for s in checked) and records and all(r['status'] in ('ok','not_detected','invalid','unsupported','no_seeds') for r in records) else 2

if __name__=='__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print('Interrupted; completed jobs preserved. Repeat the same command to resume.',file=sys.stderr); sys.exit(130)
    except Exception as exc:
        print(f'ERROR: {exc}',file=sys.stderr); sys.exit(2)
