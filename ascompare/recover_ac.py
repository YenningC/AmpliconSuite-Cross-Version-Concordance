"""AC-only recovery through the unchanged, locked Suite completed-AA route.

Original attempts are read-only. This is an explicitly composite full-pipeline
result, never a fabricated claim that the failed original Suite run succeeded.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil

from .common import atomic_json, digest, now
from .full_suite import (WHOLE_SAMPLE, _all_logs, _safe, audit_reference, build_command,
                         private_reference, read_cns_segment_count, read_seed_bed, runtime_env, verify_full_runtime)
from .parsing import aa_amplicon_count, parse_run
from .runner import check_ac_runtime, empty_record, execute
from .versions import _run, load_lock, verify_runtime


def _read(path):
    return json.loads(Path(path).read_text())


def _check_hashes(attempt, hashes):
    if not hashes:
        raise ValueError('Original failed result has no output checksum manifest')
    for name, expected in hashes.items():
        path=Path(name)
        if attempt not in path.resolve().parents:
            raise ValueError('Original output checksum points outside the failed attempt')
        if not path.is_file() or digest(path)!=expected:
            raise ValueError('Inherited output missing or changed: '+str(path))


def _required_hash(path, hashes):
    if str(path) not in hashes:
        raise ValueError('Required inherited evidence is not covered by original checksum: '+str(path))


def command_signature(argv):
    """Normalize only flag ordering; retain exact paths and parameter values."""
    switches={'--run_AA','--run_AC','--no_repo_check','--no_cstats'}
    valued={'-s','-o','-t','--bam','--normal_bam','--ref','--AA_src','--aa_python_interpreter',
            '--python3_path','--cnvkit_dir','--AA_runmode','--AA_extendmode','--downsample',
            '--cngain','--cnsize_min','--cnvkit_segmentation','--AA_solver'}
    if not isinstance(argv,list) or len(argv)<2 or not all(isinstance(a,str) for a in argv):
        raise ValueError('Invalid original execution command')
    options={};i=2
    while i<len(argv):
        flag=argv[i]
        if flag not in switches|valued or flag in options:
            raise ValueError('Unknown or duplicate inherited execution flag: '+flag)
        if flag in switches:
            options[flag]=True;i+=1
        else:
            if i+1>=len(argv) or argv[i+1] in switches|valued or (argv[i+1].startswith('-') and not re.fullmatch(r'-\d+(?:\.\d+)?',argv[i+1])):
                raise ValueError('Inherited execution flag has no valid value: '+flag)
            options[flag]=argv[i+1];i+=2
    return (tuple(argv[:2]),options)


def _verify_recovery_runtime(runtime, baseline, lock_entry, config):
    declared=lock_entry.get('supported_labels',lock_entry.get('supported_classes'))
    prior=baseline.get('supported_labels')
    if not isinstance(declared,list) or not declared or len(set(declared))!=len(declared):
        raise ValueError('Locked supported categories are absent or invalid; no inferred default allowed')
    if not isinstance(prior,list) or len(set(prior))!=len(prior) or set(prior)!=set(declared):
        raise ValueError('Baseline supported categories differ from the exact version lock')
    env=runtime_env(runtime)
    env['AA_DATA_REPO']=config['reference']['aa_data_repo']
    verified=verify_full_runtime(runtime,lock_entry,env,verify_runtime(runtime,lock_entry,env))
    if not verified.get('ready'):
        raise ValueError('Recovery runtime verification failed: '+'; '.join(verified.get('diagnostics',[])))
    # verify_full_runtime verifies executables; the registry is the authority
    # for category availability. Never let parse_run fall back to legacy labels.
    verified['supported_labels']=list(declared)
    for component in ('suite','aa','ac','cnvkit'):
        old,new=baseline.get(component,{}),verified.get(component,{})
        if not old.get('version'):
            raise ValueError('Original component version evidence missing: '+component)
        for key in ('version','manifest_sha256','sha256','package'):
            if key in old and old[key]!=new.get(key):
                raise ValueError('Recovery changed the verified component: '+component+'/'+key)
        if old.get('packages',{}).get('stdout')!=new.get('packages',{}).get('stdout'):
            raise ValueError('Recovery changed component Python dependencies: '+component)
    if lock_entry.get('bfb'):
        probe=_run([runtime['ac_python'],'-c',
            'import json;from bfbarchitect import detect_solver;print(json.dumps({"selected_backend":detect_solver(),"scope":"autodetection only; actual solve outcome retained in AC logs"}))'],env)
        verified['recovery_bfb_backend_probe']=probe
        if probe['returncode']:
            raise ValueError('BFB backend autodetection probe failed before AC recovery')
    return verified


def build_recovery_command(config, runtime, verified, out):
    out=Path(out);helptext=verified['suite'].get('help','')
    if not all(flag in helptext for flag in ('--completed_AA_runs','--run_AC','--no_repo_check')):
        raise ValueError('Recovery requires the verified native completed-AA route with explicit offline option')
    command=[runtime['suite_python'],runtime['suite_script'],'-s','sample','-o',str(out/'pipeline'),
             '-t',str(config.get('threads',2)),'--completed_AA_runs',str(out/'aa'),'--run_AC',
             '--ref',config['reference']['build'],'--AA_src',str(Path(runtime['aa_script']).parent),
             '--aa_python_interpreter',runtime['aa_python'],'--python3_path',runtime['ac_python'],'--no_repo_check']
    if '--completed_run_metadata' in helptext:
        command+=['--completed_run_metadata',str(out/'inherited/sample_run_metadata.json')]
    for value in command:_safe(value)
    return command


def recover_ac(attempt, config_path, runtime_check_path, input_check_path, out, extra_path=(), original_result=None):
    """Validate/copy inherited stages and rerun only native AC in a new directory.

    extra_path changes supporting executable lookup only (e.g. verified Bash5).
    The baseline config and original component identities remain mandatory.
    This function does not update the old result.json or the active job queue.
    """
    attempt=Path(attempt).resolve();out=Path(out).absolute()
    if out.resolve()==attempt or attempt in out.resolve().parents or out.resolve() in attempt.parents:
        raise ValueError('Recovery output must be separate from the original attempt')
    if out.exists():raise ValueError('Recovery requires a new output directory')
    _safe(out)
    config=_read(config_path);checks=_read(runtime_check_path);samples=_read(input_check_path)
    reference_source=Path(config['reference']['aa_data_repo']).resolve()
    if out.resolve()==reference_source or reference_source in out.resolve().parents or out.resolve() in reference_source.parents:
        raise ValueError('Recovery output must be separate from the reference resource directory')
    original_result=Path(original_result or attempt.parent/'result.json').resolve()
    old=_read(original_result)
    if old.get('status') not in ('failed','missing') or old.get('stage')!='ac' or old.get('endpoint')!='full_pipeline':
        raise ValueError('Only a recorded full-pipeline AC-stage failure can be recovered')
    if Path(old.get('run_dir','')).resolve()!=attempt:
        raise ValueError('Original failure record does not refer to this attempt')
    sample=next((s for s in samples if s['sample_id']==old['sample_id']),None)
    if not sample or sample.get('status')!='ok':raise ValueError('Original sample input check was not successful')
    runtime=next((dict(r) for r in config['versions'] if r['id']==old['version']),None)
    if runtime is None:raise ValueError('Version absent from original config')
    baseline=checks[old['version']]
    hashes=old.get('output_sha256',{});_check_hashes(attempt,hashes)
    execution=attempt/'execution.json';metadata_path=attempt/'pipeline/sample_run_metadata.json'
    cns=attempt/'pipeline/sample_cnvkit_output/tumor.cns';seeds_path=attempt/'pipeline/sample_AA_CNV_SEEDS.bed'
    for path in (execution,metadata_path,cns,seeds_path,attempt/'suite.log'):_required_hash(path,hashes)
    prior_execution=_read(execution)
    expected=build_command(config,runtime,baseline,attempt,str(attempt/'alignment/tumor.bam'),str(attempt/'alignment/normal.bam'))
    if command_signature(prior_execution.get('argv'))!=command_signature(expected) or prior_execution.get('endpoint')!='full_pipeline':
        raise ValueError('Original execution command differs from the supplied original config')
    metadata=_read(metadata_path)
    for key,component in (('AmpliconSuite-pipeline_version','suite'),('AA_version','aa'),('AC_version','ac'),('cnvkit_version','cnvkit')):
        version=baseline.get(component,{}).get('version')
        if not version or not re.search(r'(?<![\w.])'+re.escape(version)+r'(?![\w.])',str(metadata.get(key,''))):
            raise ValueError('Inherited native metadata version mismatch: '+key)
    if '--normal ' not in metadata.get('cnvkit_cmd','') or any(metadata.get(k) in (None,'','NA') for k in ('amplified_intervals_cmd','AA_cmd','AC_cmd')):
        raise ValueError('Inherited metadata does not establish matched-normal CNV, seed, AA and attempted AC stages')
    cns_segment_count=read_cns_segment_count(cns)
    seeds=read_seed_bed(seeds_path)
    if not seeds:raise ValueError('AC recovery requires nonempty inherited seeds and reconstruction')
    aa=attempt/'pipeline/sample_AA_results'
    essentials=[*aa.glob('*_summary.txt'),*aa.glob('*_amplicon*_graph.txt'),*aa.glob('*_amplicon*_cycles.txt')]
    if not essentials:raise ValueError('Inherited AA outputs are missing')
    for path in essentials:_required_hash(path,hashes)
    if any(p.is_symlink() for p in aa.rglob('*')):
        raise ValueError('Inherited AA directory contains links; refusing a writable recovery alias')
    # Verify the reference provenance read-only; audit_reference writes, so use it
    # only on the NEW recovery directory after the native command.
    reference_before=attempt/'reference_audit_before.json'
    _required_hash(reference_before,hashes)
    for entry in _read(reference_before):
        p=Path(entry['path']);st=p.stat()
        if st.st_size!=entry['size'] or st.st_mtime_ns!=entry['mtime_ns'] or (entry.get('sha256') and digest(p)!=entry['sha256']):
            raise ValueError('Reference resource changed since the inherited reconstruction')
    if not isinstance(runtime.get('extra_path',[]),list):raise ValueError('Runtime extra_path must be a list')
    runtime['extra_path']=[*map(str,extra_path),*runtime.get('extra_path',[])]
    lock=load_lock(Path(config_path).parent/'versions.lock.json' if (Path(config_path).parent/'versions.lock.json').is_file() else None)
    entry=next(v for v in lock['versions'] if v['id']==runtime['id'])
    verified=_verify_recovery_runtime(runtime,baseline,entry,config)
    command=build_recovery_command(config,runtime,verified,out)
    out.mkdir(parents=True);inherited=out/'inherited';inherited.mkdir()
    shutil.copytree(aa,out/'aa')
    for path in (metadata_path,cns,seeds_path,execution,original_result,Path(config_path),Path(runtime_check_path),Path(input_check_path)):
        target=inherited/('original_result.json' if path==original_result else path.name)
        shutil.copy2(path,target)
    stages={'input_check':{'status':'complete','source':'inherited'},'seed_processing':{'status':'complete','source':'inherited'},
            'aa':{'status':'complete','source':'inherited','evidence':str(aa)},'ac':{'status':'pending','source':'rerun'}}
    atomic_json(out/'stages.json',stages)
    count=aa_amplicon_count(out)
    pre=parse_run(out,sample['sample_id'],WHOLE_SAMPLE,runtime['id'],verified.get('supported_labels'),
                  reference_lengths=sample.get('reference_lengths'),allowed_missing_contigs=sample.get('missing_reference_contigs',[]))
    if not count or pre.get('status')!='missing' or pre.get('stage')!='ac':
        raise ValueError('Inherited AA evidence is not a complete nonempty reconstruction: '+pre.get('reason',''))
    aa_hashes={str(p.relative_to(out/'aa')):digest(p) for p in (out/'aa').rglob('*') if p.is_file()}
    manifest={'created_utc':now(),'status':'running','recovery_kind':'recovered_ac','endpoint':'full_pipeline',
              'interpretation':'Composite result: original CNVkit, seed selection and AA; native AC rerun. Original full Suite attempt remains failed.',
              'inherited_attempt':str(attempt),'inherited_result':str(original_result),'inherited_result_sha256':digest(original_result),
              'inherited_output_sha256':hashes,'inherited_stages':['input_check','cnvkit','seed_processing','aa'],
              'rerun_stages':['ac'],'argv':command,'support_extra_path':list(map(str,extra_path)),
              'baseline_runtime':baseline,'recovery_runtime':verified,'aa_copy_sha256':aa_hashes}
    manifest['inherited_cns_segment_count']=cns_segment_count
    atomic_json(out/'recovery_manifest.json',manifest)
    env=runtime_env(runtime);env.update(AA_DATA_REPO=str(private_reference(config['reference']['aa_data_repo'],out/'reference')),
        AA_SRC=str(Path(runtime['aa_script']).parent),AC_SRC=str(Path(runtime['ac_script']).parent),MPLCONFIGDIR=str(out/'mplcache'),
        OMP_NUM_THREADS=str(config.get('threads',2)),OPENBLAS_NUM_THREADS=str(config.get('threads',2)))
    try:
        rc=execute(command,out,env,out/'suite.log',timeout=config.get('job_timeout_seconds'))
        audit_reference(out);_check_hashes(attempt,hashes)
        if any(not (out/'aa'/p).is_file() or digest(out/'aa'/p)!=sha for p,sha in aa_hashes.items()):
            raise ValueError('Native recovery modified inherited AA evidence')
        log=_all_logs(out)
        if rc or re.search(r'Traceback \(most recent call last\)|declare: -A: invalid option|Failed to make AmpliconClassifier input',log):
            raise RuntimeError('Native AC recovery failed; inspect its separate suite.log')
        finish=out/'pipeline/sample_finish_flag.txt'
        if not finish.is_file() or finish.read_text().strip()!='All stages completed':
            raise RuntimeError('Native completed-AA route did not produce its own successful finish flag')
        ac=out/'pipeline/sample_classification'
        (out/'ac').symlink_to(ac.resolve(),target_is_directory=True)
        input_rows=[line.split() for line in (ac/'sample.input').read_text().splitlines() if line.strip()]
        expected_cycles={str(p.resolve()) for p in (out/'aa').glob('*_amplicon*_cycles.txt')}
        if len(input_rows)!=count or any(len(row)<3 for row in input_rows) or {str(Path(row[1]).resolve()) for row in input_rows}!=expected_cycles:
            raise ValueError('Native AC input does not cover every inherited AA amplicon exactly once')
        if any(Path(row[2]).resolve()!=Path(row[1]).resolve().with_name(Path(row[1]).name.replace('_cycles.txt','_graph.txt')) for row in input_rows):
            raise ValueError('Native AC input does not pair each cycle file with its inherited graph')
        check_ac_runtime(log,verified['ac']['version'])
        if verified['ac']['version'].startswith('2.'):
            ac_log=ac/'sample.log'
            if not ac_log.is_file() or not re.search(r' - INFO - Complete\s*$',ac_log.read_text(errors='replace')):
                raise RuntimeError('AC2 final completion log marker absent')
        stages['ac']={'status':'complete','returncode':rc,'source':'rerun'};atomic_json(out/'stages.json',stages)
        record=parse_run(out,sample['sample_id'],WHOLE_SAMPLE,runtime['id'],verified.get('supported_labels'),
                         reference_lengths=sample.get('reference_lengths'),allowed_missing_contigs=sample.get('missing_reference_contigs',[]))
    except Exception as exc:
        stages['ac']={'status':'failed','source':'rerun','reason':str(exc)};atomic_json(out/'stages.json',stages)
        record=empty_record(sample['sample_id'],WHOLE_SAMPLE,runtime['id'],'failed','ac',str(exc),supported_labels=verified.get('supported_labels',[]))
    record.update(endpoint='full_pipeline',analysis_scope='full_suite',recovery_kind='recovered_ac',run_dir=str(out),
                  seed_count=len(seeds),seed_intervals=seeds,seed_stage_evaluable=True,aa_executed=True,ac_executed=True,
                  aa_reconstruction_performed=True,detection_stage='ac',inherited_attempt=str(attempt),
                  recovery_manifest=str(out/'recovery_manifest.json'),inherited_stages=manifest['inherited_stages'],rerun_stages=['ac'],
                  reference_scope='human_only_incomplete_reference' if sample.get('missing_reference_contigs') else 'complete_reference',
                  cnvkit_version=baseline['cnvkit']['version'],suite_version=baseline['suite']['version'],
                  seed_thresholds=prior_execution.get('seed_thresholds',{}),runtime_variants=verified.get('runtime_variants',[]))
    manifest.update(status=record['status'],finished_utc=now(),reason=record.get('reason',''))
    atomic_json(out/'recovery_manifest.json',manifest)
    record.setdefault('evidence',[]).extend([str(out/'recovery_manifest.json'),str(attempt/'suite.log'),str(out/'suite.log')])
    record['output_sha256']={str(p):digest(p) for p in out.rglob('*') if p.is_file() and not p.is_symlink() and not any(x in p.relative_to(out).parts for x in ('reference','mplcache'))}
    atomic_json(out/'result.json',record)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('attempt','config','runtime-check','input-check','out'):parser.add_argument('--'+name,required=True)
    parser.add_argument('--original-result');parser.add_argument('--extra-path',action='append',default=[])
    args=parser.parse_args()
    result=recover_ac(args.attempt,args.config,args.runtime_check,args.input_check,args.out,args.extra_path,args.original_result)
    print(json.dumps({k:result.get(k) for k in ('status','recovery_kind','run_dir','reason')},ensure_ascii=False))
    return 0 if result['status'] in ('ok','not_detected','invalid','unsupported') else 2


if __name__=='__main__':raise SystemExit(main())
