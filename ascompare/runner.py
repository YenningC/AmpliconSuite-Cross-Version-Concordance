"""Serial, resumable AA+AC jobs; all inputs remain local and read-only."""
import json
import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path
from .common import atomic_json, digest, now, table


def _terminate(proc):
    """Terminate the launched group even if its leader just exited."""
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


def execute(command, cwd, env, logfile, heartbeat=30, timeout=None):
    """Never invokes a shell; interrupts terminate the whole child process group."""
    with open(logfile,'w') as stream:
        stream.write(json.dumps({'argv':command,'start':now()})+'\n'); stream.flush()
        proc = subprocess.Popen(command,cwd=cwd,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            started = last = time.monotonic()
            while proc.poll() is None:
                time.sleep(0.2)
                current = time.monotonic()
                if timeout is not None and current - started > timeout:
                    raise TimeoutError(f'Command exceeded configured {timeout} seconds; see {logfile}')
                if current-last >= heartbeat:
                    print(f'  running {Path(logfile).name}; elapsed {int(current-started)}s',flush=True)
                    last = current
            return proc.returncode
        except BaseException:
            _terminate(proc)
            raise


def local_reference(source, target):
    """Private root prevents old AA coverage.stats writes into shared resources."""
    target.mkdir(parents=True,exist_ok=True)
    for child in Path(source).iterdir():
        if child.name == 'coverage.stats':
            continue
        dest = target/child.name
        if not dest.exists():
            dest.symlink_to(child.resolve(),target_is_directory=child.is_dir())
    return target


def prepare_alignment(row, config, staging):
    import pysam
    staging.mkdir(parents=True,exist_ok=True)
    source = Path(row['alignment'])
    dest = staging / 'input.bam'
    if source.suffix.lower() == '.cram':
        raise ValueError('Whole-file CRAM conversion is disabled in AmpliconSuite Benchmarks; supply an existing indexed BAM')
    elif source.suffix.lower() == '.bam':
        dest.symlink_to(source.resolve())
        suffix = '.csi' if row['index'].endswith('.csi') else '.bai'
        Path(str(dest)+suffix).symlink_to(Path(row['index']).resolve())
    else:
        raise ValueError('Alignment extension must be .bam or .cram')
    return str(dest)


def empty_record(sample,seed,version,status,stage,reason,**extra):
    return dict(sample_id=sample,seed_id=seed,version=version,status=status,stage=stage,reason=reason,
                labels=[],intervals=[],evidence=[],**extra)


def _help(verified,component):
    return str(verified.get(component,{}).get('help',''))


def check_ac_runtime(log, version):
    """Distinguish disabled integration from ordinary algorithmic candidate skips."""
    if not version.startswith('2.'):
        return
    errors=('skipping bfbarchitect integration','bfbarchitect failed for',
            'bfbarchitect is not installed','centromeres could not be loaded')
    if any(x in log.lower() for x in errors) or not re.search(r'BFBArchitect\s+\S+\s+integration enabled',log):
        raise RuntimeError('Required BFBArchitect integration was not confirmed or failed; inspect ac.log')


def run_jobs(config, samples, runtimes, out):
    if not config.get('synthetic'):
        raise ValueError('Legacy unguarded AA execution is disabled; use the bounded benchmark workflow')
    from .parsing import parse_run, aa_amplicon_count
    out = Path(out).resolve()
    out.mkdir(parents=True,exist_ok=True)
    records=[]
    total=sum(len(s['seeds']) for s in samples if s['status']=='ok')*len(config['versions'])
    done=0
    for sample in samples:
        if sample['status']=='no_seeds':
            for runtime in config['versions']:
                records.append(empty_record(sample['sample_id'],'__no_seeds__',runtime['id'],'no_seeds','seed_processing',sample.get('reason','No seeds supplied; AA/AC not evaluated'),supported_labels=runtimes[runtime['id']].get('supported_labels',[])))
            continue
        if sample['status']!='ok':
            for runtime in config['versions']:
                for seed in sample.get('seeds') or [{'seed_id':'__input_check__'}]:
                    records.append(empty_record(sample['sample_id'],seed['seed_id'],runtime['id'],
                                                'failed','input_check',sample.get('reason','Input validation failed')))
            continue
        for runtime in config['versions']:
            version=runtime['id']; verified=runtimes[version]
            supported=verified.get('supported_labels',runtime.get('supported_labels',['ecDNA','BFB','Complex-non-cyclic','Linear']))
            for seed in sample['seeds']:
                done+=1
                print(f'[{done}/{total}] {sample["sample_id"]} / {seed["seed_id"]} / {version}',flush=True)
                job=out/'jobs'/sample['sample_id']/seed['seed_id']/version
                job.mkdir(parents=True,exist_ok=True)
                saved=job/'result.json'
                if saved.exists():
                    try:
                        old=json.loads(saved.read_text())
                    except (OSError, ValueError):
                        old={}
                    hashes=old.get('output_sha256',{})
                    if old.get('status') in {'ok','not_detected','invalid','unsupported'} and old.get('run_dir') and hashes and all(Path(p).is_file() and digest(p)==h for p,h in hashes.items()):
                        print('  verified completed job: reuse',flush=True)
                        # Re-read verified raw artifacts, not editable cached classifications.
                        fresh=parse_run(old['run_dir'],sample['sample_id'],seed['seed_id'],version,supported,
                                        reference_lengths=sample.get('reference_lengths',config.get('_reference_lengths')),
                                        allowed_missing_contigs=sample.get('missing_reference_contigs',[]))
                        fresh.update(run_dir=old['run_dir'],output_sha256=hashes,seed=seed)
                        fresh['reference_scope']='human_only_incomplete_reference' if sample.get('missing_reference_contigs') else 'complete_reference'
                        atomic_json(saved,fresh)
                        records.append(fresh); continue
                if not verified.get('ready'):
                    record=empty_record(sample['sample_id'],seed['seed_id'],version,'failed','input_check','Runtime not verified',diagnostics=verified.get('diagnostics',[]),supported_labels=supported)
                    atomic_json(saved,record); records.append(record); continue
                previous=[int(p.name.split('_')[-1]) for p in job.glob('attempt_*') if p.is_dir() and p.name.split('_')[-1].isdigit()]
                attempt=job/f'attempt_{max(previous,default=0)+1:04d}'
                attempt.mkdir(); (attempt/'aa').mkdir(); (attempt/'ac').mkdir()
                stages={}; stage='seed_processing'
                record=None
                try:
                    bed=attempt/'seed.bed'
                    bed.write_text(f'{seed["chrom"]}\t{seed["start"]}\t{seed["end"]}\t{seed["seed_id"]}\n')
                    bam=prepare_alignment(sample,config,attempt/'alignment')
                    env=os.environ.copy()
                    runtime_env=runtime.get('env',{})
                    if not isinstance(runtime_env,dict) or any(not isinstance(k,str) or not isinstance(v,str) for k,v in runtime_env.items()):
                        raise ValueError('runtime.env must contain string environment values')
                    env.update(runtime_env)
                    env.update({'AA_DATA_REPO':str(local_reference(config['reference']['aa_data_repo'],attempt/'reference')),
                                'PYTHONDONTWRITEBYTECODE':'1','MPLCONFIGDIR':str(attempt/'mplcache'),
                                'OMP_NUM_THREADS':str(config.get('threads',1)), 'OPENBLAS_NUM_THREADS':str(config.get('threads',1))})
                    # Resolving a venv Python symlink changes it into the base interpreter
                    # and silently loses the selected site-packages and executable PATH.
                    bins=[str(Path(runtime[k]).expanduser().absolute().parent) for k in ('aa_python','ac_python')]
                    extra_path=runtime.get('extra_path',[])
                    if not isinstance(extra_path,list) or any(not isinstance(x,str) for x in extra_path):
                        raise ValueError('runtime.extra_path must be a list of executable directories')
                    extra_bins=[str(Path(x).expanduser().resolve()) for x in extra_path]
                    if any(not Path(x).is_dir() for x in extra_bins):
                        raise ValueError('runtime.extra_path contains a missing directory')
                    env['PATH']=os.pathsep.join(bins+extra_bins+[env.get('PATH','')])
                    ac_env=env.copy()
                    ac_env['PATH']=os.pathsep.join([bins[1],bins[0]]+extra_bins+[runtime_env.get('PATH',os.environ.get('PATH',''))])
                    atomic_json(attempt/'seed_provenance.json',{'seed':seed,'coordinate_system':'BED0 half-open','meaning':'candidate region, not confirmed amplification','scope':'seed-anchored full connected reconstruction','source':sample['seeds_bed']})
                    stages[stage]={'status':'complete'}
                    aa=[runtime['aa_python'],runtime['aa_script'],'--bam',bam,'--bed',str(bed),'--out',str(attempt/'aa'/'reconstruction'),
                        '--ref',config['reference']['build'],'--runmode','FULL','--extendmode','EXPLORE','--downsample',str(config.get('downsample',-1))]
                    if '--no_cstats' in _help(verified,'aa'):
                        aa+=['--no_cstats']
                    solver=runtime.get('solver',config.get('solver','mosek'))
                    if '--solver' in _help(verified,'aa'):
                        aa+=['--solver',solver]
                    elif solver!='mosek':
                        raise ValueError('Historical AA does not support requested solver; no fallback allowed')
                    stage='aa'; stages[stage]={'status':'running','argv':aa,'started':now()}; atomic_json(attempt/'stages.json',stages)
                    rc=execute(aa,attempt,env,attempt/'aa.log',timeout=config.get('job_timeout_seconds'))
                    stages[stage].update(status='complete' if rc==0 else 'failed',returncode=rc,ended=now()); atomic_json(attempt/'stages.json',stages)
                    if rc:
                        raise RuntimeError(f'AA exited {rc}; see aa.log')
                    cycles=sorted((attempt/'aa').glob('*_cycles.txt'))
                    amp_count=aa_amplicon_count(attempt)
                    graphs=sorted((attempt/'aa').glob('*_graph.txt'))
                    if amp_count is None or amp_count != len(cycles) or amp_count != len(graphs):
                        stage='result_reading'
                        raise FileNotFoundError('AA summary missing/unrecognized or graph/cycles count mismatch; AC not started')
                    if cycles:
                        rows=[]
                        for cycle in cycles:
                            graph=Path(str(cycle).replace('_cycles.txt','_graph.txt'))
                            if not graph.exists() or not graph.stat().st_size or not cycle.stat().st_size:
                                stage='result_reading'; raise FileNotFoundError('AA graph/cycles pair missing or empty')
                            # Official AC input: sample, cycles, graph (sample may include amplicon prefix).
                            rows.append('\t'.join(['reconstruction',str(cycle),str(graph)]))
                        inputfile=attempt/'ac.input'; inputfile.write_text('\n'.join(rows)+'\n')
                        ac=[runtime['ac_python'],runtime['ac_script'],'-i',str(inputfile),'--ref',config['reference']['build'],'-o',str(attempt/'ac'/'classification'),
                            '--min_size',str(config.get('ac_min_size',50000)),'--min_flow',str(config.get('ac_min_flow',1.0))]
                        stage='ac'; stages[stage]={'status':'running','argv':ac,'started':now()}; atomic_json(attempt/'stages.json',stages)
                        rc=execute(ac,attempt,ac_env,attempt/'ac.log',timeout=config.get('job_timeout_seconds'))
                        stages[stage].update(status='complete' if rc==0 else 'failed',returncode=rc,ended=now()); atomic_json(attempt/'stages.json',stages)
                        if rc:
                            raise RuntimeError(f'AC exited {rc}; see ac.log')
                        log=(attempt/'ac.log').read_text(errors='replace')
                        # A successful exit must not hide a disabled required new classifier.
                        check_ac_runtime(log, verified.get('ac',{}).get('version',''))
                    else:
                        stages['ac']={'status':'skipped_zero_amplicons'}; atomic_json(attempt/'stages.json',stages)
                    stage='result_reading'
                    record=parse_run(attempt,sample['sample_id'],seed['seed_id'],version,supported,
                                     reference_lengths=sample.get('reference_lengths',config.get('_reference_lengths')),
                                     allowed_missing_contigs=sample.get('missing_reference_contigs',[]))
                except KeyboardInterrupt:
                    stages[stage]={'status':'interrupted','ended':now()}; atomic_json(attempt/'stages.json',stages)
                    atomic_json(saved,empty_record(sample['sample_id'],seed['seed_id'],version,'failed',stage,'Interrupted; resume creates new attempt',supported_labels=supported))
                    raise
                except Exception as exc:
                    stages[stage]={**stages.get(stage,{}),'status':'failed','reason':str(exc)}; atomic_json(attempt/'stages.json',stages)
                    record=empty_record(sample['sample_id'],seed['seed_id'],version,'missing' if isinstance(exc,FileNotFoundError) else 'failed',stage,str(exc),supported_labels=supported)
                record['run_dir']=str(attempt)
                record['supported_labels']=supported
                record['seed']=seed
                record['reference_scope']='human_only_incomplete_reference' if sample.get('missing_reference_contigs') else 'complete_reference'
                # Check output integrity on resume; never reuse a truncated successful output.
                paths=[p for p in attempt.rglob('*') if p.is_file() and not p.is_symlink() and 'alignment' not in p.parts and 'reference' not in p.parts and 'mplcache' not in p.parts]
                record['output_sha256']={str(p):digest(p) for p in paths}
                record.setdefault('evidence',[]).append(str(attempt/'stages.json'))
                atomic_json(saved,record); records.append(record)
                atomic_json(out/'records.json',records)
    table(out/'anomalies.tsv',[{'sample_id':r['sample_id'],'seed_id':r['seed_id'],'version':r['version'],'status':r['status'],'stage':r.get('stage',''),'reason':r.get('reason',''),'evidence':r.get('evidence',[]),'diagnostics':r.get('diagnostics',[])} for r in records if r['status'] not in ('ok','not_detected')])
    return records
