"""Run the exact official Suite wrapper: CNVkit -> seeds -> AA -> AC.

The full-pipeline endpoint is per sample. It must never be pooled with fixed-seed
AA/AC endpoints. No external CNV/seed file, upload, reference fetch, or reuse of
historical scientific outputs is supplied to the wrapper.
"""
from __future__ import annotations
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
from .common import atomic_json, digest, now, table
from .runner import execute, local_reference, empty_record, check_ac_runtime
from .versions import _check_source, _run

SAFE = re.compile(r'^[A-Za-z0-9_./:+-]+$')
WHOLE_SAMPLE = '__whole_sample__'
# Exact audited source commits: ordinary BAM flow only reads local last_updated.txt;
# network downloads/uploads are confined to options this adapter never supplies.
OFFLINE_LEGACY_COMMITS = {
    '543b7ebc8b0d251eefbbcd14e5c1e27635e09a6a': 'suite-1.2.2',
    '306ec9b87ac169c5f9bacc034f4f220a6cf179c0': 'suite-1.3.3',
    '8e71026d799b5ca82604e53838162818b076386d': 'suite-1.3.4',
    '9f3629430829be260072053bc7de1a3e0a0a1093': 'suite-1.3.7',
    '3817f61043d526fdd57e20db1b26fd1ebbe0e664': 'suite-1.4.0',
    'a6457444f4889ca640c4fd573244524b442832b2': 'suite-1.5.0',
    '3dfe637eb3de1c522f2b209cbe158034c72be32a': 'suite-1.5.2',
}


def offline_policy(commit, helptext):
    if '--no_repo_check' in helptext:
        return 'explicit_no_repo_check'
    if commit in OFFLINE_LEGACY_COMMITS:
        return 'audited_intrinsically_offline_BAM_flow'
    raise ValueError('No audited offline BAM path for this exact Suite source')


def _safe(value):
    value = str(value)
    # Official historical wrappers interpolate arguments into shell commands.
    # Private staged links hide arbitrary source filenames safely.
    if not SAFE.fullmatch(value):
        raise ValueError('Official Suite shell arguments require paths without spaces/metacharacters: ' + value)
    return value


def runtime_env(runtime, base=None):
    env = dict(os.environ if base is None else base)
    env.update(runtime.get('env', {}))
    extra = runtime.get('extra_path', [])
    if not isinstance(extra, list):
        raise ValueError('extra_path must be a list')
    bins = [str(Path(runtime.get(k, runtime['aa_python'])).absolute().parent)
            for k in ('suite_python', 'aa_python', 'ac_python')]
    env['PATH'] = os.pathsep.join(bins + extra + [env.get('PATH', '')])
    env['AS_NO_REPO_CHECK'] = '1'
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    return env


def supporting_tools(env):
    """Lightweight provenance for tools outside pip; no patient data read."""
    records={}
    for name,flag in (('samtools','--version'),('Rscript','--version'),('bash','--version')):
        executable=shutil.which(name,path=env.get('PATH'))
        if not executable:
            records[name]={'available':False}
            continue
        record={'available':True,'path':str(Path(executable).absolute()),
                'resolved_path':str(Path(executable).resolve()),'sha256':digest(executable)}
        record['version_probe']=_run([executable,flag],env)
        if name=='bash':
            record['associative_array_probe']=_run([executable,'-c','declare -A probe || exit 1; probe[key]=ok; test "${probe[key]}" = ok'],env)
        records[name]=record
    if records.get('Rscript',{}).get('available'):
        # Match official wrapper invocation; no --vanilla that would change its R paths.
        code='cat("DNAcopy_version=",as.character(packageVersion("DNAcopy")),"\\nDNAcopy_path=",find.package("DNAcopy"),"\\nR_home=",R.home(),"\\n",sep="")'
        records['DNAcopy']={'probe':_run([records['Rscript']['path'],'-e',code],env)}
    return records


def require_ac_shell(ac_source, support_tools):
    script=Path(ac_source)/'make_input.sh'
    if not script.is_file():
        raise ValueError('Official AC make_input.sh is missing')
    if re.search(r'\bdeclare\s+-A\b',script.read_text()):
        bash=support_tools.get('bash',{})
        if not bash.get('available') or bash.get('associative_array_probe',{}).get('returncode')!=0:
            raise ValueError('AC make_input.sh requires Bash associative arrays (Bash 4+); configure a verified modern Bash in runtime extra_path')


class PipelineStageError(RuntimeError):
    def __init__(self,stage,message):
        super().__init__(message)
        self.stage=stage


def verify_full_runtime(runtime, lock_entry, env=None, verified=None):
    """Extend AA/AC verification with exact Suite and actual CNVkit environment."""
    result = dict(verified or {'ready': True, 'diagnostics': []})
    result['diagnostics'] = list(result.get('diagnostics', []))
    env = runtime_env(runtime, env)
    result['support_tools'] = supporting_tools(env)
    try:
        require_ac_shell(Path(runtime['ac_script']).parent,result['support_tools'])
        component = lock_entry.get('suite')
        if not component:
            raise ValueError('No official Suite bundle exists for this AA/AC-only combination')
        script = Path(runtime['suite_script']).resolve()
        source = script.parent
        info = _check_source(source, component)
        if script.name != component['entrypoint']:
            raise ValueError('Suite script does not match locked entrypoint')
        py = str(Path(runtime['suite_python']).absolute())
        if py != str(Path(runtime['ac_python']).absolute()):
            raise ValueError('Suite python must equal AC python: wrapper launches AC and CNVkit with --python3_path')
        for value in (py, script, runtime['aa_python'], runtime['aa_script'], runtime['ac_script'], runtime['cnvkit_script']):
            _safe(value)
        for flag, key in (('--version', 'version_probe'), ('--help', 'help_probe')):
            info[key] = _run([py, str(script), flag], env)
            if info[key]['returncode']:
                raise ValueError(f'Suite {flag} failed: {info[key]["stderr"][-1200:]}')
        text = info['version_probe']['stdout'] + info['version_probe']['stderr']
        if not re.search(r'(?<![\w.])' + re.escape(component['version']) + r'(?![\w.])', text):
            raise ValueError('Actual Suite version differs from lock: ' + text)
        info['version'] = component['version']
        info['help'] = info['help_probe']['stdout'] + info['help_probe']['stderr']
        info['offline_policy'] = offline_policy(component['commit'], info['help'])
        result['suite'] = info
        ck = Path(runtime['cnvkit_script']).resolve()
        if ck.name != 'cnvkit.py' or not ck.is_file():
            raise ValueError('cnvkit_script must point to the actual installed cnvkit.py')
        probe = _run([py, str(ck), 'version'], env)
        actual = probe['stdout'].strip()
        expected = runtime.get('cnvkit_version')
        if not expected or actual != expected or probe['returncode']:
            raise ValueError(f'CNVkit expected explicit {expected!r}, actual {actual!r}; {probe["stderr"][-500:]}')
        code = "import cnvlib,importlib.metadata as m,json,pathlib,hashlib; root=pathlib.Path(cnvlib.__file__).parent; print(json.dumps({'origin':cnvlib.__file__,'version':m.version('CNVkit'),'source_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob('*.py')}}))"
        origin = _run([py, '-c', code], env)
        imported = json.loads(origin['stdout'])
        if imported['version'] != expected:
            raise ValueError('CNVkit imported package version differs from command')
        result['cnvkit'] = {'version': actual, 'script': str(ck), 'sha256': digest(ck), 'package': imported,
                            'version_probe': probe, 'interpretation': 'Explicit configured version; not inferred from Suite tag'}
    except (OSError, ValueError, KeyError, TypeError) as exc:
        result['ready'] = False
        result['diagnostics'].append('full_suite: ' + str(exc))
    return result


def private_reference(source, target):
    """Private directory tree and private indexes/control files; immutable large data links.

    Official Suite writes reduced FASTA/CNV references under its output directory.
    Existing local FASTA indexes prevent faidx writes through a reference link.
    """
    source, target = Path(source), Path(target)
    target.mkdir(parents=True, exist_ok=True)
    audit = []
    def walk(directory, ancestors=()):
        real=directory.resolve()
        if real in ancestors:
            raise ValueError('Cyclic reference directory links are unsupported')
        for child in directory.iterdir():
            yield child
            if child.is_dir():
                yield from walk(child, ancestors+(real,))
    for item in walk(source):
        relative = item.relative_to(source)
        if item.name == 'coverage.stats':
            continue
        dest = target/relative
        if item.is_dir():
            dest.mkdir(parents=True, exist_ok=True)
        elif item.is_file():
            dest.parent.mkdir(parents=True, exist_ok=True)
            if item.name == 'file_list.txt' or item.suffix in ('.fai','.dict','.json','.txt'):
                shutil.copy2(item,dest)
            else:
                dest.symlink_to(item.resolve())
            st=item.stat()
            audit.append({'path':str(item.resolve()),'size':st.st_size,'mtime_ns':st.st_mtime_ns,
                          'sha256':digest(item) if st.st_size<1024*1024 else None})
    atomic_json(target.parent/'reference_audit_before.json',audit)
    return target


def audit_reference(attempt):
    path=Path(attempt)/'reference_audit_before.json'
    if not path.exists():
        return
    changed=[]
    for old in json.loads(path.read_text()):
        item=Path(old['path']);st=item.stat()
        if st.st_size!=old['size'] or st.st_mtime_ns!=old['mtime_ns'] or (old['sha256'] and digest(item)!=old['sha256']):
            changed.append(str(item))
    atomic_json(Path(attempt)/'reference_audit_after.json',{'unchanged':not changed,'changed':changed,
               'method':'size+mtime for all files; full SHA256 for files smaller than1MiB'})
    if changed:
        raise RuntimeError('Original reference resource changed during execution; see audit')


def stage_alignment(source, index, target, fasta):
    import pysam
    source, index, target = Path(source), Path(index), Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.suffix.lower() == '.cram':
        pysam.view('-b', '-T', fasta, '-o', str(target), str(source), catch_stdout=False)
        pysam.index(str(target))
    else:
        target.symlink_to(source.resolve())
        suffix = '.csi' if index.suffix == '.csi' else '.bai'
        Path(str(target) + suffix).symlink_to(index.resolve())
    return str(target)


def build_command(config, runtime, verified, attempt, tumor, normal):
    attempt = Path(attempt)
    params = config.get('full_suite', {})
    command = [runtime['suite_python'], runtime['suite_script'], '-s', 'sample', '-o', str(attempt/'pipeline'),
               '-t', str(config.get('threads', 2)), '--bam', tumor, '--normal_bam', normal,
               '--ref', config['reference']['build'], '--run_AA', '--run_AC',
               '--AA_src', str(Path(runtime['aa_script']).parent),
               '--aa_python_interpreter', runtime['aa_python'], '--python3_path', runtime['ac_python'],
               '--cnvkit_dir', runtime['cnvkit_script'], '--AA_runmode', 'FULL', '--AA_extendmode', 'EXPLORE',
               '--downsample', str(config.get('downsample', -1)),
               '--cngain', str(params.get('cngain', 4.5)), '--cnsize_min', str(params.get('cnsize_min', 50000)),
               '--cnvkit_segmentation', params.get('segmentation', 'cbs')]
    helptext = verified.get('suite', {}).get('help', '')
    if '--no_repo_check' in helptext:
        command += ['--no_repo_check']
    elif verified.get('suite',{}).get('offline_policy') != 'audited_intrinsically_offline_BAM_flow':
        raise ValueError('No verified offline policy for Suite invocation')
    if '--no_cstats' in helptext:
        command += ['--no_cstats']
    solver = runtime.get('solver', config.get('solver', 'clarabel'))
    if '--AA_solver' in helptext:
        command += ['--AA_solver', solver]
    elif solver != 'mosek':
        raise ValueError('Historical Suite cannot pass requested AA solver; no substitution allowed')
    # Deliberately no cnv_bed/no_filter/completed_AA_runs/force/upload options.
    for value in command:
        _safe(value)
    if config['reference']['build'] == 'GRCh38_viral':
        raise ValueError('Official wrapper ignores matched normal for GRCh38_viral; paired-normal full-pipeline scope unsupported')
    return command


def read_seed_bed(path):
    """A present empty bed is meaningful only after verified upstream success."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError('Filtered seed BED is missing')
    seeds = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
            continue
        fields = line.split()
        if len(fields) < 3:
            raise ValueError('Invalid filtered seed BED row')
        start, end = int(fields[1]), int(fields[2])
        if not 0 <= start < end:
            raise ValueError('Invalid filtered seed BED coordinates')
        seeds.append([fields[0], start, end])
    return seeds


def read_cns_segment_count(path):
    """A nonempty but malformed CNV table cannot establish completed CNV calling."""
    with Path(path).open() as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        if not {'chromosome','start','end','log2'}.issubset(reader.fieldnames or []):
            raise ValueError('CNVkit CNS schema invalid')
        count=0
        for row in reader:
            try:
                start,end=int(row['start']),int(row['end'])
                valid=bool(row['chromosome']) and 0<=start<end and math.isfinite(float(row['log2']))
            except (ValueError,TypeError,KeyError):
                valid=False
            if not valid:
                raise ValueError(f'CNVkit CNS has malformed segment at data row {count+1}')
            count+=1
        if not count:
            raise ValueError('CNVkit CNS has no segments; cannot classify as pipeline negative')
        return count


def infer_stage(log):
    """Locate the latest observed stage; do not infer an unsupported root cause."""
    markers = [('cnvkit', r'Running CNVKit|cnvkit\.py batch|cnvkit\.py segment'),
               ('seed_processing', r'Running amplified_intervals|amplified_intervals\.py --ref'),
               ('aa', r'Running AA\b|AmpliconArchitect\.py --ref'),
               ('ac', r'Running AC\b|amplicon_classifier\.py -i'),
               ('result_reading', r'Running AC make_results_table')]
    observations = [(m.start(), stage) for stage, pattern in markers for m in re.finditer(pattern, log)]
    return max(observations)[1] if observations else 'input_check'


def _all_logs(attempt):
    paths = [Path(attempt)/'suite.log']
    paths += list((Path(attempt)/'pipeline').glob('*log*'))
    paths += list((Path(attempt)/'pipeline'/'sample_classification').glob('*.log'))
    return '\n'.join(p.read_text(errors='replace') for p in paths if p.is_file())


def parse_full_result(attempt, sample, runtime, verified):
    from .parsing import parse_run, aa_amplicon_count, LABELS
    attempt = Path(attempt); pipeline = attempt/'pipeline'
    finish = pipeline/'sample_finish_flag.txt'
    metadata_path = pipeline/'sample_run_metadata.json'
    cns = pipeline/'sample_cnvkit_output'/'tumor.cns'
    seedfile = pipeline/'sample_AA_CNV_SEEDS.bed'
    log = _all_logs(attempt)
    if re.search(r'declare: -A: invalid option|make_input\.sh:.*syntax error: operand expected',log):
        raise PipelineStageError('ac','AC input generation failed: make_input.sh ran with a shell lacking Bash associative arrays; empty AC output is not a negative call')
    if 'Traceback (most recent call last)' in log:
        raise RuntimeError('A Python traceback appears in official pipeline logs; wrapper exit zero alone cannot establish successful subprocess completion')
    if not finish.is_file() or finish.read_text().strip() != 'All stages completed':
        raise RuntimeError('Official Suite successful finish flag absent; latest observed stage: ' + infer_stage(log))
    if not metadata_path.is_file():
        raise FileNotFoundError('Official run metadata missing')
    metadata = json.loads(metadata_path.read_text())
    if verified['suite']['version']=='1.2.2' and 'AmpliconSuite-pipeline_version' not in metadata:
        metadata['AmpliconSuite-pipeline_version']=metadata.get('PAA_version')
    checks = {'AmpliconSuite-pipeline_version': verified['suite']['version'],
              'AA_version': verified['aa']['version'], 'AC_version': verified['ac']['version'],
              'cnvkit_version': verified['cnvkit']['version']}
    for key, expected in checks.items():
        if not re.search(r'(?<![\w.])' + re.escape(expected) + r'(?![\w.])', str(metadata.get(key, ''))):
            raise ValueError(f'Official run metadata {key} does not match the verified runtime')
    for key in ('cnvkit_cmd', 'amplified_intervals_cmd'):
        if metadata.get(key) in (None, '', 'NA'):
            raise ValueError('Required full-pipeline step was not executed: ' + key)
    if '--normal ' not in metadata['cnvkit_cmd']:
        raise ValueError('CNVkit did not record matched-normal analysis')
    if not cns.is_file():
        raise FileNotFoundError('CNVkit segmentation output is missing')
    read_cns_segment_count(cns)
    seeds = read_seed_bed(seedfile)
    for name, target in [('aa', pipeline/'sample_AA_results'), ('ac', pipeline/'sample_classification')]:
        link = attempt/name
        if not link.exists() and target.exists():
            link.symlink_to(target.resolve(), target_is_directory=True)
    aa_called = metadata.get('AA_cmd') not in (None, '', 'NA')
    ac_called = metadata.get('AC_cmd') not in (None, '', 'NA')
    count = aa_amplicon_count(attempt) if aa_called else None
    stages = {'input_check': {'status':'complete'}, 'cnvkit': {'status':'complete'},
              'seed_processing': {'status':'complete','seed_count':len(seeds)},
              'aa': {'status':'complete' if aa_called else 'not_run_no_seeds'},
              'ac': {'status':'complete' if ac_called else 'not_run_no_amplicons'}}
    atomic_json(attempt/'stages.json', stages)
    supported = verified.get('supported_labels', [])
    if not seeds:
        if count not in (None, 0):
            raise ValueError('Empty filtered seeds conflict with nonzero AA reconstruction')
        if aa_called and count is None:
            raise FileNotFoundError('AA was invoked with empty seeds but its completion summary is missing')
        result = empty_record(sample['sample_id'], WHOLE_SAMPLE, runtime['id'], 'not_detected', 'seed_processing',
                              'Completed full pipeline produced no qualifying CNVkit-derived AA seeds; not a fixed-seed AA sensitivity result',
                              detected=False, supported_labels=supported,
                              category_status={x:'not_detected' if x in supported else 'unsupported' for x in LABELS})
    else:
        if count:
            check_ac_runtime(log, verified['ac']['version'])
            if verified['ac']['version'].startswith('2.'):
                ac_log=pipeline/'sample_classification'/'sample.log'
                if not ac_log.is_file() or not re.search(r' - INFO - Complete\s*$',ac_log.read_text(errors='replace')):
                    raise RuntimeError('AC2 final completion log marker is absent; profile presence alone is insufficient')
        result = parse_run(attempt, sample['sample_id'], WHOLE_SAMPLE, runtime['id'], supported,
                           reference_lengths=sample.get('reference_lengths'),
                           allowed_missing_contigs=sample.get('missing_reference_contigs', []))
    execution=json.loads((attempt/'execution.json').read_text()) if (attempt/'execution.json').exists() else {}
    result.update(endpoint='full_pipeline', analysis_scope='full_suite', seed_count=len(seeds), seed_intervals=seeds, seed_stage_evaluable=True,
                  aa_executed=aa_called, ac_executed=ac_called, aa_reconstruction_performed=bool(seeds and count),
                  detection_stage='seed_processing' if not seeds else ('aa' if count==0 else 'ac'),
                  reference_scope='human_only_incomplete_reference' if sample.get('missing_reference_contigs') else 'complete_reference',
                  cnvkit_version=verified['cnvkit']['version'], seed_thresholds=execution.get('seed_thresholds',{'status':'not_recorded_in_fixture'}), suite_version=verified['suite']['version'],
                  runtime_variants=verified.get('runtime_variants', []))
    result.setdefault('evidence', []).extend(map(str, [finish, metadata_path, cns, seedfile, attempt/'suite.log']))
    return result


def partial_seed_evidence(attempt):
    """Evaluate seed output only when official control flow reached AA afterward."""
    pipeline=Path(attempt)/'pipeline'
    log=_all_logs(attempt)
    if not re.search(r'Running AA\b',log):
        return {'seed_stage_evaluable':False}
    try:
        cns=pipeline/'sample_cnvkit_output'/'tumor.cns'
        read_cns_segment_count(cns)
        seedfile=pipeline/'sample_AA_CNV_SEEDS.bed'
        seeds=read_seed_bed(seedfile)
        return {'seed_stage_evaluable':True,'seed_count':len(seeds),'seed_intervals':seeds,
                'seed_stage_evidence':[str(cns),str(seedfile),str(Path(attempt)/'suite.log')]}
    except (OSError,ValueError):
        return {'seed_stage_evaluable':False}


def run_full_suite(config, samples, runtimes, out):
    if any(s.get('status') == 'ok' for s in samples):
        raise ValueError('Whole-genome Suite execution is disabled in AmpliconSuite Cross-Version Concordance; use local candidate validation')
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=True)
    records = []; total = len(samples)*len(config['versions']); done = 0
    for sample in samples:
        for runtime in config['versions']:
            done += 1; version = runtime['id']; verified = runtimes[version]
            print(f'[{done}/{total}] full Suite / {sample["sample_id"]} / {version}', flush=True)
            job = out/'jobs'/sample['sample_id']/WHOLE_SAMPLE/version;job.mkdir(parents=True, exist_ok=True)
            saved = job/'result.json'
            if saved.is_file():
                try:
                    old = json.loads(saved.read_text())
                except (OSError, ValueError):
                    old = {}
                hashes = old.get('output_sha256', {})
                if old.get('status') in ('ok','not_detected','invalid','unsupported') and hashes and all(Path(p).is_file() and digest(p)==h for p,h in hashes.items()):
                    fresh = parse_full_result(old['run_dir'], sample, runtime, verified)
                    fresh.update(run_dir=old['run_dir'],output_sha256=hashes);atomic_json(saved,fresh)
                    records.append(fresh);print('  verified complete full-pipeline job: reuse',flush=True);continue
            if sample.get('status') != 'ok' or not verified.get('ready'):
                record = empty_record(sample['sample_id'],WHOLE_SAMPLE,version,'failed','input_check',
                                      sample.get('reason') or 'Full runtime verification failed',
                                      diagnostics=verified.get('diagnostics',[]),endpoint='full_pipeline',
                                      supported_labels=verified.get('supported_labels',[]))
                record['reference_scope']='human_only_incomplete_reference' if sample.get('missing_reference_contigs') else 'complete_reference'
                record['seed_stage_evaluable']=False
                atomic_json(saved,record);records.append(record);atomic_json(out/'records.json',records);continue
            attempts = [int(p.name.split('_')[-1]) for p in job.glob('attempt_*') if p.is_dir()]
            attempt = job/f'attempt_{max(attempts,default=0)+1:04d}';attempt.mkdir()
            stage='input_check'
            try:
                _safe(attempt)
                tumor=stage_alignment(sample['alignment'],sample['index'],attempt/'alignment'/'tumor.bam',config['reference']['fasta'])
                normal=stage_alignment(sample['normal_alignment'],sample['normal_index'],attempt/'alignment'/'normal.bam',config['reference']['fasta'])
                env=runtime_env(runtime)
                env.update({'AA_DATA_REPO':str(private_reference(config['reference']['aa_data_repo'],attempt/'reference')),
                            'AA_SRC':str(Path(runtime['aa_script']).parent),'AC_SRC':str(Path(runtime['ac_script']).parent),
                            'MPLCONFIGDIR':str(attempt/'mplcache'),'OMP_NUM_THREADS':str(config.get('threads',2)),
                            'OPENBLAS_NUM_THREADS':str(config.get('threads',2)), 'NUMEXPR_MAX_THREADS':str(config.get('threads',2))})
                command=build_command(config,runtime,verified,attempt,tumor,normal)
                execution={'argv':command,'started':now(),'endpoint':'full_pipeline',
                            'scope':'Fresh paired-normal CNVkit, native seed processing, AA FULL/EXPLORE, AC',
                            'sources':{k:sample[k] for k in ('alignment','index','normal_alignment','normal_index')},
                            'runtime_variants':verified.get('runtime_variants',[]),
                            'seed_thresholds':{'cngain':config.get('full_suite',{}).get('cngain',4.5),'cnsize_min':config.get('full_suite',{}).get('cnsize_min',50000)},
                            'CNVkit':verified.get('cnvkit'), 'controlled_CNVkit':True,
                            'support_tools_at_start':supporting_tools(env)}
                require_ac_shell(Path(runtime['ac_script']).parent,execution['support_tools_at_start'])
                atomic_json(attempt/'execution.json',execution)
                stage='cnvkit'
                rc=execute(command,attempt,env,attempt/'suite.log',timeout=config.get('job_timeout_seconds'))
                execution.update(ended=now(),returncode=rc,support_tools_at_end=supporting_tools(env))
                atomic_json(attempt/'execution.json',execution)
                stage=infer_stage(_all_logs(attempt))
                audit_reference(attempt)
                if rc:
                    raise RuntimeError(f'Official Suite exited {rc}; see suite.log; root cause pending unless logged explicitly')
                record=parse_full_result(attempt,sample,runtime,verified)
            except KeyboardInterrupt:
                stage=infer_stage(_all_logs(attempt));atomic_json(saved,empty_record(sample['sample_id'],WHOLE_SAMPLE,version,'failed',stage,'Interrupted; completed sample/version jobs preserved',endpoint='full_pipeline'))
                atomic_json(attempt/'stages.json',{stage:{'status':'interrupted'}});raise
            except Exception as exc:
                stage=getattr(exc,'stage',stage)
                record=empty_record(sample['sample_id'],WHOLE_SAMPLE,version,'missing' if isinstance(exc,FileNotFoundError) else 'failed',stage,str(exc),endpoint='full_pipeline',supported_labels=verified.get('supported_labels',[]),evidence_paths=[str(attempt/'suite.log')])
                atomic_json(attempt/'stages.json',{stage:{'status':'failed','reason':str(exc)}})
            record['run_dir']=str(attempt)
            record['reference_scope']='human_only_incomplete_reference' if sample.get('missing_reference_contigs') else 'complete_reference'
            if record['status'] in ('failed','missing'):
                record.update(partial_seed_evidence(attempt))
            # Preserve scientific outputs; hash essential text evidence for safe resume.
            files=[p for p in attempt.rglob('*') if p.is_file() and not p.is_symlink() and p.suffix in ('.json','.txt','.tsv','.bed','.cns','.log') and not any(x in p.parts for x in ('reference','alignment','mplcache'))]
            record['output_sha256']={str(p):digest(p) for p in files}
            atomic_json(saved,record);records.append(record);atomic_json(out/'records.json',records)
    table(out/'anomalies.tsv',[{k:r.get(k,'') for k in ('sample_id','seed_id','version','status','stage','reason','diagnostics')} for r in records if r['status'] not in ('ok','not_detected')])
    return records
