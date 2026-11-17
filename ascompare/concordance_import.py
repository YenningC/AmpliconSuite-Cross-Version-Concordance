"""Local small-artifact import and AC-only reclassification for cross-version concordance.

No BAM/CRAM is opened. Imported component identities are declarations, never
runtime verification. AC-only results inherit AA graphs and cannot measure AA.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil

from .common import atomic_json, digest, now
from .full_suite import WHOLE_SAMPLE, private_reference, audit_reference, read_seed_bed
from .parsing import (LABELS, OutputError, _amplicon_files, _associate_beds, _profile_labels,
                      _profile_rows, _summary_count, aa_amplicon_count, parse_run)
from .runner import check_ac_runtime, empty_record, execute, prepare_ac_input
from .versions import _check_source, _run

MAX_FILE_BYTES=64*1024**2
MAX_TOTAL_BYTES=256*1024**2


def _new_output(destination, sources):
    out=Path(destination).absolute()
    for source in sources:
        p=Path(source).resolve()
        if p==out.resolve() or p in out.resolve().parents or out.resolve() in p.parents:
            raise ValueError('Concordance output must be separate from original inputs/resources')
    if out.exists():raise ValueError('Concordance import/reclassification requires a new output directory')
    out.mkdir(parents=True)
    return out


def _artifacts(aa, ac=None):
    """Bounded, explicit globs only; never recursively walk source data trees."""
    paths=[]
    for pattern in ('*_summary.txt','*_amplicon*_graph.txt','*_amplicon*_cycles.txt'):
        paths.extend((p,Path('aa')/p.name) for p in sorted(Path(aa).glob(pattern)))
    if ac is not None:
        ac=Path(ac)
        paths.extend((p,Path('ac')/p.name) for p in sorted(ac.glob('*_amplicon_classification_profiles.tsv')))
        paths.extend((p,Path('ac')/p.name) for p in sorted(ac.glob('*_features_to_graph.txt')))
        for directory in sorted(ac.glob('*_classification_bed_files')):
            if directory.is_symlink():raise ValueError('Linked classification directory is not an importable small artifact')
            if directory.is_dir():
                paths.extend((p,Path('ac')/directory.name/p.name) for p in sorted(directory.glob('*_intervals.bed')))
    return paths


def _validate_artifacts(paths):
    total=0
    for path,_ in paths:
        if path.is_symlink() or not path.is_file():raise ValueError('Linked or non-regular result artifact refused')
        size=path.stat().st_size;total+=size
        if size>MAX_FILE_BYTES or total>MAX_TOTAL_BYTES:
            raise ValueError('Selected result text exceeds the small-artifact copy limit')


def _copy_artifacts(paths, out):
    _validate_artifacts(paths)
    provenance=[]
    for path,relative in paths:
        before=digest(path);target=out/relative;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,target)
        if digest(target)!=before or digest(path)!=before:raise ValueError('Input artifact changed during import')
        provenance.append(dict(source=str(path.resolve()),relative_path=str(relative),sha256=before,bytes=target.stat().st_size))
    return provenance


def _source_paths(spec):
    base = Path(spec['run_dir']) if spec.get('run_dir') else None
    aa = Path(spec['aa_dir']) if spec.get('aa_dir') else base / 'aa' if base else None
    ac = Path(spec['ac_dir']) if spec.get('ac_dir') else base / 'ac' if base else None
    stage = Path(spec['stages_json']) if spec.get('stages_json') else base / 'stages.json' if base else None
    return base, aa, ac, stage


def inspect_import_format(spec, supported_labels=None):
    """Inspect small native output headers and associations without writing files.

    Counts describe file structure. They do not validate interval contents,
    runtime success, BAMs, or sample comparability.
    """
    _, aa, ac, stage_path = _source_paths(spec)
    result = dict(status='missing', format_only=True, summary_amplicons=None, graph_files=0,
                  cycles_files=0, profile_columns=[], profile_rows=0, matched_amplicons=0,
                  feature_bed_files=0, issues=[])

    def issue(status, stage, reason, amp=None):
        item = dict(status=status, stage=stage, reason=reason)
        if amp is not None:
            item['amplicon_id'] = amp
        result['issues'].append(item)

    try:
        if aa is None or not aa.is_dir():
            raise OutputError('AA result directory missing', 'missing', 'aa')
        if ac is None:
            raise OutputError('AC result directory not specified', 'missing', 'ac')
        artifact_paths = [path for path, _ in _artifacts(aa, ac)]
        if stage_path and stage_path.is_file():
            artifact_paths.append(stage_path)
        elif spec.get('stages_json'):
            issue('missing', 'result_reading', 'Explicit imported stage evidence is missing')
        _validate_artifacts([(p, None) for p in artifact_paths])
        stages = json.loads(stage_path.read_text()) if stage_path in artifact_paths else spec.get('declared_completion')
        if stages is not None:
            if not isinstance(stages, dict) or any(not isinstance(stages.get(name, {}), dict)
                                                    for name in ('input_check', 'seed_processing', 'aa', 'ac')):
                raise OutputError('Imported stage evidence is not a stage-record mapping')
            for name in ('input_check', 'seed_processing', 'aa', 'ac'):
                state = stages.get(name, {})
                if state.get('status') in ('failed', 'interrupted', 'unsupported'):
                    issue(state['status'], name, state.get('reason') or state.get('error') or 'Execution did not complete')
        summaries = [p for p in artifact_paths if p.parent == aa and p.name.endswith('_summary.txt')]
        graphs = [p for p in artifact_paths if p.parent == aa and p.name.endswith('_graph.txt')]
        cycles = [p for p in artifact_paths if p.parent == aa and p.name.endswith('_cycles.txt')]
        profiles = [p for p in artifact_paths if p.parent == ac and p.name.endswith('_amplicon_classification_profiles.tsv')]
        result.update(graph_files=len(graphs), cycles_files=len(cycles),
                      feature_bed_files=sum(p.name.endswith('_intervals.bed') for p in artifact_paths))
        expected = _summary_count(summaries[0]) if len(summaries) == 1 else None
        result['summary_amplicons'] = expected
        if expected is None:
            issue('missing', 'aa', 'AA summary missing, duplicated, or count unrecognized')
        elif expected != len(graphs) or expected != len(cycles):
            issue('missing', 'aa', 'AA summary count does not match graph/cycles file counts')
        if expected == 0:
            aa_stage, ac_stage = (stages or {}).get('aa', {}), (stages or {}).get('ac', {})
            if aa_stage.get('status') != 'complete' or ac_stage.get('status') != 'skipped_zero_amplicons' \
                    or aa_stage.get('returncode', 0) != 0:
                issue('missing', 'aa', 'Zero amplicons require explicit AA completion and AC skipped-zero evidence')
        else:
            _inspect_amplicon_names(ac, supported_labels, result, issue, graphs, cycles, profiles, expected)
    except OutputError as exc:
        issue(exc.status, exc.stage, str(exc))
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as exc:
        issue('invalid', 'result_reading', f'{type(exc).__name__}: {exc}')
    issues = result['issues']
    result['status'] = ('partial' if result['matched_amplicons'] else issues[0]['status']) if issues else 'ok'
    return result


def _inspect_amplicon_names(ac, supported_labels, result, issue, graphs, cycles, profiles, expected):
    graphs = _amplicon_files(graphs, 'graph')
    cycles = _amplicon_files(cycles, 'cycles')
    if len(profiles) != 1:
        raise OutputError('Expected one AC classification profile', 'missing', 'ac')
    profile = profiles[0]
    with profile.open(newline='') as handle:
        result['profile_columns'] = next(csv.reader(handle, delimiter='\t'), [])
    supported = set(supported_labels or LABELS)
    if not supported_labels and 'FAN+' not in result['profile_columns']:
        supported.remove('FAN')
    fields, rows = _profile_rows(profile, supported)
    result['profile_rows'] = len(rows)
    if len(rows) != expected or set(rows) != set(graphs):
        issue('missing', 'ac', 'AC rows do not match reconstructed amplicons exactly')
    beds, aliases, map_issues, _ = _associate_beds(ac, profile, graphs)
    for reason in map_issues:
        issue('invalid', 'ac', reason)
    for amp in sorted(set(graphs) | set(cycles) | set(rows)):
        graph, cycle, row = graphs.get(amp), cycles.get(amp), rows.get(amp)
        if not graph or not cycle or graph[1] != cycle[1] or not all(f[0].stat().st_size for f in (graph, cycle)):
            issue('missing', 'aa', 'Missing or empty paired graph/cycles artifact', amp)
        elif row is None:
            issue('missing', 'ac', 'Missing AC profile row', amp)
        elif None in row or any(row.get(field) is None for field in fields):
            issue('invalid', 'ac', 'Malformed AC profile row', amp)
        else:
            prefix = Path(row['sample_name']).name
            if prefix != graph[1] and (amp, prefix) not in aliases:
                issue('invalid', 'ac', 'AC sample prefix does not match AA graph or native mapping', amp)
                continue
            try:
                labels = _profile_labels(row, fields, supported)
                supported_beds = [(feature.rsplit('_', 1)[0], bed) for feature, bed in beds.get(amp, [])
                                  if not feature.lower().startswith('unknown')]
                if supported_beds and row['amplicon_decomposition_class'] == 'No-FSCNA':
                    raise OutputError('No-FSCNA profile contradicts classified feature BED')
                if {label for label, _ in supported_beds} - supported:
                    raise OutputError('Unrecognized class on a feature BED', 'unsupported')
                available = {label for label, bed in supported_beds if bed.stat().st_size}
                if labels - available:
                    raise OutputError('Positive class lacks a nonempty named feature BED', 'missing')
                result['matched_amplicons'] += 1
            except OutputError as exc:
                issue(exc.status, 'ac', str(exc), amp)


def _aa_identity(out):
    hashes={p.name:digest(p) for p,_ in _artifacts(Path(out)/'aa')}
    identity=hashlib.sha256(json.dumps(hashes,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return hashes,identity


def _component_versions(spec):
    supplied=spec.get('component_versions',{})
    result={}
    for key in ('aa','ac'):
        value=spec.get(key+'_version',supplied.get(key,supplied.get(key.upper())))
        if isinstance(value,dict):value=value.get('version')
        result[key]=str(value) if value is not None else None
    return result


def _save(out, record, provenance):
    record['run_dir']=str(out);record['input_artifacts']=provenance
    hashes,identity=_aa_identity(out)
    record['aa_artifact_sha256']=hashes;record['source_graph_identity']=identity
    atomic_json(out/'provenance.json',{'created_utc':now(),'analysis_kind':record.get('analysis_kind'),
        'provenance_kind':record.get('provenance_kind'),'completion_evidence':record.get('completion_evidence'),
        'input_artifacts':provenance,'source_graph_identity':identity,'recorded_components':record.get('component_versions'),
        'source_aa_version':record.get('source_aa_version'),'original_execution_verified':False,
        'ac_runtime_verified':record.get('ac_runtime_verified',False),'ac_executed':record.get('ac_executed',False),
        'source_metadata':record.get('source_metadata')})
    record['output_sha256']={str(p):digest(p) for p in out.rglob('*') if p.is_file() and not p.is_symlink()
                             and not any(x in p.relative_to(out).parts for x in ('reference','mplcache'))}
    atomic_json(out/'result.json',record)
    return record


def import_result(spec, sample_id, version, destination, supported_labels, reference_lengths=None):
    """Import selected AA/AC artifacts, preserving execution uncertainty.

    spec: aa_dir/ac_dir OR run_dir; optional seeds_bed/reference_build,
    aa_version/ac_version or component_versions, stages_json, and
    declared_completion (a stages-shaped mapping). Validated declared component
    matching is supplied by the caller as validated_declared_components.
    """
    base, aa, ac, stage_path = _source_paths(spec)
    sources=[p for p in (base,aa,ac,spec.get('seeds_bed'),spec.get('stages_json')) if p]
    out=_new_output(destination,sources)
    support_known=bool(supported_labels)
    support=list(supported_labels or [])
    components=_component_versions(spec);provenance=[]
    declaration_check=spec.get('validated_declared_components')
    declarations_valid=(declaration_check is True or (isinstance(declaration_check,dict)
        and all(components.values()) and _component_versions({'component_versions':declaration_check})==components))
    common=dict(analysis_kind='aa_ac',provenance_kind='imported_declared',component_versions=components,
        aa_version=components['aa'],ac_version=components['ac'],source_aa_version=components['aa'],
        reference_build=spec.get('reference_build'),reference_lengths=reference_lengths,original_execution_verified=False,
        ac_runtime_verified=False,
        declared_components_validated=declarations_valid,
        category_support_known=support_known,comparison_eligible=False)
    try:
        if aa is None or not aa.is_dir():raise FileNotFoundError('AA result directory missing')
        if ac is None:raise FileNotFoundError('AC result directory not specified')
        provenance=_copy_artifacts(_artifacts(aa,ac),out)
        (out/'aa').mkdir(exist_ok=True);(out/'ac').mkdir(exist_ok=True)
        if spec.get('stages_json') and not stage_path.is_file():
            raise FileNotFoundError('Explicit imported stage evidence is missing')
        if spec.get('declared_completion') is not None and not isinstance(spec['declared_completion'],dict):
            raise ValueError('declared_completion must be an explicit AA/AC stage mapping')
        if stage_path and stage_path.is_file():
            provenance+=_copy_artifacts([(stage_path,Path('source_stages.json'))],out)
            stages=json.loads((out/'source_stages.json').read_text())
            completion='imported_stage_record_unverified'
        elif isinstance(spec.get('declared_completion'),dict):
            stages=spec['declared_completion'];completion='external_completion_declaration'
        else:
            stages={'aa':{'status':'complete','basis':'artifact_consistency_only'},
                    'ac':{'status':'complete','basis':'artifact_consistency_only'}}
            completion='artifact_consistency_only'
        if not isinstance(stages,dict) or any(not isinstance(stages.get(name,{}),dict) for name in ('input_check','seed_processing','aa','ac')):
            raise ValueError('Imported stages must be a mapping of named stage records')
        common['completion_evidence']=completion
        atomic_json(out/'stages.json',stages)
        count=aa_amplicon_count(out)
        if count==0 and (completion=='artifact_consistency_only' or stages.get('aa',{}).get('status')!='complete'
                or stages.get('ac',{}).get('status')!='skipped_zero_amplicons'):
            raise FileNotFoundError('Zero-amplicon import needs explicit AA completion and AC skipped-zero stage evidence')
        # Unknown capability supports descriptive parsing only, not a declaration
        # that the version implements every category present in current schemas.
        parse_support=support if support_known else [x for x in LABELS if x!='FAN']
        if not support_known:
            profiles=list((out/'ac').glob('*_amplicon_classification_profiles.tsv'))
            if any('FAN+' in p.read_text().splitlines()[0].split('\t') for p in profiles if p.stat().st_size):
                parse_support.append('FAN')
        record=parse_run(out,sample_id,spec.get('seed_id',WHOLE_SAMPLE),version,parse_support,reference_lengths)
        if not support_known:
            record.update(supported_labels=[],unsupported_labels=[],category_status={x:'support_unknown' for x in LABELS})
        if spec.get('seeds_bed'):
            seed_path=Path(spec['seeds_bed'])
            provenance+=_copy_artifacts([(seed_path,Path('seeds.bed'))],out)
            record['seed_intervals']=read_seed_bed(out/'seeds.bed')
            if reference_lengths is not None and any(c not in reference_lengths or e>reference_lengths[c] for c,s,e in record['seed_intervals']):
                raise ValueError('Imported seed is outside the declared reference coordinates')
        record.update(common)
        record['comparison_eligible']=bool(record['status'] in ('ok','not_detected') and support_known
            and all(components.values()) and declarations_valid)
        record['comparison_limitations']=['Imported component identities are declarations matched to a registry, not a verified historical execution.',
            'Compatibility applies only to imported output metrics; parameters/resources may be incompletely documented.']
        if not record['comparison_eligible']:record['comparison_limitations'].append('Component combination, category support or complete output evidence is insufficient for comparison.')
    except (OSError,ValueError,KeyError,TypeError) as exc:
        record=empty_record(sample_id,spec.get('seed_id',WHOLE_SAMPLE),version,'missing' if isinstance(exc,FileNotFoundError) else 'invalid',
            'result_reading',str(exc),supported_labels=support,detected=None,**common)
    return _save(out,record,provenance)


def _ac_environment(runtime):
    env=dict(os.environ);env.update(runtime.get('env',{}))
    extra=runtime.get('extra_path',[])
    if not isinstance(extra,list):raise ValueError('AC extra_path must be a list')
    env['PATH']=os.pathsep.join([str(Path(runtime['ac_python']).absolute().parent),*map(str,extra),env.get('PATH','')])
    env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1')
    return env


def _source_metadata(source_record, graph_identity=None):
    aa=source_record.get('aa_version') or source_record.get('source_aa_version') or source_record.get('component_versions',{}).get('aa')
    return {**{k:source_record.get(k) for k in ('sample_id','source_version','reference_build','reference_scope',
            'component_identity_level','provenance_kind','declared_components_validated','category_support_known','comparison_eligible')},
            'source_aa_version':aa,'source_graph_identity':graph_identity}


def _ac_parameters(config):
    requested=config.get('ac_parameters',{})
    if not isinstance(requested,dict) or set(requested)-{'min_size','min_flow'}:
        raise ValueError('Unsupported AC parameter configuration; no parameter is silently ignored')
    parameters={'min_size':requested.get('min_size',config.get('ac_min_size',50000)),
                'min_flow':requested.get('min_flow',config.get('ac_min_flow',1.0))}
    if any(not math.isfinite(float(v)) or float(v)<0 for v in parameters.values()):
        raise ValueError('AC classification parameters must be finite and nonnegative')
    if not float(parameters['min_size']).is_integer():raise ValueError('AC minimum size must be an integer')
    return parameters


def _probe_json(stdout):
    """Solver imports may print license/status lines around the JSON payload."""
    for line in reversed(stdout.splitlines()):
        try:value=json.loads(line)
        except ValueError:continue
        if isinstance(value,dict):return value
    raise ValueError('Dependency probe produced no JSON object')


def _verify_ac_runtime(runtime, entry, env):
    """Verify AC/BFB only; no AA interpreter, BAM or AA solver is required."""
    component=entry['ac'];script=Path(runtime['ac_script']).resolve();root=script
    for _ in Path(component['entrypoint']).parts:root=root.parent
    if script!=root/component['entrypoint']:raise ValueError('AC entrypoint differs from locked source')
    info=_check_source(root,component);python=str(Path(runtime['ac_python']).absolute())
    info.update(version=component['version'],python=python,sha256=digest(script))
    for flag,key in (('--version','version_probe'),('--help','help_probe')):
        info[key]=_run([python,str(script),flag],env)
        if info[key]['returncode']:raise ValueError('AC '+flag+' failed: '+info[key]['stderr'][-1200:])
    text=info['version_probe']['stdout']+info['version_probe']['stderr']
    if not re.search(r'(?<![\w.])v?'+re.escape(component['version'])+r'(?![\w.])',text):
        raise ValueError('Actual AC version differs from the lock')
    info['help']=info['help_probe']['stdout']+info['help_probe']['stderr']
    info['packages']=_run([python,'-m','pip','freeze','--all'],env)
    if info['packages']['returncode']:raise ValueError('Cannot record AC dependency versions')
    supported=entry.get('supported_labels',entry.get('supported_classes'))
    if not isinstance(supported,list) or not supported:raise ValueError('Locked AC category support is unknown')
    verified={'ac':info,'supported_labels':list(supported),'validation_level':'AC launch and dependency checks only'}
    if entry.get('bfb'):
        code='''import importlib.metadata as m,json,pathlib,hashlib,base64
d=m.distribution('bfbarchitect');import bfbarchitect
from bfbarchitect._version import __version__
from bfbarchitect import detect_solver
r={'version':d.version,'module_version':__version__,'changed':[],'origin_matches':pathlib.Path(bfbarchitect.__file__).resolve()==pathlib.Path(d.locate_file('bfbarchitect/__init__.py')).resolve(),'selected_backend':detect_solver()}
for f in d.files or []:
 if str(f).endswith('.py') and f.hash and f.hash.mode=='sha256':
  p=d.locate_file(f)
  if base64.urlsafe_b64encode(hashlib.sha256(p.read_bytes()).digest()).decode().rstrip('=')!=f.hash.value:r['changed'].append(str(f))
print(json.dumps(r))'''
        probe=_run([python,'-c',code],env);verified['bfb_probe']=probe
        if probe['returncode']:raise ValueError('Required BFBArchitect cannot be imported/verified')
        bfb=_probe_json(probe['stdout'])
        if bfb.get('version')!=entry['bfb']['version'] or bfb.get('module_version')!=entry['bfb']['version'] or bfb.get('changed') or not bfb.get('origin_matches'):
            raise ValueError('BFBArchitect version/source integrity differs from the lock')
        probe=_run([python,'-c','from bfbarchitect.BFBSolver import reconstruct_BFB_cbc; print(reconstruct_BFB_cbc([2],[0],[1],1,max_time=15,max_threads=1))'],env)
        verified['bfb_fallback_probe']=probe
        if probe['returncode']:raise ValueError('Required BFB CBC fallback cannot execute')
    return verified


def reclassify_result(source_record, runtime, entry, config, destination):
    """Run verified AC on immutable inherited AA graphs; AA is never executed."""
    source=Path(source_record['run_dir']);reference=config.get('reference',{})
    out=_new_output(destination,[source,reference['aa_data_repo']] if reference.get('aa_data_repo') else [source])
    sample=source_record['sample_id'];seed=source_record.get('seed_id',WHOLE_SAMPLE);version=runtime.get('id',entry.get('id'))
    source_aa=source_record.get('aa_version') or source_record.get('source_aa_version') or source_record.get('component_versions',{}).get('aa')
    common=dict(analysis_kind='ac_only',provenance_kind='ac_only_attempt',source_aa_version=source_aa,
        aa_version=source_aa,ac_version=entry.get('ac',{}).get('version'),source_record_dir=str(source.resolve()),
        component_versions={'aa':source_aa,'ac':entry.get('ac',{}).get('version')},
        reference_build=reference.get('build'),comparison_eligible=False,
        completion_evidence='inherited_AA_artifacts_and_current_AC_execution',aa_executed=False,ac_executed=False,ac_runtime_verified=False,
        source_metadata=_source_metadata(source_record))
    provenance=[];support=entry.get('supported_labels',entry.get('supported_classes',[])) or [];stage='input_check'
    try:
        if runtime.get('id')!=entry.get('id'):raise ValueError('AC runtime and version lock identifiers differ')
        if not reference.get('aa_data_repo') or not reference.get('build'):raise ValueError('AC needs an explicit local reference resource and build')
        if source_record.get('reference_build') and source_record['reference_build']!=reference['build']:
            raise ValueError('Source AA and target AC reference builds differ')
        paths=_artifacts(source/'aa')
        recorded=source_record.get('aa_artifact_sha256',{})
        raw_hashes={str(Path(p).resolve()):v for p,v in source_record.get('output_sha256',{}).items()}
        if not paths:raise FileNotFoundError('No inherited AA artifacts')
        for path,_ in paths:
            expected=recorded.get(path.name,raw_hashes.get(str(path.resolve())))
            if not expected or digest(path)!=expected:raise ValueError('Inherited AA artifact has missing/mismatched provenance checksum')
        provenance=_copy_artifacts(paths,out);(out/'ac').mkdir()
        common['source_metadata']=_source_metadata(source_record,_aa_identity(out)[1])
        stages={'aa':{'status':'complete','basis':'inherited_verified_artifacts'},'ac':{'status':'pending'}}
        atomic_json(out/'stages.json',stages)
        pre=parse_run(out,sample,seed,version,support,source_record.get('reference_lengths'))
        count=aa_amplicon_count(out)
        if count==0:raise ValueError('AC-only comparison unsupported: source AA has no reconstructed graph to classify')
        if pre.get('status')!='missing' or pre.get('stage')!='ac':
            error=FileNotFoundError if pre.get('status')=='missing' else ValueError
            raise error('Inherited AA reconstruction is incomplete/invalid: '+pre.get('reason',''))
        stage='ac';inputfile=prepare_ac_input(out);env=_ac_environment(runtime)
        env.update(AA_DATA_REPO=str(private_reference(reference['aa_data_repo'],out/'reference')),MPLCONFIGDIR=str(out/'mplcache'))
        verified=_verify_ac_runtime(runtime,entry,env);support=verified['supported_labels']
        common['ac_runtime_verified']=True
        common['provenance_kind']='verified_ac_on_inherited_aa'
        atomic_json(out/'runtime_check.json',verified)
        parameters=_ac_parameters(config)
        if any('--'+k not in verified['ac']['help'] for k in parameters):raise ValueError('Locked AC does not support requested classification parameters')
        command=[str(Path(runtime['ac_python']).absolute()),str(Path(runtime['ac_script']).resolve()),'-i',str(inputfile),
            '--ref',reference['build'],'-o',str(out/'ac/classification')]
        for key,value in parameters.items():command+=['--'+key,str(value)]
        atomic_json(out/'execution.json',{'created_utc':now(),'argv':command,'analysis_kind':'ac_only','parameters':parameters,
            'source_aa_version':source_aa,'source_graph_identity':_aa_identity(out)[1],'AA_was_not_rerun':True})
        common['ac_executed']=True
        rc=execute(command,out,env,out/'ac.log',timeout=config.get('job_timeout_seconds'))
        audit_reference(out)
        if rc:raise RuntimeError(f'AC exited {rc}; inspect ac.log')
        log=(out/'ac.log').read_text(errors='replace')
        native_log=out/'ac/classification.log'
        if native_log.is_file():log+='\n'+native_log.read_text(errors='replace')
        if 'Traceback (most recent call last)' in log:raise RuntimeError('AC traceback despite process exit status')
        check_ac_runtime(log,verified['ac']['version'])
        if verified['ac']['version'].startswith('2.') and (not native_log.is_file() or not re.search(r' - INFO - Complete\s*$',native_log.read_text(errors='replace'))):
            raise RuntimeError('AC2 final completion marker missing')
        for item in provenance:
            if digest(item['source'])!=item['sha256'] or digest(out/item['relative_path'])!=item['sha256']:
                raise ValueError('Source/inherited AA artifact changed during AC execution')
        stages['ac']={'status':'complete','returncode':rc,'basis':'executed'};atomic_json(out/'stages.json',stages)
        record=parse_run(out,sample,seed,version,support,source_record.get('reference_lengths'))
        record.update(common)
        record['comparison_eligible']=record['status'] in ('ok','not_detected') and bool(source_aa) and source_record.get('comparison_eligible') is True
        record['comparison_limitations']=['AC-only classification of identical inherited AA graphs; no AA-version or end-to-end Suite equivalence is tested.']
    except KeyboardInterrupt:
        record=empty_record(sample,seed,version,'failed',stage,'AC-only execution interrupted',supported_labels=list(support),detected=None,**common)
        atomic_json(out/'stages.json',{'aa':{'status':'inherited'},'ac':{'status':'interrupted'}})
        _save(out,record,provenance)
        raise
    except (OSError,ValueError,KeyError,TypeError,RuntimeError,TimeoutError) as exc:
        status='missing' if isinstance(exc,FileNotFoundError) else 'failed' if common['ac_executed'] else 'unsupported'
        record=empty_record(sample,seed,version,status,stage,str(exc),supported_labels=list(support),detected=None,**common)
        atomic_json(out/'stages.json',{'aa':{'status':'inherited'},'ac':{'status':status,'reason':str(exc)}})
    return _save(out,record,provenance)


def _stable_hash(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def _resource_fingerprint(path):
    """Reference/software identity without reading a huge FASTA in full."""
    p=Path(path).resolve(strict=True)
    if p.suffix.lower() in ('.bam','.cram','.bai','.crai'):
        raise ValueError('Alignment files cannot be used as AC reference/cache resources')
    stat=p.stat();row={'path':str(p),'bytes':stat.st_size,'mtime_ns':stat.st_mtime_ns}
    if stat.st_size<=16*1024**2:
        row.update(sha256=digest(p),method='full_SHA256')
    else:
        h=hashlib.sha256()
        with p.open('rb') as handle:
            for offset in sorted({0,max(0,stat.st_size//2-32768),max(0,stat.st_size-65536)}):
                handle.seek(offset);h.update(str(offset).encode());h.update(handle.read(65536))
        row.update(sampled_sha256=h.hexdigest(),method='size_mtime_three_64KiB_blocks_not_full_integrity')
    return row


def _actual_reference_identity(config):
    ref=config['reference'];root=Path(ref['aa_data_repo'])/ref['build']
    if not root.is_dir():raise FileNotFoundError('Configured AC reference build directory missing')
    files={}
    def walk(directory,ancestors=()):
        actual=directory.resolve()
        if actual in ancestors:raise ValueError('Cyclic reference directory link')
        for p in sorted(directory.iterdir()):
            if p.name=='coverage.stats':continue
            if p.is_dir():walk(p,ancestors+(actual,))
            elif p.is_file():files[str(p.resolve())]=_resource_fingerprint(p)
    walk(root)
    if ref.get('fasta'):
        for p in (Path(ref['fasta']),Path(str(ref['fasta'])+'.fai')):
            files[str(p.resolve())]=_resource_fingerprint(p)
    if not files:raise ValueError('AC reference resources are empty')
    return sorted(files.values(),key=lambda r:r['path'])


def _source_aa_hashes(source_record):
    paths=_artifacts(Path(source_record['run_dir'])/'aa')
    expected=source_record.get('aa_artifact_sha256',{})
    raw={str(Path(p).resolve()):sha for p,sha in source_record.get('output_sha256',{}).items()}
    hashes={};total=0
    for p,_ in paths:
        if p.is_symlink() or not p.is_file():raise ValueError('Inherited AA must contain regular small artifacts')
        total+=p.stat().st_size
        if p.stat().st_size>MAX_FILE_BYTES or total>MAX_TOTAL_BYTES:raise ValueError('Inherited AA exceeds small-artifact limit')
        sha=digest(p)
        if sha!=expected.get(p.name,raw.get(str(p.resolve()))):raise ValueError('Inherited AA checksum missing or changed')
        hashes[p.name]=sha
    if not hashes:raise FileNotFoundError('Inherited AA artifacts missing')
    return hashes


def _raw_ac_cache_hashes(run):
    """Known raw scientific evidence only, never normalized result.json."""
    run=Path(run)
    paths=[p for p,_ in _artifacts(run/'aa',run/'ac')]
    for name in ('stages.json','execution.json','runtime_check.json','ac.input','ac_summary_map.txt','ac.log','ac/classification.log',
                 'provenance.json','reference_audit_before.json','reference_audit_after.json'):
        if (run/name).is_file():paths.append(run/name)
    result={}
    for p in paths:
        if p.is_symlink() or p.stat().st_size>MAX_FILE_BYTES:raise ValueError('Unsafe/oversized cached raw result')
        result[str(p.relative_to(run))]=digest(p)
    return result


def _next_ac_attempt(job):
    existing=[int(p.name[8:]) for p in Path(job).glob('attempt_*') if p.is_dir() and p.name[8:].isdigit()]
    return Path(job)/f'attempt_{max(existing,default=0)+1:04d}'


def _runtime_cache_identity(verified):
    return {'actual_ac':{k:verified['ac'].get(k) for k in ('version','commit','manifest_sha256','sha256','python')},
            'actual_packages':sorted(verified['ac'].get('packages',{}).get('stdout','').splitlines()),
            'actual_bfb':_probe_json(verified['bfb_probe']['stdout']) if verified.get('bfb_probe') else None}


def reclassify_cached(source_record,runtime,entry,config,job_dir):
    """Verify the actual AC runtime every call; reuse only hashed, reparsed output.

    A manifest cannot supply a trusted verification result. Signature fields are
    derived from local files and executable probes, not cached normalized JSON.
    """
    job=Path(job_dir).absolute();reference=config.get('reference',{})
    originals=[Path(source_record['run_dir']),Path(runtime['ac_script']).parent]
    if reference.get('aa_data_repo'):originals.append(Path(reference['aa_data_repo']))
    for p in originals:
        p=p.resolve()
        if p==job.resolve() or p in job.resolve().parents or job.resolve() in p.parents:
            raise ValueError('AC cache job must be separate from inputs/reference/source installation')
    job.mkdir(parents=True,exist_ok=True)
    verified=None;signature=None
    try:
        if runtime.get('id')!=entry.get('id'):raise ValueError('AC runtime and source-lock identifiers differ')
        if not reference.get('aa_data_repo') or not reference.get('build'):raise ValueError('Explicit AC reference resources required')
        env=_ac_environment(runtime)
        env.update(AA_DATA_REPO=str(Path(reference['aa_data_repo']).resolve()),MPLCONFIGDIR=str(job/'probe_mplcache'))
        # Required on every invocation, including a prospective cache hit.
        verified=_verify_ac_runtime(runtime,entry,env)
        hashes=_source_aa_hashes(source_record);graph_identity=_stable_hash(hashes)
        source_metadata=_source_metadata(source_record,graph_identity)
        if source_record.get('reference_build') and source_record['reference_build']!=reference['build']:
            raise ValueError('Source AA and target AC reference builds differ')
        resources=_actual_reference_identity(config)
        binaries={'python':_resource_fingerprint(Path(runtime['ac_python']).resolve())}
        cbc=shutil.which('cbc',path=env.get('PATH'))
        if cbc:binaries['cbc']=_resource_fingerprint(cbc)
        for item in runtime.get('dependency_overrides',[]):
            actual=_resource_fingerprint(item['path'])
            if item.get('sha256') and digest(item['path'])!=item['sha256']:raise ValueError('Declared dependency binary changed')
            binaries[item['path']]=actual
        # Store only a digest: environment values may contain credentials.
        effective_env={k:v for k,v in env.items() if k not in ('PWD','OLDPWD','SHLVL','_','MPLCONFIGDIR','AA_DATA_REPO')}
        identity={'source_aa_files':hashes,'source_metadata':source_metadata,'target_id':runtime['id'],
            'ac_lock':entry['ac'],'bfb_lock':entry.get('bfb'),'supported_labels':verified['supported_labels'],
            'ac_input_writer_sha256':digest(Path(__file__).with_name('runner.py')),
            **_runtime_cache_identity(verified),'binaries':binaries,'parameters':_ac_parameters(config),
            'reference_build':reference['build'],'reference_resources':resources,'environment_sha256':_stable_hash(effective_env)}
        signature=_stable_hash(identity)
        atomic_json(job/'preflight_latest.json',{'created_utc':now(),'signature':signature,'identity':identity,'verified':verified})
        cache_path=job/'cache.json'
        try:cache=json.loads(cache_path.read_text()) if cache_path.is_file() else {}
        except (OSError,ValueError):cache={}
        if not isinstance(cache,dict):cache={}
        run=Path(cache['run_dir']) if isinstance(cache.get('run_dir'),str) else Path('')
        valid_location=run.is_dir() and run.resolve().parent==job.resolve() and re.fullmatch(r'attempt_\d+',run.name)
        if cache.get('signature')==signature and valid_location:
            actual=_raw_ac_cache_hashes(run)
            if actual and actual==cache.get('raw_sha256') and _aa_identity(run)[0]==hashes:
                fresh=parse_run(run,source_record['sample_id'],source_record.get('seed_id',WHOLE_SAMPLE),runtime['id'],
                    verified['supported_labels'],source_record.get('reference_lengths') or config.get('_reference_lengths'))
                if fresh['status'] in ('ok','not_detected'):
                    aa_version=source_metadata['source_aa_version']
                    fresh.update(run_dir=str(run),analysis_kind='ac_only',provenance_kind='verified_ac_on_inherited_aa',
                        source_aa_version=aa_version,aa_version=aa_version,ac_version=entry['ac']['version'],
                        source_record_dir=str(Path(source_record['run_dir']).resolve()),source_metadata=source_metadata,
                        component_versions={'aa':aa_version,'ac':entry['ac']['version']},reference_build=reference['build'],
                        reference_scope=source_record.get('reference_scope','declared_reference_'+reference['build']),
                        completion_evidence='revalidated_hashed_AC_execution_and_inherited_AA',aa_executed=False,ac_executed=True,
                        ac_executed_this_invocation=False,ac_runtime_verified=True,source_graph_identity=graph_identity,
                        aa_artifact_sha256=hashes,cache_hit=True,cache_signature=signature,
                        comparison_eligible=bool(aa_version) and source_record.get('comparison_eligible') is True,
                        output_sha256={str(run/p):sha for p,sha in actual.items()},
                        comparison_limitations=['Revalidated cached AC-only classification; AA was not rerun. Source metadata remains declaration-scoped.'])
                    atomic_json(job/'result.json',fresh)
                    return fresh
    except (OSError,ValueError,KeyError,TypeError,RuntimeError,TimeoutError) as exc:
        attempt=_next_ac_attempt(job);attempt.mkdir()
        status='missing' if isinstance(exc,FileNotFoundError) else 'unsupported'
        record=empty_record(source_record['sample_id'],source_record.get('seed_id',WHOLE_SAMPLE),runtime.get('id',entry.get('id')),
            status,'ac',str(exc),analysis_kind='ac_only',provenance_kind='ac_only_attempt',comparison_eligible=False,
            source_metadata=_source_metadata(source_record),source_aa_version=_source_metadata(source_record)['source_aa_version'],
            aa_executed=False,ac_executed=False,ac_executed_this_invocation=False,ac_runtime_verified=verified is not None,
            supported_labels=entry.get('supported_labels',entry.get('supported_classes',[])),cache_hit=False,run_dir=str(attempt))
        atomic_json(attempt/'preflight_failure.json',{'created_utc':now(),'reason':str(exc),'verified':verified})
        atomic_json(attempt/'result.json',record);atomic_json(job/'result.json',record)
        return record
    attempt=_next_ac_attempt(job)
    result=reclassify_result(source_record,runtime,entry,config,attempt)
    result.update(cache_hit=False,cache_signature=signature,ac_executed_this_invocation=result.get('ac_executed',False))
    if result['status'] in ('ok','not_detected'):
        # The fresh execution probes again. Never label a changed backend or
        # changed reference resource with the earlier cache identity.
        try:
            actual_verified=json.loads((attempt/'runtime_check.json').read_text())
            result['cache_saved']=(_runtime_cache_identity(actual_verified)==_runtime_cache_identity(verified)
                                  and _actual_reference_identity(config)==resources)
        except (OSError,ValueError,KeyError,TypeError):result['cache_saved']=False
        if result['cache_saved']:
            atomic_json(job/'cache.json',{'signature':signature,'run_dir':str(attempt),'raw_sha256':_raw_ac_cache_hashes(attempt),
                'record':result,'warning':'Normalized record is for review only; cache reuse reparses raw artifacts and reconstructs metadata.'})
        else:result['cache_note']='Runtime/backend or reference changed between preflight and execution; this result was not cached'
    atomic_json(job/'result.json',result)
    return result
