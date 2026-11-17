"""Optional local adapters to unchanged official comparison tools.

These tools describe feature/path differences. Their scores never establish
biological truth or statistical compatibility. Always retain primary paired
status, cluster inference, and prespecified margins from ascompare.statistics.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path
from .full_suite import WHOLE_SAMPLE

from .common import atomic_json, digest, now, table
from .inputs import safe_id
from .runner import execute
from .versions import load_lock, _check_source


def _output(directory, inputs):
    out = Path(directory).absolute()
    for source in inputs:
        p = Path(source).resolve()
        if out.resolve() == p or out.resolve() in p.parents or p.is_dir() and p in out.resolve().parents:
            raise ValueError('Diagnostic output must not contain/replace source inputs')
    if out.exists() and any(out.iterdir()):
        raise ValueError('Use an empty output directory for official diagnostics')
    out.mkdir(parents=True, exist_ok=True)
    return out


def _runtime(python, env, out):
    """Record actual interpreter, dependency versions, and external bedtools."""
    probe=subprocess.run([str(Path(python).absolute()),'-c',
        'import sys,json,importlib.metadata as m;print(json.dumps({"executable":sys.executable,"prefix":sys.prefix,"python":sys.version,"packages":{d.metadata["Name"]:d.version for d in m.distributions()}}))'],
        env=env,capture_output=True,text=True,timeout=60)
    (out/'runtime.log').write_text(probe.stdout+probe.stderr)
    if probe.returncode:
        raise RuntimeError('Cannot inspect official-tool Python environment')
    runtime=json.loads(probe.stdout)
    bedtools=shutil.which('bedtools',path=env.get('PATH'))
    runtime['bedtools']={'path':bedtools,'sha256':digest(bedtools)} if bedtools else None
    atomic_json(out/'runtime.json',runtime)
    return runtime


def valid_ec_metrics(data):
    values=data.get('distances',{}) if isinstance(data,dict) else {}
    if not isinstance(values,dict):return False
    required=('cn_hamming_dist','cn_hamming_norm_dist','cn_cos_dist','cn_jc_dist',
              'breakpoint_dist','fragments_dist','cycles_dist','total_cost')
    # total_cost_description is an upstream dictionary of weights, not a score.
    return all(isinstance(values.get(key),(float,int)) and math.isfinite(values[key]) for key in required)


def export_cycles(paths, destination):
    """Transport AA cycles into official EC BED schema, preserving path semantics.

    Unlike upstream cycles_to_bed_modified.py, use explicit iscyclic=False for
    source/sink paths and convert inclusive AA right endpoints to BED end+1.
    Repeated segments and traversal order are retained, not merged or sorted.
    """
    fields = ['#chr','start','end','strand','circ_id','estimated_cn','iscyclic']
    rows, sources = [], []
    for file_number, path in enumerate(paths, 1):
        path = Path(path).resolve()
        segments, cycles = {}, []
        for line in path.read_text().splitlines():
            if line.startswith('Segment'):
                f = line.split()
                sid, start, end = int(f[1]), int(f[3]), int(f[4])
                if sid in segments or start < 0 or end < start:
                    raise ValueError('Invalid/duplicate AA segment')
                segments[sid] = (f[2], start, end+1)
            elif line.startswith('Cycle='):
                info = dict(item.split('=',1) for item in line.split(';') if '=' in item)
                cycles.append(info)
        seen = set()
        for info in cycles:
            cycle_id = info['Cycle']
            if cycle_id in seen:
                raise ValueError('Duplicate cycle ID in input')
            seen.add(cycle_id)
            tokens = info['Segments'].split(',')
            if any(not re.fullmatch(r'\d+[+-]', token) for token in tokens):
                raise ValueError('Unrecognized AA path token')
            source_sink = [i for i,t in enumerate(tokens) if int(t[:-1]) == 0]
            if source_sink and source_sink != [0, len(tokens)-1]:
                raise ValueError('AA source/sink not at both terminal positions; path conversion unsupported')
            cyclic = not source_sink
            cn = float(info['Copy_count'])
            if not math.isfinite(cn) or cn < 0:
                raise ValueError('Invalid AA cycle copy count')
            for token in tokens:
                sid = int(token[:-1])
                if not sid:
                    continue
                if sid not in segments:
                    raise ValueError('AA cycle references missing segment')
                chrom, start, end = segments[sid]
                rows.append(dict(zip(fields,[chrom,start,end,token[-1],f'{file_number}_{cycle_id}',cn,cyclic])))
        sources.append({'path':str(path),'sha256':digest(path)})
    if not rows:
        raise ValueError('No reconstructed segments: EC structural diagnostic is not evaluable')
    with Path(destination).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
    return {'input_files':sources,'amplicon_file_count':len(sources),'inter_amplicon_edges_added':False,
            'collection_scope':'Independent reconstructed paths/cycles; not one assembled physical molecule; overlapping path CN contributions retained',
            'transport':'AA inclusive end +1; explicit iscyclic; ordered traversals including repeats; unique file-prefixed cycle IDs',
            'rows':len(rows),'structures':len({r['circ_id'] for r in rows})}


def _validate_results_table(path):
    with open(path,newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not {'Sample name','Feature ID','Classification','Location'} <= set(reader.fieldnames or []):
            raise ValueError('Official result table required, not classification_profiles.tsv')
        rows=list(reader)
    if not rows:
        raise ValueError('Empty result table is not evidence of successful negative analysis')
    for row in rows:
        location=(row.get('Location') or '').strip()
        feature_id=row.get('Feature ID') or ''
        if location in ('','[]',"['']") and feature_id.endswith('_NA'):
            continue
        # The native tool silently drops malformed coordinates, including just
        # one invalid member of an otherwise valid list. Reject partial reads.
        tokens=location.removeprefix('[').removesuffix(']').split(',')
        for token in tokens:
            match=re.fullmatch(r'([^\s:\[\]\'\"]+):(\d+)-(\d+)',token.strip().strip('\'\"'))
            if not match or int(match[3])<=int(match[2]):
                raise ValueError('Feature has missing/invalid coordinates; upstream could ignore it')
        if '_NA' in feature_id:
            raise ValueError('Coordinate-bearing feature ID contains reserved _NA; upstream would ignore it')
    return rows


def results_comparator(table1, table2, source, python, out, reference_version, candidate_version, version_lock=None):
    """Run unchanged official results_comparator: run1=reference, run2=candidate.

    Inputs must be completed native AC result tables from the same sample(s)
    and seed(s). This function cannot infer execution success from empty rows.
    """
    lock=load_lock(version_lock)
    entry=next(v for v in lock['versions'] if v['id']==lock['reference_version'])
    verified=_check_source(source,entry['ac'])
    source=Path(source).resolve();script=source/'scripts/results_comparator.py'
    if not script.is_file():
        raise ValueError('Pinned AC source has no official results_comparator.py')
    rows1,rows2=_validate_results_table(table1),_validate_results_table(table2)
    if {r['Sample name'] for r in rows1}!={r['Sample name'] for r in rows2}:
        raise ValueError('Official comparison requires matching analyzed sample identities')
    out=_output(out,[table1,table2,source]);local=[]
    for index,path in enumerate((table1,table2),1):
        target=out/f'run{index}_native_result_table.tsv';shutil.copyfile(path,target);local.append(target)
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    runtime=_runtime(python,env,out)
    command=[str(Path(python).absolute()),str(script),*map(str,local),'-o',str(out/'official')]
    metadata={'created_utc':now(),'tool':'AmpliconClassifier/scripts/results_comparator.py','source':verified,
              'tool_version':entry['ac']['version'],'reference_version':reference_version,'candidate_version':candidate_version,
              'orientation':'run1 reference; run2 candidate','inputs':[{'path':str(Path(p).resolve()),'sha256':digest(p)} for p in (table1,table2)],
              'argv':command,'scope':'feature-overlap diagnostic; >=1000bp upstream matching; no compatibility decision',
              'status':'running','runtime':runtime}
    atomic_json(out/'provenance.json',metadata)
    try:
        rc=execute(command,out,env,out/'official.log',timeout=600)
    except (OSError,TimeoutError) as exc:
        metadata.update(status='failed',reason=str(exc));atomic_json(out/'provenance.json',metadata);return metadata
    required=[out/f'official_{x}' for x in ('comparison_report.txt','comparison_table.csv','classification_summary.txt')]
    metadata.update(returncode=rc,status='complete' if rc==0 and all(p.is_file() for p in required) else 'failed')
    if metadata['status']=='failed':metadata['reason']='Official command failed or required raw outputs are absent; inspect official.log'
    metadata['output_sha256']={p.name:digest(p) for p in required if p.is_file()}
    atomic_json(out/'provenance.json',metadata)
    return metadata


def ec_comparator(cycles1,cycles2,source,python,out,reference_version,candidate_version,extra_path=()):
    """Run source-pinned official EC comparator on all AA paths for one seed pair."""
    lock=json.loads(Path(__file__).with_name('official_tools.lock.json').read_text())['ec_comparator']
    source=Path(source).resolve()
    for name,sha in lock['files'].items():
        if not (source/name).is_file() or digest(source/name)!=sha:
            raise ValueError('Official EC-comparator source differs from pinned commit: '+name)
    extra={str(p.relative_to(source)) for p in source.rglob('*.py')}-set(lock['files'])
    if extra:
        raise ValueError('Unexpected Python files in pinned EC-comparator source: '+','.join(sorted(extra)[:3]))
    out=_output(out,[*cycles1,*cycles2,source])
    one,two=out/'reference_paths.bed',out/'candidate_paths.bed'
    meta={'created_utc':now(),'tool':'EC-comparator','tool_version':lock['version'],'commit':lock['commit'],
          'repository':lock['repository'],'reference_version':reference_version,'candidate_version':candidate_version,
          'scope':'Collection of independent AA reconstructed paths/cycles, not one physical molecule; overlapping CN contributions retained; not confirmed ecDNA or compatibility probability',
          'reference_transport':export_cycles(cycles1,one),'candidate_transport':export_cycles(cycles2,two)}
    env=os.environ.copy();env.update(PYTHONPATH=str(source),PYTHONNOUSERSITE='1',PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(out/'mplcache'))
    env['PATH']=os.pathsep.join([str(Path(python).absolute().parent),*map(str,extra_path),env.get('PATH','')])
    meta['runtime']=_runtime(python,env,out)
    mismatches=[]
    for line in (source/'requirements.txt').read_text().splitlines():
        if '==' in line:
            name,expected=line.strip().split('==')
            installed={k.lower().replace('_','-'):v for k,v in meta['runtime']['packages'].items()}.get(name.lower().replace('_','-'))
            if installed!=expected:
                mismatches.append(f'{name}: expected {expected}, found {installed}')
    if mismatches or not meta['runtime']['bedtools']:
        meta.update(status='unsupported',reason='Pinned optional dependencies unavailable',dependency_mismatches=mismatches)
        atomic_json(out/'provenance.json',meta)
        return meta
    raw=out/'official_raw'
    command=[str(Path(python).absolute()),'-m','eccomparator.main','-a',str(one),'-b',str(two),'-d',str(raw),
             '--alabel',reference_version,'--blabel',candidate_version,'--no-plot','--no-report','--min-cn','0',
             '--breakpoint-selection-threshold','1000','--breakpoint-cost-function-threshold','3000',
             '--breakpoint-dist-calc','breakpoint_match_unweighted','--gap','1000000']
    meta.update(argv=command,status='running');atomic_json(out/'provenance.json',meta)
    try:
        rc=execute(command,out,env,out/'official.log',timeout=600)
    except (OSError,TimeoutError) as exc:
        meta.update(status='failed',reason=str(exc));atomic_json(out/'provenance.json',meta);return meta
    metrics=raw/'metrics.json'
    # Official main catches exceptions without nonzero exit: require actual metrics.
    try:
        data=json.loads(metrics.read_text())
        valid=valid_ec_metrics(data)
    except (OSError,ValueError):
        valid=False
    meta.update(returncode=rc,status='complete' if rc==0 and valid else 'failed')
    if meta['status']=='failed':meta['reason']='Official execution failed or finite metrics.json unavailable, including upstream exceptions with exit0; inspect official.log'
    meta['output_sha256']={str(p.relative_to(out)):digest(p) for p in raw.rglob('*') if p.is_file()} if raw.exists() else {}
    atomic_json(out/'provenance.json',meta)
    return meta


def compare_completed_records(records, out, toolconfig, reference_version):
    """Optional diagnostics on matched completed sample/seed records only.

    toolconfig={'results':{'source':..., 'python':...}, 'ec':{...}}.
    Missing native tables, incomplete runs and absent cycles remain explicit
    entries in index.tsv; never create fake empty native result tables.
    Whole-sample comparisons support independent amplicon structure sets.
    Seed-specific multi-amplicon aggregation requires an explicit EC opt-in.
    """
    root=Path(out)/'official_comparisons';root.mkdir(parents=True,exist_ok=True)
    lookup={(r['version'],r['sample_id'],r['seed_id']):r for r in records}
    versions=sorted({r['version'] for r in records}-{reference_version})
    keys=sorted({(r['sample_id'],r['seed_id']) for r in records})
    index=[]
    for candidate in versions:
        for sample,seed in keys:
            a=lookup.get((reference_version,sample,seed));b=lookup.get((candidate,sample,seed))
            for tool in ('results','ec'):
                options=toolconfig.get(tool)
                if not options:
                    continue
                row={'reference_version':reference_version,'candidate_version':candidate,'sample_id':sample,
                     'seed_id':seed,'tool':tool,'status':'not_evaluable','reason':'','diagnostic_only':True}
                index.append(row)
                if not a or not b or a['status'] not in ('ok','not_detected') or b['status'] not in ('ok','not_detected'):
                    row['reason']='Missing or unevaluable primary run; no negative call inferred';continue
                if not a.get('run_dir') or not b.get('run_dir'):
                    row.update(status='missing',reason='Native run directory unavailable');continue
                ar,br=Path(a['run_dir']),Path(b['run_dir'])
                if tool=='results':
                    first=list((ar/'ac').glob('*_result_table.tsv'));second=list((br/'ac').glob('*_result_table.tsv'))
                    if len(first)!=1 or len(second)!=1:
                        row.update(status='missing',reason='Native official result table absent/ambiguous in one or both versions; no fake table generated');continue
                else:
                    first=sorted((ar/'aa').glob('*_amplicon*_cycles.txt'));second=sorted((br/'aa').glob('*_amplicon*_cycles.txt'))
                    if not first or not second:
                        row.update(status='not_evaluable',reason='One or both versions have no reconstructed cycle/path input');continue
                    if (len(first)>1 or len(second)>1) and not options.get('allow_multiple_amplicons',seed in (WHOLE_SAMPLE,'whole_sample')):
                        row.update(status='unsupported',reason='Multiple amplicons require explicit aggregate structure-set opt-in');continue
                try:
                    seed_directory='whole_sample' if seed in (WHOLE_SAMPLE,'whole_sample') else safe_id(seed)
                    base=root/safe_id(sample)/seed_directory/safe_id(candidate)/tool
                    attempt=base/f'attempt_{len(list(base.glob("attempt_*")))+1:04d}'
                    kwargs=dict(source=options['source'],python=options['python'],out=attempt,
                                reference_version=reference_version,candidate_version=candidate)
                    result=(results_comparator(first[0],second[0],**kwargs) if tool=='results' else
                            ec_comparator(first,second,extra_path=options.get('extra_path',[]),**kwargs))
                    row.update(status=result['status'],reason=result.get('reason','Official raw output retained; not a compatibility decision'),
                               directory=str(attempt.resolve()),provenance=str((attempt/'provenance.json').resolve()))
                except Exception as exc:
                    row.update(status='failed',reason=f'{type(exc).__name__}: {exc}')
    table(root/'index.tsv',index)
    atomic_json(root/'index.json',index)
    summarize_official_metrics(index,root)
    # Plot from the retained native scalar diagnostics; a plotting failure must
    # not alter the execution status or imply a compatibility decision.
    # A failed refresh must not link older figures as current observations.
    table(root/'official_plots.tsv',[],['sample_id','seed_id','reference_version','figure','path'])
    atomic_json(root/'plots_status.json',{'status':'running','diagnostic_only':True})
    try:
        from .official_plots import plot_official_metrics
        figures=plot_official_metrics(out,records=records)
        atomic_json(root/'plots_status.json',{'status':'complete','figure_count':len(figures),
                                             'diagnostic_only':True})
    except Exception as exc:
        atomic_json(root/'plots_status.json',{'status':'failed','reason':f'{type(exc).__name__}: {exc}',
                                             'diagnostic_only':True})
    return index


def summarize_official_metrics(index, root):
    """Scalar exports retain native direction, units and unavailable statuses."""
    output=[]
    for pair in index:
        base={k:pair.get(k,'') for k in ('reference_version','candidate_version','sample_id','seed_id','tool','status','reason','directory')}
        base.update(diagnostic_only=True,compatibility_cutoff=None,compatibility_cutoff_basis='No cross-version acceptance cutoff identified in inspected pinned official documentation')
        if pair['status']!='complete':
            output.append({**base,'metric':'','value':None,'unit':'','direction':''});continue
        directory=Path(pair['directory'])
        try:
            provenance=json.loads((directory/'provenance.json').read_text())
            base['tool_version']=provenance.get('tool_version','')
            if pair['tool']=='ec':
                metrics=json.loads((directory/'official_raw/metrics.json').read_text())['distances']
                for metric,value in metrics.items():
                    if not isinstance(value,(int,float)):
                        continue
                    lo,hi,unit=(0,None,'bp') if metric=='cn_hamming_dist' else (0,5,'weighted_cost') if metric=='total_cost' else (0,1,'distance')
                    within=math.isfinite(value) and value>=lo and (hi is None or value<=hi)
                    output.append({**base,'metric':metric,'value':value,'unit':unit,'direction':'lower_is_closer',
                                   'documented_min':lo,'documented_max':hi,'within_documented_range':within,
                                   'range_source':'EC-comparator v0.1.3 README Distances explained; total cost for this fixed all-enabled configuration',
                                   'status':'complete' if within else 'out_of_documented_range'})
            else:
                with (directory/'official_comparison_table.csv').open(newline='') as f:
                    events=list(csv.DictReader(f))
                counts=Counter(r['Event_Type'] for r in events)
                for event in ('GAIN','LOSS','CLASS_CHANGE','SPLIT','MERGE','COMPLEX_MERGE_SPLIT','STABLE'):
                    output.append({**base,'metric':event+'_count','value':counts[event],'unit':'overlap_components',
                                   'direction':'descriptive_count','interpretation':'Official overlap relationship; overlapping multi-label features may produce merge/split on unchanged inputs'})
                for event in events:
                    try:value=float(event['Jaccard_Overlap'])
                    except (TypeError,ValueError):continue
                    output.append({**base,'metric':'component_interval_jaccard','value':value,'unit':'fraction',
                                   'direction':'higher_is_closer','documented_min':0,'documented_max':1,
                                   'reference_features':event['Run1_Features'],'candidate_features':event['Run2_Features'],
                                   'official_event':event['Event_Type'],'interpretation':'Component-level interval overlap, not sample-wide structural truth'})
        except (OSError,ValueError,KeyError) as exc:
            output.append({**base,'metric':'','value':None,'status':'failed','reason':str(exc),'unit':'','direction':''})
    table(Path(root)/'official_metrics.tsv',output)
    atomic_json(Path(root)/'official_metrics.json',output)
    return output


def main():
    parser=argparse.ArgumentParser(description='Optional unchanged official comparison diagnostics; local only')
    sub=parser.add_subparsers(dest='tool',required=True)
    for tool in ('results','ec'):
        p=sub.add_parser(tool);p.add_argument('--source',required=True);p.add_argument('--python',required=True)
        p.add_argument('--out',required=True);p.add_argument('--reference-version',required=True);p.add_argument('--candidate-version',required=True)
        p.add_argument('--first',required=True,nargs='+' if tool=='ec' else None);p.add_argument('--second',required=True,nargs='+' if tool=='ec' else None)
        if tool=='ec': p.add_argument('--extra-path',action='append',default=[])
    args=parser.parse_args()
    kwargs=dict(source=args.source,python=args.python,out=args.out,reference_version=args.reference_version,candidate_version=args.candidate_version)
    result=results_comparator(args.first,args.second,**kwargs) if args.tool=='results' else ec_comparator(args.first,args.second,extra_path=args.extra_path,**kwargs)
    print(json.dumps(result,ensure_ascii=False,indent=2));return 0 if result['status']=='complete' else 2


if __name__=='__main__':
    raise SystemExit(main())
