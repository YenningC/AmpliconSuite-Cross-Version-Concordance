"""Exploratory seed-count sensitivity against a complete fixed-seed baseline.

This resamples cached independent seed reconstructions; it does not run AA,
recreate a full Suite workflow, or establish population-level adequacy.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

from .statistics import merge_intervals, subtract_intervals, _pair
from .reporting import write_tsv

UNION_METRICS = ("detection_loss", "category_loss", "interval_loss")
COMPARISON_METRICS = ("comparison_error", "comparison_event_loss")
COMPARISON_COMPONENTS = ("detection_difference_mean_absolute_deviation", "label_discordance_mean_absolute_deviation", "interval_distance_mean_absolute_deviation", "class_footprint_mean_absolute_deviation")
PAIRED_METRIC_SCALES=(('detection_difference',2.),('label_discordance',1.),('interval_distance',1.),('class_footprint_distance',1.))
METRICS = UNION_METRICS + COMPARISON_METRICS


def _comparison(rows):
    signal = [r for r in rows if r['signal_positive']]
    metrics = {}
    for name,_ in PAIRED_METRIC_SCALES:
        values = [r[name] for r in signal if r[name] is not None]
        metrics[name] = float(np.mean(values)) if values and len(values)==len(signal) else None
    events = set()
    if any(r['extra_detection'] for r in rows):events.add('extra_detection')
    if any(r['missed_detection'] for r in rows):events.add('missed_detection')
    if any(r['class_change'] for r in rows):events.add('class_change')
    if any(r.get('class_footprint_distance') is not None and r['class_footprint_distance']>0
           and not r['extra_detection'] and not r['missed_detection'] for r in rows):
        events.add('spatial_class_change')
    return metrics,events


def _comparison_loss(full, subset):
    metrics,events=full
    selected,selected_events=subset
    deviations={}
    no_signal=False
    for metric,scale in PAIRED_METRIC_SCALES:
        if metrics[metric] is None:
            deviations[metric]=None
        elif selected[metric] is None:
            # Losing all evaluable signal is maximal information loss, not zero discrepancy.
            deviations[metric]=scale
            no_signal=True
        else:
            deviations[metric]=abs(metrics[metric]-selected[metric])
    normalized=[deviations[m]/scale for m,scale in PAIRED_METRIC_SCALES if deviations[m] is not None]
    error=max(normalized) if len(normalized)==len(PAIRED_METRIC_SCALES) else np.nan
    lost=sorted(events-selected_events)
    event_loss=len(lost)/len(events) if events else np.nan
    return [error,event_loss],deviations,lost,no_signal


def _union(records):
    positive = [r for r in records if r["status"] == "ok" and r.get("detected", True)]
    labels = sorted({label for r in positive for label in r["labels"]})
    intervals = merge_intervals([i for r in positive for i in r["intervals"]])
    return bool(positive), labels, intervals


def _loss(full, selected):
    detected, labels, intervals = full
    sd, sl, si = selected
    lost_labels = sorted(set(labels)-set(sl))
    lost_intervals = subtract_intervals(intervals,si)
    bp = sum(e-s for _,s,e in intervals)
    lost_bp = sum(e-s for _,s,e in lost_intervals)
    values = [float(detected and not sd) if detected else np.nan,
              len(lost_labels)/len(labels) if labels else np.nan,
              lost_bp/bp if bp else np.nan]
    return values,lost_labels,lost_bp,lost_intervals


def _summary(values, versions, size, sample_indices, bootstrap, rng, metric_names=UNION_METRICS):
    """Array is repeat x sample x version x metric; resample samples as clusters."""
    rows = []
    for j,version in enumerate(versions):
        for m,metric in enumerate(metric_names):
            data = values[:,sample_indices,j,m]
            eligible = np.isfinite(data).all(axis=0)
            data = data[:,eligible]
            n = data.shape[1]
            row = dict(size=size,version=version,metric=metric,n_informative_samples=n,
                       mean_loss=None,min_subset_loss=None,max_subset_loss=None,
                       bootstrap_envelope_low=None,bootstrap_envelope_high=None,
                       bootstrap_degenerate=False)
            if n:
                averages = data.mean(axis=1)
                # Shared sample draws across repeats retain each sample's entire selected seed set.
                indices = rng.integers(0,n,size=(bootstrap,n))
                boot = data[:,indices].mean(axis=2)
                lo = float(np.quantile(boot,.025,axis=1).min())
                hi = float(np.quantile(boot,.975,axis=1).max())
                row.update(mean_loss=float(averages.mean()),min_subset_loss=float(averages.min()),
                           max_subset_loss=float(averages.max()),bootstrap_envelope_low=lo,
                           bootstrap_envelope_high=hi,bootstrap_degenerate=lo==hi)
            rows.append(row)
    return rows


def _decision(rows,tolerances):
    rows=[r for r in rows if r['metric'] in METRICS]
    if any(tolerances.get(m) is None for m in METRICS):
        return False,"tolerances_not_predefined"
    if any(r['n_informative_samples']==0 for r in rows if r['metric']!='comparison_event_loss'):
        return False,"no_informative_positive_baseline_for_some_versions_or_metrics"
    if not any(r['metric']=='comparison_event_loss' and r['n_informative_samples'] for r in rows):
        return False,"no_baseline_version_difference_events_to_assess_retention"
    rows=[r for r in rows if r['n_informative_samples']]
    if all(r['max_subset_loss'] <= tolerances[r['metric']] and
           r['bootstrap_envelope_high'] <= tolerances[r['metric']] for r in rows):
        return True,"all_sampled_subsets_and_bootstrap_envelopes_within_tolerances"
    return False,"loss_or_sample_bootstrap_envelope_exceeds_tolerance"


def _choose(rows, sizes, tolerances):
    # Require stability at this and every larger tested size, never choose a lone dip.
    passes = {k:_decision([r for r in rows if r['size']==k],tolerances)[0] for k in sizes}
    candidates = [k for i,k in enumerate(sizes) if all(passes[j] for j in sizes[i:])]
    return min(candidates) if candidates else None


def analyze_seed_sensitivity(records, options):
    """Return baseline readiness, repeated-subset losses and exploratory suggestion.

    Required: versions, reference_version, candidate_seeds={sample:[seed_ids]}, sizes, baseline_design.
    baseline_design must be independent_fixed_seed_runs. tolerances default empty.
    mode is random_nested (default) or prespecified_nested (requires seed_orders).
    """
    result = dict(ready=False,baseline_errors=[],summary_rows=[],subset_rows=[],comparison_subset_rows=[],size_rows=[],
                  leave_one_sample_out=[],recommendation=None,
                  inference="exploratory_only_not_population_validated",options=options)
    if options.get('baseline_design') != 'independent_fixed_seed_runs':
        result['baseline_errors'].append('Requires independent_fixed_seed_runs; joint/full-Suite catalog is not interchangeable')
        return result
    versions = list(options.get('versions',[]))
    result['reference_version']=options.get('reference_version')
    result['reference_scopes']=sorted({r.get('reference_scope','declared_reference') for r in records})
    candidates = options.get('candidate_seeds',{})
    if not versions or len(set(versions))!=len(versions) or not candidates or any(len(s)!=len(set(s)) for s in candidates.values()):
        result['baseline_errors'].append('Provide unique versions and explicit complete candidate_seeds for every sample')
        return result
    if options.get('reference_version') not in versions:
        result['baseline_errors'].append('Provide the fixed reference_version among versions')
        return result
    empty_samples=sorted(s for s in candidates if not candidates[s])
    result['selected_sample_count']=len(candidates)
    result['no_seed_samples']=empty_samples
    samples = sorted(s for s in candidates if candidates[s])
    for sample in empty_samples:
        sentinels=[r for r in records if r['sample_id']==sample]
        if len(sentinels)!=len(versions) or {r['version'] for r in sentinels}!=set(versions) or any(r['status']!='no_seeds' or r['seed_id']!='__no_seeds__' for r in sentinels):
            result['baseline_errors'].append(f'Explicit empty candidate sample requires one no_seeds sentinel per version: {sample}')
    records=[r for r in records if r['sample_id'] not in empty_samples]
    if not samples:
        result['baseline_errors'].append('No samples have input seeds; no seed-count sensitivity can be estimated')
    if result['baseline_errors']:
        return result
    expected = {(s,t,v) for s in samples for t in candidates[s] for v in versions}
    lookup = {}
    for r in records:
        key = (r['sample_id'],r['seed_id'],r['version'])
        if key in lookup:
            result['baseline_errors'].append(f'Duplicate baseline record: {key}')
        lookup[key] = r
    for key in sorted(expected - lookup.keys()):
        result['baseline_errors'].append(f'Missing baseline record: {key}')
    if lookup.keys() - expected:
        result['baseline_errors'].append('Baseline has records outside the explicitly supplied sample/seed/version universe')
    for key in sorted(expected & lookup.keys()):
        row = lookup[key]
        if row['status'] not in ('ok','not_detected'):
            result['baseline_errors'].append(f'Baseline not complete/evaluable: {key}: {row["status"]}')
        elif row['status']=='ok' and (not row.get('labels') or not row.get('intervals')):
            result['baseline_errors'].append(f'Positive baseline lacks labels or intervals: {key}')
        elif row['status']=='ok' and not row.get('detected',True):
            result['baseline_errors'].append(f'Contradictory positive status/detected flag: {key}')
        elif row['status']=='not_detected' and (row.get('labels') or row.get('intervals')):
            result['baseline_errors'].append(f'Negative baseline contains positive labels/intervals: {key}')
        if row.get('supported_labels') is None:
            result['baseline_errors'].append(f'Unknown category support: {key}')
        elif set(row.get('labels',[]))-set(row['supported_labels']):
            result['baseline_errors'].append(f'Labels outside supported categories: {key}')
    if result['baseline_errors']:
        return result
    sizes = sorted(set(int(x) for x in options.get('sizes',[])))
    cap = max(len(seeds) for seeds in candidates.values())
    if not sizes or sizes[0] < 1 or sizes[-1] > cap or sizes[-1] != cap:
        raise ValueError('sizes must be positive, <= maximum candidate count and include that full-baseline count')
    tolerances = options.get('tolerances',{})
    if any(k not in METRICS or v is not None and not 0 <= float(v) <= 1 for k,v in tolerances.items()):
        raise ValueError('tolerances must use known union/comparison loss names with values in [0,1]')
    tolerances = {k:float(v) if v is not None else None for k,v in tolerances.items()}
    mode = options.get('mode','random_nested')
    repeats = int(options.get('repeats',100)) if mode == 'random_nested' else 1
    bootstrap = int(options.get('bootstrap',1000))
    if mode not in ('random_nested','prespecified_nested') or repeats < 2 and mode=='random_nested' or bootstrap<100:
        raise ValueError('random_nested requires >=2 repeats; bootstrap >=100; mode must be supported')
    if mode=='prespecified_nested':
        orders = options.get('seed_orders',{})
        if any(s not in orders or len(orders[s])!=len(candidates[s]) or set(orders[s])!=set(candidates[s]) for s in samples):
            raise ValueError('prespecified_nested requires a complete seed_orders list per sample')
    rng = np.random.default_rng(int(options.get('seed',20261117)))
    full = {(s,v):_union([lookup[s,t,v] for t in candidates[s]]) for s in samples for v in versions}
    ref=options['reference_version']
    pairs=[(v,ref) for v in versions if v!=ref]
    pairs+=list(itertools.combinations([v for v in versions if v!=ref],2))
    pair_ids=[f'{a} vs {b}' for a,b in pairs]
    pair_lookup={(s,t,a,b):_pair(lookup[s,t,a],lookup[s,t,b],ref,'seed_sensitivity')
                 for s in samples for t in candidates[s] for a,b in pairs}
    result['comparison_scope_issues']=[dict(sample_id=s,seed_id=t,version=a,comparator=b,reason=row['class_footprint_reason'])
        for (s,t,a,b),row in pair_lookup.items() if row['signal_positive'] and row['class_footprint_distance'] is None]
    result['comparison_scope_complete']=not result['comparison_scope_issues']
    full_comparisons={(s,a,b):_comparison([pair_lookup[s,t,a,b] for t in candidates[s]])
                       for s in samples for a,b in pairs}
    orderings = [{s:list(rng.permutation(candidates[s])) if mode=='random_nested' else list(options['seed_orders'][s])
                  for s in samples} for _ in range(repeats)]
    arrays,comparison_arrays = {},{}
    for size in sizes:
        values = np.full((repeats,len(samples),len(versions),len(UNION_METRICS)),np.nan)
        comparisons=np.full((repeats,len(samples),len(pairs),len(COMPARISON_METRICS+COMPARISON_COMPONENTS)),np.nan)
        for repeat,ordering in enumerate(orderings):
            for i,sample in enumerate(samples):
                selected_seeds = ordering[sample][:size]
                for j,version in enumerate(versions):
                    subset = _union([lookup[sample,t,version] for t in selected_seeds])
                    losses,lost_labels,lost_bp,lost_intervals = _loss(full[sample,version],subset)
                    values[repeat,i,j] = losses
                    result['subset_rows'].append(dict(size=size,repeat=repeat+1,sample_id=sample,version=version,
                        selected_seed_count=len(selected_seeds),candidate_seed_count=len(candidates[sample]),
                        selected_seed_ids=selected_seeds,baseline_detected=full[sample,version][0],selected_detected=subset[0],
                        lost_categories=lost_labels,lost_interval_bp=lost_bp,lost_intervals=lost_intervals,
                        **{m:float(x) if np.isfinite(x) else None for m,x in zip(UNION_METRICS,losses)}))
                for j,(a,b) in enumerate(pairs):
                    full_pair=full_comparisons[sample,a,b]
                    sub_pair=_comparison([pair_lookup[sample,t,a,b] for t in selected_seeds])
                    losses,deviations,lost_events,no_signal=_comparison_loss(full_pair,sub_pair)
                    components=[deviations[k] if deviations[k] is not None else np.nan for k,_ in PAIRED_METRIC_SCALES]
                    comparisons[repeat,i,j]=losses+components
                    result['comparison_subset_rows'].append(dict(size=size,repeat=repeat+1,sample_id=sample,
                        version=a,comparator=b,selected_seed_ids=selected_seeds,
                        baseline_signal_regions=sum(pair_lookup[sample,t,a,b]['signal_positive'] for t in candidates[sample]),
                        selected_signal_regions=sum(pair_lookup[sample,t,a,b]['signal_positive'] for t in selected_seeds),
                        baseline_metrics=full_pair[0],selected_metrics=sub_pair[0],absolute_metric_deviations=deviations,
                        signed_detection_difference_shift=(sub_pair[0]['detection_difference']-full_pair[0]['detection_difference'])
                            if sub_pair[0]['detection_difference'] is not None and full_pair[0]['detection_difference'] is not None else None,
                        baseline_difference_events=sorted(full_pair[1]),lost_difference_events=lost_events,
                        selected_signal_missing=no_signal,
                        **{m:float(x) if np.isfinite(x) else None for m,x in zip(COMPARISON_METRICS,losses)}))
        arrays[size]=values
        comparison_arrays[size]=comparisons
        result['summary_rows'] += _summary(values,versions,size,list(range(len(samples))),bootstrap,rng)
        result['summary_rows'] += _summary(comparisons,pair_ids,size,list(range(len(samples))),bootstrap,rng,COMPARISON_METRICS+COMPARISON_COMPONENTS)
    for size in sizes:
        rows = [r for r in result['summary_rows'] if r['size']==size]
        passed,reason = _decision(rows,tolerances)
        if not result['comparison_scope_complete']:
            passed,reason=False,'class_footprint_information_missing_or_unavailable'
        result['size_rows'].append(dict(size=size,observed_criteria_passed=passed,reason=reason,
            n_samples=len(samples),repeats=repeats,min_selected_per_sample=min(min(size,len(candidates[s])) for s in samples),
            max_selected_per_sample=max(min(size,len(candidates[s])) for s in samples),
            min_unique_subsets_per_sample=min(len({tuple(sorted(o[s][:size])) for o in orderings}) for s in samples)))
    minimum = _choose(result['summary_rows'],sizes,tolerances) if result['comparison_scope_complete'] else None
    loo_minima = []
    if len(samples)>=2:
        for omitted in range(len(samples)):
            rows = []
            for size in sizes:
                rows += _summary(arrays[size],versions,size,[i for i in range(len(samples)) if i!=omitted],bootstrap,rng)
                rows += _summary(comparison_arrays[size],pair_ids,size,[i for i in range(len(samples)) if i!=omitted],bootstrap,rng,COMPARISON_METRICS+COMPARISON_COMPONENTS)
            chosen = _choose(rows,sizes,tolerances)
            loo_minima.append(chosen)
            result['leave_one_sample_out'].append(dict(omitted_sample=samples[omitted],minimum_observed_size=chosen,
                reason='descriptive_selection_stability_not_external_validation',training_samples=len(samples)-1))
    else:
        result['leave_one_sample_out'].append(dict(omitted_sample=samples[0],minimum_observed_size=None,
            reason='insufficient_samples_for_leave_one_out',training_samples=0))
    result.update(ready=True,n_samples=len(samples),versions=versions,sizes=sizes,repeats=repeats,
                  baseline_candidate_total=sum(map(len,candidates.values())),minimum_observed_size=minimum)
    if minimum is not None and mode=='random_nested' and repeats>=20 and len(samples)>=2 and all(x is not None for x in loo_minima):
        stable_minimum=max([minimum]+loo_minima)
        larger=[k for k in sizes if k>stable_minimum]
        recommended=min(larger) if larger else cap
        result['recommendation']=dict(seed_count_per_sample_cap=recommended,
            minimum_observed_size=minimum,maximum_leave_one_out_minimum=max(loo_minima),
            safety_rule='next larger tested size above full-data and leave-one-out minima, capped by available candidates',
            at_candidate_cap=recommended==cap,
            interpretation='exploratory suggestion for these samples and candidate seeds; not a sufficient number for general studies')
    else:
        result['recommendation_reason']=('class_footprint_information_missing_or_unavailable' if not result['comparison_scope_complete'] else
            'tolerances/multisample positive evidence/repeated-subset stability insufficient; no adequacy claim')
    for table in ('summary_rows','subset_rows','comparison_subset_rows','size_rows','leave_one_sample_out'):
        for row in result[table]:
            row['reference_version']=result['reference_version']
    return result


def write_sensitivity_report(outdir, result):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    (out/'seed_sensitivity.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    for key,name in [('summary_rows','size_metrics.tsv'),('subset_rows','subset_losses.tsv'),
                     ('comparison_subset_rows','comparison_retention.tsv'),
                     ('size_rows','size_decisions.tsv'),('leave_one_sample_out','leave_one_sample_out.tsv')]:
        write_tsv(out/name,result[key])
    lines=['# 探索性种子数量敏感性分析','',
      '**这不是充分种子数或总体兼容性验证。** 仅重采样既有独立固定种子重建的并集，不重新运行 AA，不验证全流程 CNVkit 或联合种子运行的交互。',
      '',f"基线可分析：{result['ready']}。",'']
    lines += [f"固定参考版本：`{result.get('reference_version')}`；所含版本均分别对自己的完整固定候选基线评估保留能力。",'',
              '参考范围：`'+','.join(result.get('reference_scopes',[]))+'`。若仅人类参考可用，不能外推病毒相关结构或完整原始参考。','']
    if not result['ready']:
        lines+=['基线失败、缺失、类别不明确或候选清单不完整时停止数量判断：','']+[f'- {e}' for e in result['baseline_errors']]
    else:
        lines += [f"有输入种子的独立样本数：{result['n_samples']}；原选择样本数：{result.get('selected_sample_count',result['n_samples'])}；无种子样本数：{len(result.get('no_seed_samples',[]))}；基线候选总数：{result['baseline_candidate_total']}；每档随机嵌套子集重复：{result['repeats']}。",
         '', '各重复中同一样本使用相同种子集合比较全部版本。先合并同样本、同版本的结构区段与类别，再对有信号的样本等权统计。全阴性基线不提供阳性保留能力证据。',
         '', '表中 bootstrap envelope 是所有已抽取子集的样本整簇95%百分位区间包络，不是同时置信区间或新样本预测区间。5例时尤其只能描述当前样本；零损失的退化区间不能证明未见样本也无损失。',
         '', '版本对保留指标还检查：正信号区域内配对检出差、标签失配、总体区段距离和分类特异区段距离的样本内均值偏离；以及基线增检/漏检/分类改变/分类位置改变事件类型是否丢失。任何应有feature位置资料缺失都阻止数量推荐，不能只保留已知部分。仅保留每个版本的结构并集不足以通过。没有基线版本差异事件时，不宣称已验证差异事件保留能力。',
         '', '容忍度：`'+json.dumps(result['options'].get('tolerances',{}),ensure_ascii=False)+'`。必须预设 detection_loss、category_loss、interval_loss、comparison_error、comparison_event_loss 五项容忍度才可能推荐数量。', '']
        if result['recommendation']:
            r=result['recommendation']
            lines += [f"**探索性建议：每样本最多 {r['seed_count_per_sample_cap']} 个候选种子**（候选不足时取全部）。当前样本首次满足阈值的档位为 {r['minimum_observed_size']}；留一例选择中最大的最低档位为 {r['maximum_leave_one_out_minimum']}。建议取更高一档并受候选总数约束。",'',
                '这只是固定候选集合内的保守执行数量，不称真实充分数；推荐位于候选上限时，无法额外增加余量。']
        else:
            lines += ['**证据不足：不给出充分性或推荐数量结论。** 可检查预设容忍度、完整基线与有信号样本数；不可通过删去差异种子来提高通过率。']
    (out/'README.md').write_text('\n'.join(lines)+'\n')
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records',required=True)
    parser.add_argument('--options',required=True)
    parser.add_argument('--out',required=True)
    args=parser.parse_args(argv)
    records=json.loads(Path(args.records).read_text())
    options=json.loads(Path(args.options).read_text())
    result=analyze_seed_sensitivity(records,options)
    result['input_provenance']=[dict(path=str(Path(path).resolve()),sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest()) for path in (args.records,args.options)]
    write_sensitivity_report(args.out,result)
    print(json.dumps({k:result.get(k) for k in ['ready','n_samples','minimum_observed_size','recommendation','recommendation_reason']},ensure_ascii=False))
    return 0 if result['ready'] else 2


if __name__=='__main__':
    raise SystemExit(main())
