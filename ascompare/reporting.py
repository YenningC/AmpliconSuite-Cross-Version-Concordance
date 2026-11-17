"""Local, static scientific reports; no network calls or patient data export."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import tempfile
import re
from collections import Counter
from pathlib import Path


def _cell(value):
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return value


def write_tsv(path, rows):
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows({key: _cell(value) for key,value in row.items()} for row in rows)


def _plots(out, records, stats, compatibility, ref, synthetic, analysis_scope='fixed_seed'):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "ascompare-matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import numpy as np
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":9, "axes.spines.top":False,
                         "axes.spines.right":False, "pdf.fonttype":42, "svg.fonttype":"none"})
    versions = list(dict.fromkeys([ref]+[r["version"] for r in records]+[r["version"] for r in compatibility]))
    others = [v for v in versions if v != ref]
    title = f"Fixed reference: {ref}" + ("  |  SYNTHETIC TEST DATA" if synthetic else "")
    title += '\nFull Suite, including CNVkit; whole-sample endpoints' if analysis_scope=='full_suite' else '\nFixed seeds: AA + AC; CNVkit outside scope'
    incomplete_reference = any(r.get('reference_scope') == 'human_only_incomplete_reference' for r in records)
    if incomplete_reference:
        title += "\nScope: human reference only; viral contigs unavailable"
    metrics = [("detection_agreement", "Detection agreement"), ("label_discordance", "Label-set agreement"),
               ("interval_jaccard", "Interval Jaccard")]
    fig, axes = plt.subplots(1,3, figsize=(12,max(4.0,.43*len(versions)+2.4)), sharey=True)
    for ax,(metric,label) in zip(axes,metrics):
        for y,version in enumerate(versions):
            if version == ref:
                ax.scatter([1],[y],marker="D",s=30,color="#333333",zorder=4)
                continue
            row = next((r for r in stats if r["version"] == version and r["comparator"] == ref
                        and r["subset"] == "signal_positive" and r["metric"] == metric), None)
            if not row or row["estimate"] is None:
                ax.text(.5,y,"not evaluable",color="#777777",ha="center",va="center",fontsize=8)
                continue
            value,lo,hi = row["estimate"],row["conservative_ci_low"],row["conservative_ci_high"]
            if metric == "label_discordance":
                value,lo,hi = 1-value,1-hi,1-lo
            ax.errorbar(value,y,xerr=[[max(0,value-lo)],[max(0,hi-value)]],fmt="o",capsize=3,
                        color="#4178A9",elinewidth=1.1,markersize=5)
            counts=f"n={row['n_samples']} samples" if analysis_scope=='full_suite' else f"n={row['n_samples']}; regions={row['n_regions']}"
            ax.text(.02,y-.22,counts,fontsize=6.5,color="#666666")
        ax.set_xlim(-.03,1.05)
        ax.set_xticks([0,.25,.5,.75,1], ["0","25","50","75","100"])
        ax.set_xlabel(label+" (%)")
        ax.grid(axis="x",color="#e7e7e7",linewidth=.6)
        ax.set_axisbelow(True)
    axes[0].set_yticks(range(len(versions)), [v+(" [reference]" if v == ref else "") for v in versions])
    # Text-only unavailable rows do not contribute to autoscaling; reserve a
    # full half-row at both ends even when only the reference has a point.
    axes[0].set_ylim(len(versions)-.5,-.5)
    unit='samples' if analysis_scope=='full_suite' else 'regions'
    fig.suptitle(f"Version ladder — {unit} with amplification signal in either version\n"+title,fontsize=11)
    confidence = stats[0]['confidence_level'] if stats else .95
    fig.text(.5,.015,f"Points: equal-weight sample means. Bars: conservative {100*confidence:g}% bounds. Reference is an identity anchor.\nLabel comparison may be restricted to shared supported categories; see differences.tsv. No biological truth probability is estimated.",ha="center",fontsize=7)
    fig.tight_layout(rect=[0,.08,1,1])
    for ext in ("png","pdf","svg"):
        fig.savefig(out/f"version_ladder.{ext}",dpi=180,bbox_inches="tight")
    plt.close(fig)
    lookup = {frozenset((r["version"],r["comparator"])): r for r in compatibility}
    matrix = np.full((len(versions),len(versions)),np.nan)
    text = {}
    for i,a in enumerate(versions):
        for j,b in enumerate(versions):
            if i == j:
                text[i,j] = "identity"
            else:
                row = lookup.get(frozenset((a,b)))
                status = row["decision"] if row else "insufficient_evidence"
                matrix[i,j] = {"insufficient_evidence":0,"within_margin":1,"exceeds_margin":2}[status]
                text[i,j] = {"insufficient_evidence":"insufficient", "within_margin":"within margin", "exceeds_margin":"exceeds margin"}[status]
    fig,ax = plt.subplots(figsize=(max(6,len(versions)*1.5),max(4,len(versions)*.8)))
    ax.imshow(matrix,cmap=ListedColormap(["#e8e8e8","#c7d9e9","#ecc6a4"]),vmin=-.5,vmax=2.5)
    for (i,j),label in text.items():
        ax.text(j,i,label,ha="center",va="center",fontsize=8)
    ticks = [v+(" *" if v==ref else "") for v in versions]
    ax.set_xticks(range(len(versions)),ticks,rotation=35,ha="right")
    ax.set_yticks(range(len(versions)),ticks)
    ax.set_title("Compatibility within prespecified loss margins\n"+title,fontsize=11)
    fig.text(.5,.01,"* Fixed reference. Decisions apply only to tested conditions and evaluable metrics. See compatibility.tsv for exclusions and margins.",ha="center",fontsize=7)
    fig.tight_layout(rect=[0,.04,1,1])
    for ext in ("png","pdf","svg"):
        fig.savefig(out/f"compatibility.{ext}",dpi=180,bbox_inches="tight")
    plt.close(fig)


def _sample_cell(row,metric):
    state=row['comparison_state']
    if state=='no_seeds':return None,'NS'
    if state=='both_not_detected':return None,'ND/ND'
    if state=='partially_evaluable':return None,'PART'
    if state=='not_evaluable':
        codes={'failed':'F','missing':'M','unsupported':'U','invalid':'I','interrupted':'X','blocked':'B','no_seeds':'NS','not_run':'NR'}
        own=[codes[s] for s in row['observed_statuses'] if s in codes]
        reference=[codes[s] for s in row['reference_statuses'] if s in codes]
        return None,('/'.join(own)+(' REF:'+('/'.join(reference)) if reference else '')).strip() or 'NA'
    value=row.get(metric)
    return (None,'NA') if value is None else (value,f'{100*value:.0f}%'+('*' if metric in ('label_set_agreement','class_footprint_jaccard') and not row['full_label_comparable'] else ''))


def _sample_plots(out,rows,ref,synthetic,analysis_scope):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    samples=list(dict.fromkeys(r['sample_id'] for r in rows))
    versions=list(dict.fromkeys(r['version'] for r in rows))
    if not samples or not versions:return []
    lookup={(r['sample_id'],r['version']):r for r in rows}
    metrics=[('detection_agreement','Detection agreement'),('label_set_agreement','Label-set agreement'),('interval_jaccard','Interval Jaccard')]
    mode='Full Suite including CNVkit' if analysis_scope=='full_suite' else 'Fixed-seed AA + AC'
    suffix=' | SYNTHETIC TEST DATA' if synthetic else ''
    legend='Gray: ND/ND both not detected; F failed; M missing; U unsupported; I invalid; NS no seeds; NR not run; PART partial; NA unavailable. REF: fixed reference.'
    caveat='Descriptive similarity, no per-sample confidence interval or compatibility decision. * Shared supported label scope only. Reference self-comparison is an identity anchor.'
    cmap=LinearSegmentedColormap.from_list('observed_similarity',['#D18A50','#FAFAFA','#497CA5']).with_extremes(bad='#dedede')
    n_panels=4 if analysis_scope=='full_suite' else 3
    fig,axes=plt.subplots(n_panels,1,figsize=(max(9,1.1*len(versions)+4),max(7,n_panels*(.27*len(samples)+1.3))),squeeze=False)
    for index,(metric,title) in enumerate(metrics):
        ax=axes[index,0]
        matrix=np.full((len(samples),len(versions)),np.nan)
        labels={}
        for i,sample in enumerate(samples):
            for j,version in enumerate(versions):
                row=lookup[sample,version]
                value,label=_sample_cell(row,metric)
                if row['identity_anchor'] and value is not None:
                    labels[i,j]='REF'
                else:
                    if value is not None:matrix[i,j]=value
                    labels[i,j]=label
        im=ax.imshow(matrix,cmap=cmap,vmin=0,vmax=1,aspect='auto')
        for (i,j),label in labels.items():ax.text(j,i,label,ha='center',va='center',fontsize=8,color='#222222')
        ax.set_xticks(range(len(versions)),[v+(' [ref]' if v==ref else '') for v in versions],rotation=20,ha='right')
        ax.set_yticks(range(len(samples)),samples)
        ax.set_title(title,loc='left',fontsize=10)
        fig.colorbar(im,ax=ax,fraction=.022,pad=.02,label='Observed agreement / similarity',ticks=[0,.5,1])
    if analysis_scope=='full_suite':
        ax=axes[3,0]
        matrix=np.full((len(samples),len(versions)),np.nan)
        for i,sample in enumerate(samples):
            for j,version in enumerate(versions):
                count=lookup[sample,version].get('seed_count')
                if count is not None:matrix[i,j]=count
        countmap=plt.get_cmap('Blues').with_extremes(bad='#dedede')
        countmax=max(1,float(np.nanmax(matrix))) if np.isfinite(matrix).any() else 1
        ax.imshow(matrix,cmap=countmap,vmin=0,vmax=countmax,aspect='auto')
        for i in range(len(samples)):
            for j in range(len(versions)):
                ax.text(j,i,str(int(matrix[i,j])) if np.isfinite(matrix[i,j]) else 'NA',ha='center',va='center',fontsize=8,
                        color='white' if np.isfinite(matrix[i,j]) and matrix[i,j]/countmax>.65 else '#222222')
        ax.set_xticks(range(len(versions)),versions,rotation=20,ha='right');ax.set_yticks(range(len(samples)),samples)
        ax.set_title('Final seed count — recorded seed-stage output, not a positivity or AA/AC equality claim',loc='left',fontsize=10)
    fig.suptitle(f'Per-sample version comparison | Fixed reference: {ref}\n{mode}{suffix}',fontsize=12)
    fig.text(.02,.018,legend+'\n'+caveat,fontsize=7)
    fig.tight_layout(rect=[0,.07,1,.94])
    for ext in ('png','pdf','svg'):fig.savefig(out/f'sample_version_heatmap.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)
    folder=out/'sample_ladders';folder.mkdir(exist_ok=True)
    indices=[]
    for sample in samples:
        fig,axes=plt.subplots(1,n_panels,figsize=(13 if n_panels==4 else 11,max(3.8,.42*len(versions)+2)),sharey=True)
        for ax,(metric,title) in zip(axes,metrics):
            for y,version in enumerate(versions):
                row=lookup[sample,version];value,label=_sample_cell(row,metric)
                if value is None:
                    ax.axhspan(y-.32,y+.32,color='#ededed',zorder=0)
                    ax.text(.5,y,label,ha='center',va='center',fontsize=7,color='#666666')
                else:
                    ax.scatter(value,y,marker='D' if row['identity_anchor'] else 'o',s=28,color='#333333' if row['identity_anchor'] else cmap(value),edgecolor='#555555',linewidth=.4)
                    ax.annotate(label,(value,y),xytext=((5 if value<.15 else -5),5),textcoords='offset points',ha='left' if value<.15 else 'right',fontsize=7)
            ax.set_xlim(-.04,1.06);ax.set_xticks([0,.5,1],['0','50','100']);ax.set_xlabel(title+' (%)')
            ax.grid(axis='x',color='#e7e7e7',linewidth=.6);ax.set_axisbelow(True)
        if analysis_scope=='full_suite':
            ax=axes[-1];counts=[lookup[sample,v].get('seed_count') for v in versions]
            for y,(version,count) in enumerate(zip(versions,counts)):
                if count is None:ax.text(.5,y,'NA',ha='center',va='center',fontsize=7)
                else:
                    ax.scatter(count,y,marker='D' if version==ref else 'o',s=28,color='#333333' if version==ref else '#497CA5')
                    ax.annotate(str(count),(count,y),xytext=(4,5),textcoords='offset points',fontsize=7)
            maximum=max([c for c in counts if c is not None]+[1]);ax.set_xlim(-.2,maximum*1.15+.2)
            ax.set_xlabel('Final seed count');ax.grid(axis='x',color='#e7e7e7',linewidth=.6);ax.set_axisbelow(True)
        axes[0].set_yticks(range(len(versions)),[v+(' [reference]' if v==ref else '') for v in versions]);axes[0].set_ylim(len(versions)-.5,-.5)
        fig.suptitle(f'{sample} | Version ladder | Fixed reference: {ref}\n{mode}{suffix}',fontsize=11)
        fig.text(.02,.015,legend+'\n'+caveat,fontsize=6.5)
        fig.tight_layout(rect=[0,.10,1,.84])
        safe=re.sub(r'[^A-Za-z0-9_.-]+','_',str(sample)).strip('._')[:72] or 'sample'
        safe+='_'+hashlib.sha256(str(sample).encode()).hexdigest()[:8]
        item={'sample_id':sample,'reference_version':ref,'analysis_scope':analysis_scope}
        for ext in ('png','pdf','svg'):
            target=folder/f'{safe}_version_ladder.{ext}';fig.savefig(target,dpi=180,bbox_inches='tight');item[ext]=str(target)
        indices.append(item);plt.close(fig)
    write_tsv(out/'sample_plot_index.tsv',indices)
    return indices


def _class_footprint_plot(out,rows,stats,ref,synthetic,analysis_scope):
    """Separate fourth primary metric: spatial classification, without crowding ladders."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    samples=list(dict.fromkeys(r['sample_id'] for r in rows));versions=list(dict.fromkeys(r['version'] for r in rows))
    if not samples or not versions:return
    lookup={(r['sample_id'],r['version']):r for r in rows}
    fig,(ax,heat)=plt.subplots(1,2,figsize=(max(11,.8*len(versions)+6),max(4.5,.4*max(len(versions),len(samples))+2.6)),gridspec_kw={'width_ratios':[1,2]})
    for y,version in enumerate(versions):
        if version==ref:
            if any(lookup[s,ref]['class_footprint_jaccard'] is not None for s in samples):ax.scatter(1,y,marker='D',s=28,color='#333333')
            else:ax.text(.5,y,'NA',ha='center',va='center',color='#777777')
            continue
        stat=next((r for r in stats if r['version']==version and r['comparator']==ref and r['subset']=='signal_positive' and r['metric']=='class_footprint_jaccard'),None)
        if not stat or stat['estimate'] is None:
            ax.text(.5,y,'not evaluable',ha='center',va='center',fontsize=8,color='#777777');continue
        value,lo,hi=stat['estimate'],stat['conservative_ci_low'],stat['conservative_ci_high']
        ax.errorbar(value,y,xerr=[[value-lo],[hi-value]],fmt='o',capsize=3,color='#497CA5',markersize=5)
        ax.text(.02,y-.2,f"n={stat['n_samples']}",fontsize=7,color='#666666')
    ax.set_yticks(range(len(versions)),[v+(' [ref]' if v==ref else '') for v in versions]);ax.set_ylim(len(versions)-.5,-.5)
    ax.set_xlim(-.03,1.04);ax.set_xticks([0,.5,1],['0','50','100']);ax.set_xlabel('Class-footprint Jaccard (%)')
    ax.set_title('Across evaluable signal-positive samples',fontsize=10);ax.grid(axis='x',color='#eeeeee');ax.set_axisbelow(True)
    cmap=LinearSegmentedColormap.from_list('class_footprint',['#D18A50','#FAFAFA','#497CA5']).with_extremes(bad='#dedede')
    values=np.full((len(samples),len(versions)),np.nan);labels={}
    for i,sample in enumerate(samples):
        for j,version in enumerate(versions):
            row=lookup[sample,version];value,label=_sample_cell(row,'class_footprint_jaccard')
            if row['identity_anchor'] and value is not None:labels[i,j]='REF'
            else:
                labels[i,j]=label
                if value is not None:values[i,j]=value
    heat.imshow(values,cmap=cmap,vmin=0,vmax=1,aspect='auto')
    for (i,j),label in labels.items():heat.text(j,i,label,ha='center',va='center',fontsize=7)
    heat.set_xticks(range(len(versions)),versions,rotation=30,ha='right');heat.set_yticks(range(len(samples)),samples)
    heat.set_title('Per sample: observed spatial classification similarity',fontsize=10)
    mode='Full Suite including CNVkit' if analysis_scope=='full_suite' else 'Fixed-seed AA + AC'
    fig.suptitle(f'Category-specific genomic footprints | Fixed reference: {ref}\n{mode}'+(' | SYNTHETIC TEST DATA' if synthetic else ''),fontsize=11)
    confidence=stats[0]['confidence_level'] if stats else .95
    fig.text(.02,.015,f'Macro-mean Jaccard across labels present in either version, within shared supported classes. Left: sample means and conservative {confidence*100:g}% bounds.\nGray/NA: unavailable class locations; ND/ND: both negative; PART: incomplete. No per-sample confidence intervals. This does not establish cycle order or breakpoint topology.',fontsize=7)
    fig.tight_layout(rect=[0,.10,1,1])
    for ext in ('png','pdf','svg'):fig.savefig(out/f'class_footprint_comparison.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)


def _official_report_links(out):
    root=out/'official_comparisons'
    tables=[(name,title) for name,title in (
        ('official_metrics.tsv','官方原始指标'),('index.tsv','比较状态与原始结果索引'),
        ('official_plots.tsv','全部官方距离图索引')) if (root/name).is_file()]
    if not tables:return []
    lines=['','## 官方结构比较与距离天梯图','',
           '官方距离用于描述已观察到的结构差异，不是兼容接受阈值。无重建结构、失败、缺失或不可评价项不计为相似度，也不填成零距离。固定参考及指标单位见各图和原始指标表。','',
           '；'.join(f'[{title}](official_comparisons/{name})' for name,title in tables)+'。']
    index=root/'official_plots.tsv'
    if not index.is_file():return lines
    with index.open(newline='',encoding='utf-8') as handle:
        figures=list(csv.DictReader(handle,delimiter='\t'))
    linked=set()
    for row in figures:
        path=Path(row.get('path',''))
        if not path.is_absolute():path=root/path
        try:relative=path.resolve().relative_to(out.resolve())
        except ValueError:continue
        if path.suffix.lower()!='.png' or not path.is_file() or relative in linked:continue
        if not linked:lines += ['','以下为可用图片的快捷入口；PDF、SVG及更多图片见完整图索引。','']
        title=' / '.join(str(row.get(key,'')) for key in ('sample_id','seed_id','figure'))
        title=title.replace('[','\\[').replace(']','\\]').replace('\n',' ')
        lines.append(f'- [{title}](<{relative.as_posix()}>)')
        linked.add(relative)
        if len(linked)>=10:break
    return lines


def write_report(outdir, records, pair_rows, stats_rows, compat_rows, ref, synthetic=False, analysis_scope=None):
    analysis_scope=analysis_scope or next((r.get('analysis_scope') for r in stats_rows+compat_rows if r.get('analysis_scope')),'fixed_seed')
    full=analysis_scope=='full_suite'
    from .statistics import build_sample_similarity
    sample_rows=build_sample_similarity(records,ref,pair_rows,analysis_scope)
    out = Path(outdir)
    out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/"normalized_results.tsv",records)
    write_tsv(out/"differences.tsv",pair_rows)
    write_tsv(out/"statistics.tsv",stats_rows)
    write_tsv(out/"compatibility.tsv",compat_rows)
    write_tsv(out/'sample_similarity.tsv',sample_rows)
    footprint_rows=[]
    for row in pair_rows:
        for detail in row.get('class_footprint_details',[]):
            footprint_rows.append({**{k:row.get(k) for k in ('sample_id','seed_id','version','comparator','reference_version','analysis_scope')},**detail})
    write_tsv(out/'class_footprint_differences.tsv',footprint_rows)
    anomalies = [r for r in records if r["status"] not in {"ok","not_detected","no_seeds"} or r.get("warnings")]
    write_tsv(out/"anomalies.tsv",anomalies)
    write_tsv(out/'no_seed_samples.tsv',[r for r in records if r['status']=='no_seeds'])
    if full:
        seed_fields=('reference_version','analysis_scope','version','comparator','sample_id','status','comparator_status',
                     'seed_count','comparator_seed_count','seed_count_difference','seed_interval_jaccard','seed_interval_distance',
                     'seed_intervals_added','seed_intervals_removed','detection_stage','comparator_detection_stage',
                     'aa_executed','comparator_aa_executed','aa_reconstruction_performed','comparator_aa_reconstruction_performed',
                     'ac_executed','comparator_ac_executed')
        write_tsv(out/'seed_stage_differences.tsv',[{k:r.get(k) for k in seed_fields} for r in pair_rows])
    (out/"normalized_results.json").write_text(json.dumps(records,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    status_counts = Counter((r["version"],r["status"]) for r in records)
    write_tsv(out/"status_counts.tsv",[dict(version=v,status=s,count=n) for (v,s),n in sorted(status_counts.items())])
    plot_error = ""
    try:
        _plots(out,records,stats_rows,compat_rows,ref,synthetic,analysis_scope)
        _sample_plots(out,sample_rows,ref,synthetic,analysis_scope)
        _class_footprint_plot(out,sample_rows,stats_rows,ref,synthetic,analysis_scope)
    except ImportError as exc:
        plot_error = f"绘图依赖缺失：{exc}；表格已保存。安装 requirements.txt 后重新 report。"
        (out/"plot_error.txt").write_text(plot_error+"\n",encoding="utf-8")
    scope_text=("范围：完整 Suite 流程，包括重新计算 CNVkit、自动种子处理、AA 重建与 AC 分类。流程结果差异可能来自其中任何步骤；不能直接归因为 AA 或 AC 版本。每个样本/版本只有一条整体结局；分类与区段比较使用样本内全部 feature 的并集，不按相同 amplicon 编号对齐，也不代表逐扩增子拓扑相同。" if full else
                "范围：固定种子 AA＋AC。种子是待检查区域，不是已确认扩增；允许 AA 必要的背景估计与相连区域扩展。未验证 CNVkit 版本影响，不做全基因组扩增遍历或全基因组基因注释。")
    lines = ["# AmpliconSuite "+("完整流程版本比较" if full else "固定种子版本比较"), "", f"**固定参考版本：`{ref}`。**",
             "", "**这是人工合成测试数据，仅用于程序验证，不是真实样本结论。**" if synthetic else "本报告仅描述本次输入与固定条件下的配对比较。",
             "", scope_text,
             "", "每个样本、种子和版本的状态及多标签在 `normalized_results.tsv`；逐对增检、漏检、类别改变和新增/丢失区段在 `differences.tsv`。异常步骤、原因与证据在 `anomalies.tsv` 及原日志；没有证据的原因仍待查。",
             "", "## 兼容判断", "", ("| 候选版本 | 对照版本 | 结论 | 有信号样本数 | 不可评价样本数 |" if full else "| 候选版本 | 对照版本 | 结论 | 有信号样本数 | 不可评价区域数 | 无输入种子样本数 |"), ("|---|---|---|---:|---:|" if full else "|---|---|---|---:|---:|---:|")]
    if any(r.get('reference_scope') == 'human_only_incomplete_reference' for r in records):
        lines[4:4] = ["", "**本次参考范围仅为人类参考；BAM 中的病毒 contig 未获完整参考验证。** 所有结果与兼容判断仅适用于明确保留的人类参考范围，不能外推到病毒相关结构或完整原始参考。"]
    for r in compat_rows:
        lines.append(f"| {r['version']} | {r['comparator']} | {r['conclusion']} | {r['signal_samples']} | {r['excluded_regions']} |"+("" if full else f" {r.get('no_seed_samples',0)} |"))
    lines += ["", "结论依赖预先配置的允许差异范围和独立样本数。没有设定阈值、样本不足、类别支持不完整或存在不可评价区域时，不能作整体兼容承诺。`P>0.05` 不等于兼容。",
              "", ("完整流程成功结束但最终自动种子为空时，`not_detected` 是完整流程的未检出结局，`detection_stage=seed_processing` 表示未检出由种子步骤决定。包装器仍可能调用 AA 生成零扩增子摘要；`aa_executed`/`ac_executed` 记录是否调用，`aa_reconstruction_performed` 另记有无实际扩增重建，只有执行字段明确为 false 才称未调用。这不证明 AA/AC 相同或生物学上不存在扩增。失败、缺文件仍单列为不可评价，不能当作未检出。`seed_stage_differences.tsv` 单列种子数量差和种子区段 Jaccard，双方种子都空时 Jaccard 留空，不记100%结构一致。" if full else "`no_seeds` 表示该样本没有输入候选种子，记录在 `no_seed_samples.tsv`；它不是阴性、运行失败或缺失的真实区域，不进入种子数量和配对区域分母。存在这类样本时，不能把剩余可评价样本的结果概括为整批样本已验证兼容。"),
              "", "## 统计口径", "", ("每个完整流程样本为一个独立配对单位，先合并该样本全部明确分类 feature 的标签与区段，再对样本等权比较；没有把自动种子数当作样本数。" if full else "同一样本内先对种子求均值，再对样本等权平均。")+"置信区间与置换均以整个样本为单位，不能把多个种子当作独立样本。若多个样本来自同一患者，输入前需合并统计单位或扩展为患者级聚类；目前默认样本相互独立。",
              "", "`statistics.tsv` 同时报告效应大小、样本整簇 bootstrap 区间、Hoeffding 保守区间和配对整样本交换检验；P 值对全部差异检验做 Holm 校正。交换检验需要零假设下版本标签可交换。零观察差异时 bootstrap 区间可能退化，因此兼容判定使用非退化的 Hoeffding 上界。",
              "", "一致性判定的四个主要指标是有信号"+("样本" if full else "区域")+"中的检出失配率、标签集合失配率、1−总体区段 Jaccard、分类特异区段距离（class_footprint_distance）。最后一项按共同支持类别分别计算区段Jaccard，再对任一版本出现的类别等权平均；可发现标签集合和总体区段完全相同、但各类别发生位置互换的差异。feature信息缺失时记未知，不能填零；四项允许范围均需预设。检验 H0：平均失配 ≥ 允许上限。仅当保守单侧上界低于预设上限才支持初步兼容；跨全部版本对和四指标作 Bonferroni 控制。此法较保守，少量样本通常只能得到证据不足。",
              "", "兼容表的 compatibility_metric_scope 明确限定为检出、分类集合、总体区段和分类特异区段指标；即使这些指标通过，也不证明相同路径顺序、断点拓扑、生物学真实性或软件版本全面可互换。官方EC距离仍逐项展示，除非另行预设有效接受范围，否则不并入等效性判断。",
              "", "全体可评价"+("样本与任一版本有信号的样本" if full else "区域与任一版本有信号的区域")+"分别报告，天梯图使用后者，防止全阴性结果抬高一致率。所有一致程度都不是生物学真实性概率。比较只适用于已验证样本条件、区域、参数、参考资源和指标。",
              "", "![版本天梯图](version_ladder.png)", "", "![兼容表](compatibility.png)", "", "更完整的方法、公式和依据见程序目录 `docs/statistical_methods.md`。"]
    lines += ['', '## 每个输入样本的版本差异', '',
              '`sample_similarity.tsv` 给出每个样本相对固定参考的检出一致、标签集合一致、区段Jaccard与状态。每样本单独天梯图的路径见 `sample_plot_index.tsv`；这些图不画单样本置信区间，也不给单样本统计等效结论。灰格逐项区分失败、缺失、不支持、无效、无种子、部分可评价和双方未检出。', '',
              '相似度较高只描述已观察到的结果。官方结构匹配阈值（若另行输出）用于匹配候选结构，不能直接作为本项目的等效性或整体兼容阈值。未预设兼容范围时仍为证据不足，不能用 P＞0.05 代替。', '',
              '![逐样本版本总览](sample_version_heatmap.png)']
    lines += ['', '分类所在区段的逐类别新增/丢失明细见 `class_footprint_differences.tsv`；第四项主指标单独绘图，避免与只比较标签是否出现的指标混淆。', '',
              '![分类特异区段一致性](class_footprint_comparison.png)']
    lines += _official_report_links(out)
    if plot_error:
        lines += ["",plot_error]
    (out/"report.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    return [str(p) for p in sorted(out.iterdir()) if p.is_file()]


def write_reports(output_dir, records, reference_version, options=None, synthetic=False):
    from .statistics import analyze
    pairs,stats,compat = analyze(records,reference_version,options)
    return write_report(output_dir,records,pairs,stats,compat,reference_version,synthetic,(options or {}).get('analysis_scope','fixed_seed'))
