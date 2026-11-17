"""Per-sample native EC distance ladders: observed values, not probabilities."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path

from .common import table
from .full_suite import WHOLE_SAMPLE


NORMALIZED = (
    ('cn_hamming_norm_dist','Genomic footprint'),
    ('cn_cos_dist','Copy-number cosine'),
    ('cn_jc_dist','Copy-number min-max'),
    ('fragments_dist','Fragment multiplicity'),
    ('cycles_dist','Path multiplicity'),
    ('breakpoint_dist','Breakpoint matching'),
)


def plot_groups(rows):
    """No pooling of selected seeds: each sample/seed/reference gets a panel set."""
    grouped=defaultdict(list)
    for row in rows:
        if row.get('tool')=='ec':
            grouped[(row['sample_id'],row['seed_id'],row['reference_version'])].append(row)
    output=[]
    for (sample,seed,ref), group in sorted(grouped.items()):
        versions=list(dict.fromkeys([ref]+[r['candidate_version'] for r in group]))
        points={};statuses={}
        for row in group:
            version=row['candidate_version'];metric=row.get('metric','')
            status=row.get('status','missing')
            statuses.setdefault(version,status)
            if not metric:continue
            key=(version,metric)
            if key in points:raise ValueError('Duplicate official scalar metric for sample/seed/version')
            value=None
            try:
                candidate=float(row['value']);maximum=5 if metric=='total_cost' else 1
                if status=='complete' and math.isfinite(candidate) and 0<=candidate<=maximum:
                    value=candidate
                elif status=='complete':status='out_of_range'
            except (ValueError,TypeError,KeyError):
                if status=='complete':status='missing_value'
            points[key]={'value':value,'status':status}
        output.append(dict(sample_id=sample,seed_id=seed,reference_version=ref,versions=versions,
                           points=points,statuses=statuses))
    return output


def _file_id(text):
    clean=re.sub(r'[^A-Za-z0-9_.-]','_',text)[:100]
    return clean+'_'+hashlib.sha256(text.encode()).hexdigest()[:8]


def plot_official_metrics(out,records=None):
    """Read out/official_comparisons/official_metrics.tsv and save PNG/PDF/SVG.

    Whole-sample inputs are preferred by the caller; any seed-specific inputs
    are clearly labeled and plotted separately, never averaged into a sample.
    Each point is a single observed comparison: no per-sample confidence bar.
    """
    out=Path(out);root=out/'official_comparisons';source=root/'official_metrics.tsv'
    if not source.is_file():return []
    with source.open(newline='') as stream:groups=plot_groups(list(csv.DictReader(stream,delimiter='\t')))
    if not groups:return []
    os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'ascompare-matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,
                         'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    figure_dir=root/'figures';figure_dir.mkdir(exist_ok=True)
    primary={}
    limited_samples=set();synthetic_config=False
    if (out/'config.resolved.json').is_file():
        synthetic_config=bool(json.loads((out/'config.resolved.json').read_text()).get('synthetic'))
    if records is None and (out/'records.json').is_file():
        records=json.loads((out/'records.json').read_text())
    if records is not None:
        primary={(r['sample_id'],r['seed_id'],r['version']):r for r in records}
        limited_samples={r['sample_id'] for r in records if r.get('reference_scope')=='human_only_incomplete_reference'}
    generated=[]
    for group in groups:
        sample,seed,ref=(group[k] for k in ('sample_id','seed_id','reference_version'))
        versions=group['versions'];colors=['#4178A9','#D28A45','#5D8F83','#7C6AA6']
        synthetic=synthetic_config or any('SYNTHETIC' in str(x).upper() for x in (sample,seed,ref))
        title=f'Sample: {sample} | Fixed reference: {ref}'
        title+='\nWhole-sample reconstructed path/cycle collection' if seed in (WHOLE_SAMPLE,'whole_sample') else f'\nSelected seed: {seed} | reconstructed path/cycle collection'
        if synthetic:title+=' | SYNTHETIC TEST DATA'
        if sample in limited_samples:title+='\nScope: human reference only; viral contigs unavailable'
        ref_record=primary.get((sample,seed,ref))
        if ref_record is not None:
            status=ref_record.get('status','missing')
            run_dir=ref_record.get('run_dir')
            if status not in ('ok','not_detected'):
                reference_unavailable='reference '+status
            elif not run_dir or not any((Path(run_dir)/'aa').glob('*_amplicon*_cycles.txt')):
                reference_unavailable='no reconstructed structure' if status=='not_detected' else 'native structure unavailable'
            else:
                reference_unavailable=None
        else:
            # A completed native comparison establishes an available reference
            # structure. Absent that evidence, do not invent a zero-distance dot.
            reference_unavailable=None if any(p['value'] is not None for p in group['points'].values()) else 'reference structure unavailable'
        def draw(ax,metric,maximum,label):
            for y,version in enumerate(versions):
                if version==ref:
                    if reference_unavailable:
                        ax.axhspan(y-.32,y+.32,color='#f1f1f1',zorder=0)
                        ax.text(maximum*.5,y,reference_unavailable,ha='center',va='center',color='#777777',fontsize=8)
                    else:
                        ax.scatter(0,y,marker='D',s=28,color='#333333',zorder=3)
                    continue
                missing_status=group['statuses'].get(version,'not_evaluable')
                if missing_status=='complete':missing_status='metric_missing'
                point=group['points'].get((version,metric),{'value':None,'status':missing_status})
                if point['value'] is None:
                    ax.axhspan(y-.32,y+.32,color='#f1f1f1',zorder=0)
                    ax.text(maximum*.5,y,point['status'].replace('_',' '),ha='center',va='center',color='#777777',fontsize=8)
                else:
                    value=point['value'];ax.scatter(value,y,s=32,color=colors[(y-1)%len(colors)],zorder=3)
                    align='right' if value>maximum*.8 else 'left';offset=-maximum*.025 if align=='right' else maximum*.025
                    ax.text(value+offset,y,f'{value:.3f}',ha=align,va='center',fontsize=7)
            ax.set_xlim(-.035*maximum,1.035*maximum)
            ax.set_ylim(len(versions)-.5,-.5)
            ax.set_yticks(range(len(versions)),[v+(' [reference]' if v==ref else '') for v in versions])
            ax.set_xticks([i/4 for i in range(5)] if maximum==1 else range(6))
            ax.set_title(label,fontsize=10)
            ax.set_xlabel('Distance (lower = closer)' if maximum==1 else 'Weighted total cost (lower = closer)')
            ax.grid(axis='x',color='#e6e6e6',linewidth=.6);ax.set_axisbelow(True)
        fig,axes=plt.subplots(2,3,figsize=(14,max(6.6,len(versions)*.8+3.1)))
        for index,(ax,(metric,label)) in enumerate(zip(axes.flat,NORMALIZED)):
            draw(ax,metric,1,label)
            if index%3:ax.tick_params(axis='y',labelleft=False)
        fig.suptitle('Official EC-comparator distance ladder\n'+title,fontsize=11)
        foot='Native distances; no confidence interval for one sample. Gray = unavailable, not zero. Reference diamond = identity anchor.\nIndependent path/cycle collection, not one physical molecule. No official cross-version acceptance cutoff identified in inspected sources.'
        fig.text(.5,.015,foot,ha='center',fontsize=7)
        fig.tight_layout(rect=[0,.075,1,.91])
        prefix=_file_id(sample+'__'+seed+'__'+ref)
        for extension in ('png','pdf','svg'):
            path=figure_dir/(prefix+'_normalized_distances.'+extension);fig.savefig(path,dpi=180,bbox_inches='tight')
            generated.append(dict(sample_id=sample,seed_id=seed,reference_version=ref,figure='normalized_distances',path=str(path.resolve())))
        plt.close(fig)
        fig,ax=plt.subplots(figsize=(8,max(3.8,len(versions)*.55+2.6)))
        draw(ax,'total_cost',5,'Official weighted total cost: fixed all-enabled configuration')
        fig.suptitle(title,fontsize=10)
        fig.text(.5,.015,'Weights sum to 5. Lower cost means closer reconstructed collections; it is not a compatibility probability.\nSingle observed comparison; no confidence interval. Gray = unavailable. Reference diamond = identity anchor.',ha='center',fontsize=7)
        fig.tight_layout(rect=[0,.09,1,.88])
        for extension in ('png','pdf','svg'):
            path=figure_dir/(prefix+'_total_cost.'+extension);fig.savefig(path,dpi=180,bbox_inches='tight')
            generated.append(dict(sample_id=sample,seed_id=seed,reference_version=ref,figure='total_cost',path=str(path.resolve())))
        plt.close(fig)
    table(root/'official_plots.tsv',generated)
    return generated
