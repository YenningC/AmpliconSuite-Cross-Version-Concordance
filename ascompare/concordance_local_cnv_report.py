"""Describe local CNVkit contrasts using the input sample's actual CNVkit version."""
from collections import Counter, defaultdict
from pathlib import Path

from .concordance_report import _bounds, _upstream_context, _upstream_evaluable
from .concordance_upstream_similarity import _clipped_footprint
from .reporting import write_tsv
from .statistics import _independent_means, merge_intervals


def _measurement(row, sample, version):
    if not row:
        return None, 'missing_record'
    if not _upstream_evaluable(row):
        return None, row.get('reason') or row.get('status') or 'unavailable'
    role = 'source_local_standard' if version == sample.get('source_cnvkit_version') else 'comparison'
    if (row.get('actual_cnvkit_version') != version or not version
            or row.get('source_version') != sample['source_version']
            or row.get('source_filter_suite_version') != sample['source_version']
            or row.get('source_cnvkit_version') != sample.get('source_cnvkit_version')
            or row.get('measurement_role') != role
            or not row.get('source_filter_aa_version') or not row.get('local_design_id')):
        return None, 'measurement_identity_unverified'
    context = _upstream_context(row)
    if context is None:
        return None, 'missing_candidate_context'
    pieces = _clipped_footprint(row.get('retained_intervals'), context[:3])
    if pieces is None:
        return None, 'missing_or_invalid_final_seed_coordinates'
    pieces = merge_intervals(pieces)
    fraction = sum(b-a for _, a, b in pieces)/(context[2] - context[1])
    retained = row['status'] == 'retained'
    if retained != (fraction >= context[-1]):
        return None, 'retention_status_inconsistent_with_coordinates'
    return dict(context=(row['local_design_id'], row['source_filter_aa_version'], *context),
                intervals=pieces, retained=retained, fraction=fraction), None


def _overlap(target, baseline):
    def length(parts):
        return sum(b-a for _, a, b in merge_intervals(parts))
    union = length(target + baseline)
    intersection = length(target) + length(baseline) - union
    return intersection, union, intersection/union if union else None


def _calibration(row, measurement):
    if measurement and row.get('calibration_role') == 'source_control':
        status = row.get('calibration_status')
        if status == 'failed' or status == 'passed' and measurement['retained']:
            return status
    return 'unverified'


def _seed_pair(meta, seed, sample, trow, brow):
    a, a_reason = _measurement(trow, sample, meta['target_cnvkit_version'])
    b, b_reason = _measurement(brow, sample, meta['source_cnvkit_version'])
    matching = bool(a and b and a['context'] == b['context'])
    if matching:
        intersection, union, jaccard = _overlap(a['intervals'], b['intervals'])
        status = 'ok' if union else 'both_empty'
    else:
        intersection = union = jaccard = None
        status = 'source_baseline_unavailable' if b is None else 'target_unavailable' if a is None else 'pair_context_mismatch'
    detail = dict(meta, seed_id=seed, status=status, paired=matching,
                  target_status=trow.get('status', 'missing'), baseline_status=brow.get('status', 'missing'),
                  target_reason=a_reason, baseline_reason=b_reason,
                  target_retained=a['retained'] if a else None, baseline_retained=b['retained'] if b else None,
                  target_candidate_coverage=a['fraction'] if a else None, baseline_candidate_coverage=b['fraction'] if b else None,
                  retention_delta=int(a['retained']) - int(b['retained']) if matching else None,
                  intersection_bp=intersection, union_bp=union, interval_jaccard=jaccard,
                  baseline_original_calibration_status=_calibration(brow, b),
                  baseline_original_calibration_reason=brow.get('calibration_reason'))
    for label, raw, measured in (('target', trow, a), ('baseline', brow, b)):
        for field in ('local_design_id', 'measurement_origin', 'actual_cnvkit_version', 'source_filter_aa_version', 'provenance'):
            detail[label+'_'+field] = raw.get(field)
        detail[label+'_seed'] = {k: raw.get('seed', {}).get(k) for k in ('chrom','start','end')} if raw.get('seed') else None
        detail[label+'_intervals'] = measured['intervals'] if measured else None
    return detail


def _sample_summary(meta, details, unavailable, declared):
    # A matched seed pair may have no bp in either final seed set. That pair is
    # informative for retention counts and supplies no interval similarity.
    paired = [r for r in details if r['paired']]
    n, expected = len(paired), len(details)
    complete = bool(declared and n and n == expected and not unavailable)
    summary = dict(meta, candidate_universe_declared=declared, expected_seed_count=expected if declared else None,
                   paired_seed_count=n, unavailable_seed_count=expected-n if declared else None,
                   pair_status_counts=dict(Counter(r['status'] for r in details)))
    for role in ('target', 'baseline'):
        valid = [r[role+'_retained'] for r in details if r[role+'_retained'] is not None]
        count = sum(valid)
        summary[role+'_status_counts'] = dict(Counter(r[role+'_status'] for r in details))
        versions = sorted({r[role+'_source_filter_aa_version'] for r in details if r[role+'_source_filter_aa_version']})
        designs = sorted({r[role+'_local_design_id'] for r in details if r[role+'_retained'] is not None})
        if len(versions) > 1 or len(designs) > 1:
            raise ValueError('Inconsistent source filter or local design within a sample and CNVkit version')
        summary[role+'_source_filter_aa_version'] = versions[0] if versions else None
        summary[role+'_local_design_id'] = designs[0] if designs else None
        summary[role+'_evaluable_count'] = len(valid)
        summary[role+'_retained_count'] = count if valid else None
        summary[role+'_panel_retention'] = count/expected if declared and expected and len(valid) == expected else None
        summary['paired_'+role+'_retained_count'] = sum(r[role+'_retained'] for r in paired) if n else None
    summary['baseline_original_calibration_status'] = 'unverified'
    controls = {r['baseline_original_calibration_status'] for r in details}
    if 'failed' in controls:
        summary['baseline_original_calibration_status'] = 'failed'
    elif declared and expected and controls == {'passed'}:
        summary['baseline_original_calibration_status'] = 'passed'

    extra = sum(r['retention_delta'] == 1 for r in paired)
    missed = sum(r['retention_delta'] == -1 for r in paired)
    delta = (extra - missed)/n if n else None
    left = [i for row in paired for i in row['target_intervals']]
    right = [i for row in paired for i in row['baseline_intervals']]
    intersection, union, jaccard = _overlap(left, right)
    summary.update(extra_count=extra if n else None, missed_count=missed if n else None,
                   observed_retention_delta=delta, retention_delta=delta if complete else None,
                   intersection_bp=intersection if n else None, union_bp=union if n else None,
                   observed_interval_jaccard=jaccard, interval_jaccard=jaccard if complete else None)
    if unavailable:
        status = unavailable
    elif not declared:
        status = 'candidate_universe_undeclared'
    elif not expected:
        status = 'no_candidates'
    elif complete:
        status = 'ok' if union else 'both_empty'
    elif n:
        status = 'partial'
    elif not summary['baseline_evaluable_count']:
        status = 'source_baseline_unavailable'
    elif not summary['target_evaluable_count']:
        status = 'target_unavailable'
    elif any(r['status'] == 'pair_context_mismatch' for r in details):
        status = 'pair_context_mismatch'
    else:
        status = 'no_shared_evaluable_candidates'
    summary['status'] = status
    return summary


def _cohort(samples, patient_ids, alpha):
    grouped = defaultdict(list)
    strata = ('reference_version', 'source_version', 'source_cnvkit_version', 'target_cnvkit_version',
              'reference_build', 'reference_mode')
    for sample in samples:
        grouped[tuple(sample[k] for k in strata)].append(sample)
    rows = []
    for key, group in grouped.items():
        result = dict(zip(strata, key))
        identity = group[0]['identity_anchor']
        selected = [s['sample_id'] for s in group]
        result.update(identity_anchor=identity, target_suite_versions=sorted({v for s in group for v in s['target_suite_versions']}),
                      status_counts=dict(Counter(r['status'] for r in group)),
                      n_selected_samples=len(group), n_selected_independent_units=len({patient_ids[s] if patient_ids else s for s in selected}),
                      independence_unit='patient' if patient_ids is not None else 'sample',
                      confidence_level=1-alpha, interval_method='descriptive_pointwise_two_sided_Hoeffding',
                      interpretation='descriptive_matched_local_CNVkit_only', synthetic=group[0]['synthetic'])
        for k in ('paired_seed_count', 'paired_target_retained_count', 'paired_baseline_retained_count',
                  'target_evaluable_count', 'target_retained_count', 'baseline_evaluable_count', 'baseline_retained_count',
                  'extra_count', 'missed_count'):
            values = [s[k] for s in group if s[k] is not None]
            result[k] = sum(values) if values else None
        result['expected_seed_count'] = sum(s['expected_seed_count'] for s in group) if all(
            s['expected_seed_count'] is not None for s in group) else None
        for metric in ('interval_jaccard', 'retention_delta', 'target_panel_retention', 'baseline_panel_retention'):
            values = {s['sample_id']: s[metric] for s in group if s[metric] is not None}
            independent = _independent_means(values, patient_ids).tolist()
            limits = (-1., 1.) if metric == 'retention_delta' else (0., 1.)
            bounds = _bounds(independent, alpha, limits=limits)
            if identity and metric in ('interval_jaccard', 'retention_delta'):
                bounds = None, None
            result.update({metric: sum(independent)/len(independent) if independent else None,
                           metric+'_n_samples': len(values), metric+'_n_independent_units': len(independent),
                           metric+'_ci_low': bounds[0], metric+'_ci_high': bounds[1]})
        rows.append(result)
    return rows


def write_local_cnv_comparator(out, records, config, patient_ids=None, alpha=.05, synthetic=False):
    """Write the local-only, source-CNVkit-standard contrast; no runtime/file dereferencing."""
    from .local_cnv_comparator import comparison_panel

    if not 0 < alpha < 1:
        raise ValueError('alpha must be in (0, 1)')
    out = Path(out)
    by_sample = defaultdict(dict)
    known = {s['sample_id'] for s in config['samples']}
    for row in records or []:
        sid = row['sample_id']
        key = row['cnvkit_version'], row['seed_id']
        if sid not in known or key in by_sample[sid]:
            raise ValueError('Unknown sample or duplicate source-matched local CNVkit measurement')
        by_sample[sid][key] = row
    seeds, summaries = [], []
    for sample in config['samples']:
        sid = sample['sample_id']
        source = sample.get('source_cnvkit_version')
        panel = comparison_panel(sample, config)
        local = by_sample[sid]
        if any(version not in panel for version, seed in local):
            raise ValueError('Local CNVkit measurement version outside the configured comparison panel')
        declared = 'upstream_seed_ids' in sample
        ids = sample.get('upstream_seed_ids', [])
        if declared and any(seed not in ids and seed != '__no_seeds__' for _, seed in local):
            raise ValueError('Local comparator record outside declared candidate universe')
        mode = config['mode']
        unavailable = None
        if mode == 'results_only':
            unavailable = 'not_applicable_results_only'
        elif records is None:
            unavailable = 'not_recorded_in_run'
        elif not source:
            unavailable = 'source_cnvkit_version_unknown'
        for target, suites in panel.items():
            meta = dict(sample_id=sid, source_version=sample['source_version'], reference_version=config['reference_version'],
                        reference_build=sample.get('reference_build', config.get('reference',{}).get('build')),
                        reference_mode=config.get('local_cnv',{}).get('reference_mode','matched_normal'),
                        source_cnvkit_version=source, target_cnvkit_version=target, target_suite_versions=suites,
                        identity_anchor=bool(source and target == source), synthetic=synthetic)
            details = []
            for seed in ids:
                baseline = local.get((source,seed), {}) if not unavailable else {}
                observation = local.get((target,seed), {}) if not unavailable else {}
                row = _seed_pair(meta, seed, sample, observation, baseline)
                if unavailable:
                    row['status'] = unavailable
                details.append(row)
            seeds.extend(details)
            summaries.append(_sample_summary(meta, details, unavailable, declared))
    cohort = _cohort(summaries, patient_ids, alpha)
    files = {name: 'local_cnv_comparator_'+name+'.tsv' for name in ('seeds','samples','cohort')}
    for key, values in (('seeds', seeds), ('samples', summaries), ('cohort', cohort)):
        write_tsv(out / files[key], values)
    return dict(samples=summaries, cohort=cohort, files=files)
