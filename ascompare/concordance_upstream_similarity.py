"""Descriptive comparison of native final-seed footprints inside the original candidates."""

from collections import Counter, defaultdict

from .statistics import merge_intervals


FOOTPRINT_SCOPE = 'native_final_seeds_intersected_with_original_candidate_panel'


def _clipped_footprint(intervals, candidate):
    """Return None when BED evidence is unavailable or malformed; [] is a valid negative."""
    if not isinstance(intervals, (list, tuple)):
        return None
    chrom, start, end = candidate
    pieces = []
    for interval in intervals:
        if not isinstance(interval, (list, tuple)) or len(interval) != 3:
            return None
        c, a, b = interval
        if (not isinstance(c, str) or not c or type(a) is not int or type(b) is not int
                or a < 0 or b <= a):
            return None
        if c == chrom and a < end and start < b:
            pieces.append((c, max(a, start), min(b, end)))
    return pieces


def _evidence(row):
    from .concordance_report import _upstream_context

    if not row or not row['evaluable']:
        return 'unavailable_seed_result', None, None
    context = _upstream_context(row)
    if context is None:
        return 'missing_candidate_context', None, None
    pieces = _clipped_footprint(row.get('retained_intervals'), context[:3])
    if pieces is None:
        return 'missing_or_invalid_final_seed_coordinates', context, None
    return 'ok', context, pieces


def _count_fields(sample, records):
    valid = [row for row in records if row['evaluable']]
    retained = sum(row['status'] == 'retained' for row in valid)
    n, expected = len(valid), sample['expected_seeds']
    if not sample['candidate_universe_declared']:
        status = 'candidate_universe_undeclared'
    elif not expected:
        status = 'no_candidates'
    elif not n:
        status = 'unavailable'
    elif n < expected:
        status = 'partial'
    else:
        status = 'complete'
    return dict(upstream_retained_count=retained, upstream_evaluable_count=n,
                upstream_retained_ratio=retained / n if n else None,
                upstream_panel_retained_ratio=retained / expected if status == 'complete' else None,
                upstream_unevaluable_count=expected - n, upstream_count_status=status)


def _pair_footprints(target, comparator, target_evidence, comparator_evidence):
    target_intervals, comparator_intervals = [], []
    missing = Counter()
    paired = 0
    for seed_id, (state, context, pieces) in target_evidence.items():
        other_state, other_context, other_pieces = comparator_evidence.get(
            seed_id, ('unavailable_seed_result', None, None))
        if state != 'ok' or other_state != 'ok':
            missing.update({reason: 1 for reason in {state, other_state} - {'ok'}})
        elif context != other_context:
            missing['mismatched_candidate_context'] += 1
        else:
            paired += 1
            target_intervals.extend(pieces)
            comparator_intervals.extend(other_pieces)

    def length(intervals):
        return sum(b - a for _, a, b in merge_intervals(intervals))

    target_bp = length(target_intervals)
    comparator_bp = length(comparator_intervals)
    union_bp = length(target_intervals + comparator_intervals)
    intersection_bp = target_bp + comparator_bp - union_bp
    observed = intersection_bp / union_bp if union_bp else None
    expected = target['expected_seeds']
    if not target['candidate_universe_declared']:
        status = 'candidate_universe_undeclared'
    elif not expected:
        status = 'no_candidates'
    elif not comparator:
        status = 'comparator_unavailable'
    elif not paired:
        status = next(iter(missing)) if len(missing) == 1 else 'no_paired_coordinate_evidence'
    elif paired != expected:
        status = 'partial_candidate_evidence'
    elif not union_bp:
        status = 'both_empty'
    else:
        status = 'ok'
    return dict(upstream_expected_seeds=expected, upstream_paired_seeds=paired,
                upstream_similarity_status=status, upstream_pair_issue_counts=dict(missing),
                upstream_interval_jaccard=observed if status == 'ok' else None,
                upstream_observed_interval_jaccard=observed,
                upstream_target_bp=target_bp if paired else None,
                upstream_comparator_bp=comparator_bp if paired else None,
                upstream_union_bp=union_bp if paired else None,
                upstream_intersection_bp=intersection_bp if paired else None)


def upstream_similarity_rows(seed_rows, sample_rows, reference_version):
    """Decorate upstream sample counts and return per-sample directed version pairs.

    Input is normalized by concordance_report._upstream. Compare the union of
    clipped native seed calls across identical, paired input candidate windows.
    A Jaccard of one requires nonempty signal. A full-panel Jaccard additionally
    requires complete coordinate evidence for the declared panel. BP and the
    observed-only column describe the scope of paired coordinate evidence.
    """
    by_version = defaultdict(dict)
    for row in seed_rows:
        by_version[row['sample_id'], row['target_version']][row['seed_id']] = row

    by_sample = defaultdict(dict)
    views = {}
    for row in sample_rows:
        sid, version = row['sample_id'], row['target_version']
        records = by_version[sid, version]
        row.update(_count_fields(row, records.values()))
        by_sample[sid][version] = row
        views[sid, version] = {seed: _evidence(record) for seed, record in records.items()
                               if row['expected_seeds']}

    pairs = []
    for sid, versions in by_sample.items():
        for version, target in versions.items():
            for comparator_version in dict.fromkeys([*versions, reference_version]):
                comparator = versions.get(comparator_version)
                pair = _pair_footprints(target, comparator, views[sid, version],
                                        views.get((sid, comparator_version), {}))
                pair.update(sample_id=sid, target_version=version, comparator_version=comparator_version,
                            source_version=target['source_version'], reference_version=reference_version,
                            upstream_similarity_scope=FOOTPRINT_SCOPE,
                            upstream_target_cnvkit_versions=target['cnvkit_versions'],
                            upstream_comparator_cnvkit_versions=comparator['cnvkit_versions'] if comparator else [],
                            upstream_retained_count=target['upstream_retained_count'],
                            upstream_evaluable_count=target['upstream_evaluable_count'],
                            upstream_retained_ratio=target['upstream_retained_ratio'],
                            upstream_panel_retained_ratio=target['upstream_panel_retained_ratio'],
                            upstream_unevaluable_count=target['upstream_unevaluable_count'],
                            upstream_count_status=target['upstream_count_status'])
                pairs.append(pair)
    return pairs


def upstream_similarity_cohort_rows(pairs, patient_ids, alpha=.05):
    """Pointwise descriptive summaries; counts remain separate from equally weighted means.

    Full-panel Jaccard requires a nonempty union. Retention includes complete
    panels regardless of signal. Rows with missing inputs still contribute to
    expected sample/seed totals. Within each source-version stratum a patient
    receives one mean of its evaluable samples before averaging across patients.
    """
    from .concordance_report import _bounds
    from .statistics import _independent_means

    if not 0 < alpha < 1:
        raise ValueError('alpha must be between zero and one')
    groups = defaultdict(list)
    for row in pairs:
        key = (row['reference_version'], row['source_version'], row['target_version'], row['comparator_version'])
        groups[key].append(row)
    summaries = []
    for (ref, source, target, comparator), rows in groups.items():
        samples = {row['sample_id'] for row in rows}
        if len(samples) != len(rows):
            raise ValueError('Duplicate upstream sample within a version comparison')
        ids = samples if patient_ids is None else {patient_ids[s] for s in samples}
        summary = dict(reference_version=ref, source_version=source, target_version=target,
                       comparator_version=comparator, identity_anchor=target == comparator,
                       n_selected_samples=len(samples), n_selected_independent_units=len(ids),
                       independence_unit='sample' if patient_ids is None else 'patient',
                       interpretation='descriptive', interval_method='pointwise_two_sided_Hoeffding',
                       confidence_level=1-alpha, upstream_similarity_scope=FOOTPRINT_SCOPE,
                       upstream_similarity_status_counts=dict(Counter(row['upstream_similarity_status'] for row in rows)),
                       upstream_count_status_counts=dict(Counter(row['upstream_count_status'] for row in rows)))
        for key in ('upstream_target_cnvkit_versions', 'upstream_comparator_cnvkit_versions'):
            summary[key] = sorted({version for row in rows for version in row[key]})
        for key in ('upstream_expected_seeds', 'upstream_evaluable_count', 'upstream_retained_count',
                    'upstream_unevaluable_count', 'upstream_paired_seeds'):
            summary[key] = sum(row[key] for row in rows)
        for metric, status_key, required in (
                ('upstream_interval_jaccard', 'upstream_similarity_status', 'ok'),
                ('upstream_retained_ratio', 'upstream_count_status', 'complete')):
            by_sample = {row['sample_id']: row[metric] for row in rows
                         if row[status_key] == required and row[metric] is not None}
            values = _independent_means(by_sample, patient_ids).tolist()
            # An identity comparison supplies a plotting anchor rather than a sampling interval.
            low, high = (None, None) if target == comparator and metric == 'upstream_interval_jaccard' else _bounds(values, alpha)
            summary.update({metric: sum(values) / len(values) if values else None,
                            metric + '_ci_low': low, metric + '_ci_high': high,
                            metric + '_n_samples': len(by_sample), metric + '_n_independent_units': len(values)})
        summaries.append(summary)
    return summaries
