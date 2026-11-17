"""Paired, sample-clustered comparisons. No inference of biological truth.

Each sample contributes one mean per metric, regardless of its seed count.
Compatibility uses a conservative bounded-mean agreement test, not a failure
to reject a difference. See docs/statistical_methods.md.
"""
from __future__ import annotations

import itertools
import math
from collections import defaultdict

import numpy as np

VALID = {"ok", "not_detected"}
AGREEMENT_METRICS = ("detection_discordance", "label_discordance", "interval_distance", "class_footprint_distance")


def merge_intervals(intervals):
    groups = defaultdict(list)
    for chrom, start, end in intervals:
        start, end = int(start), int(end)
        if start < 0 or end <= start:
            raise ValueError(f"Invalid BED interval: {chrom}:{start}-{end}")
        groups[str(chrom)].append((start, end))
    result = []
    for chrom in sorted(groups):
        for start, end in sorted(groups[chrom]):
            if result and result[-1][0] == chrom and start <= result[-1][2]:
                result[-1][2] = max(result[-1][2], end)
            else:
                result.append([chrom, start, end])
    return result


def subtract_intervals(first, second):
    second = merge_intervals(second)
    result = []
    for chrom, start, end in merge_intervals(first):
        cursor = start
        for c2, s2, e2 in second:
            if c2 != chrom or e2 <= cursor or s2 >= end:
                continue
            if s2 > cursor:
                result.append([chrom, cursor, min(s2, end)])
            cursor = max(cursor, e2)
            if cursor >= end:
                break
        if cursor < end:
            result.append([chrom, cursor, end])
    return result


def interval_jaccard(first, second):
    a, b = merge_intervals(first), merge_intervals(second)
    union = sum(e - s for _, s, e in merge_intervals(a + b))
    if not union:
        return None  # Both empty is not affirmative structural agreement.
    unique_a = sum(e - s for _, s, e in subtract_intervals(a, b))
    intersection = sum(e - s for _, s, e in a) - unique_a
    return intersection / union


def _detected(record):
    return bool(record.get("detected", record["status"] == "ok"))


def _feature_footprints(record):
    """Never infer class locations from an unlabelled overall interval union."""
    if record['status']=='not_detected' and not record.get('labels') and not record.get('intervals'):
        return ({},'') if not record.get('features') else (None,'negative_record_contains_features')
    features=record.get('features')
    if not isinstance(features,list) or not features:
        return None,'positive_record_missing_feature_intervals'
    footprints=defaultdict(list)
    try:
        for feature in features:
            label=feature['label']
            if not isinstance(label,str) or not feature.get('intervals'):
                return None,'feature_label_or_intervals_missing'
            footprints[label].extend(merge_intervals(feature['intervals']))
        footprints={label:merge_intervals(intervals) for label,intervals in footprints.items()}
        if set(footprints)!=set(record.get('labels',[])):
            return None,'feature_labels_differ_from_record_labels'
        union=merge_intervals([i for intervals in footprints.values() for i in intervals])
        if union!=merge_intervals(record.get('intervals',[])):
            return None,'feature_union_differs_from_record_intervals'
    except (TypeError,ValueError,KeyError):
        return None,'malformed_feature_intervals'
    return footprints,''


def _class_footprint_comparison(a,b,common_supported):
    if common_supported is None:
        return None,[],'category_support_unknown'
    labels=sorted((set(a.get('labels',[]))|set(b.get('labels',[])))&set(common_supported))
    if not labels:
        return None,[],'no_common_supported_positive_labels'
    first,error_a=_feature_footprints(a)
    second,error_b=_feature_footprints(b)
    details=[]
    if first is None or second is None:
        reason='; '.join(x for x in (('candidate: '+error_a) if error_a else '',('comparator: '+error_b) if error_b else '') if x)
        return None,[dict(label=label,status='unavailable',jaccard=None,distance=None,reason=reason) for label in labels],reason
    for label in labels:
        ia,ib=first.get(label,[]),second.get(label,[])
        jac=interval_jaccard(ia,ib)
        details.append(dict(label=label,status='complete',candidate_intervals=ia,comparator_intervals=ib,
            intervals_added=subtract_intervals(ia,ib),intervals_removed=subtract_intervals(ib,ia),
            jaccard=jac,distance=1-jac if jac is not None else None,reason=''))
    if any(r['jaccard'] is None for r in details):
        return None,details,'positive_label_has_empty_feature_footprint'
    return float(np.mean([r['jaccard'] for r in details])),details,''


def _pair(a, b, reference_version, pair_kind):
    """a is candidate, b comparator; positive differences mean extra in a."""
    valid = a["status"] in VALID and b["status"] in VALID
    la, lb = set(a.get("labels", [])), set(b.get("labels", []))
    sa, sb = a.get("supported_labels"), b.get("supported_labels")
    scope = sorted(set(sa) & set(sb)) if sa is not None and sb is not None else None
    full_labels = sa is not None and sb is not None and set(sa) == set(sb)
    row = {
        "reference_version": reference_version, "pair_kind": pair_kind,
        "version": a["version"], "comparator": b["version"],
        "sample_id": a["sample_id"], "seed_id": a["seed_id"],
        "reference_scope": a.get("reference_scope", "declared_reference"),
        "comparator_reference_scope": b.get("reference_scope", "declared_reference"),
        "status": a["status"], "comparator_status": b["status"],
        "is_seed_region": a['status']!='no_seeds' and b['status']!='no_seeds' and a['seed_id']!='__no_seeds__',
        "evaluable": valid, "signal_positive": False,
        "labels": sorted(la), "comparator_labels": sorted(lb),
        "label_scope": scope if scope is not None else "as_reported",
        "full_label_comparable": bool(full_labels),
        "labels_added": sorted(la - lb), "labels_removed": sorted(lb - la),
        "intervals": a.get("intervals", []), "comparator_intervals": b.get("intervals", []),
        "stage": a.get("stage", ""), "reason": a.get("reason", ""),
        "comparator_stage": b.get("stage", ""), "comparator_reason": b.get("reason", ""),
        "detection_stage": a.get('detection_stage',''), "comparator_detection_stage": b.get('detection_stage',''),
        "aa_executed": a.get('aa_executed'), "comparator_aa_executed": b.get('aa_executed'),
        "aa_reconstruction_performed": a.get('aa_reconstruction_performed'),
        "comparator_aa_reconstruction_performed": b.get('aa_reconstruction_performed'),
        "ac_executed": a.get('ac_executed'), "comparator_ac_executed": b.get('ac_executed'),
        "seed_count": a.get('seed_count'), "comparator_seed_count": b.get('seed_count'),
        "seed_count_difference": None,"seed_interval_jaccard": None,"seed_interval_distance": None,
        "seed_intervals_added": [],"seed_intervals_removed": [],
        "class_footprint_jaccard":None,"class_footprint_distance":None,
        "class_footprint_details":[],"class_footprint_reason":"not_evaluable",
        "changes": [], "intervals_added": [], "intervals_removed": [],
    }
    for field in ("detection_difference", "detection_discordance", "detection_agreement", "extra_detection", "missed_detection", "label_discordance", "class_change", "interval_jaccard", "interval_distance"):
        row[field] = None
    if a.get('seed_count') is not None and b.get('seed_count') is not None:
        row['seed_count_difference']=int(a['seed_count'])-int(b['seed_count'])
    if 'seed_intervals' in a and 'seed_intervals' in b:
        row['seed_interval_jaccard']=interval_jaccard(a['seed_intervals'],b['seed_intervals'])
        row['seed_interval_distance']=1-row['seed_interval_jaccard'] if row['seed_interval_jaccard'] is not None else None
        row['seed_intervals_added']=subtract_intervals(a['seed_intervals'],b['seed_intervals'])
        row['seed_intervals_removed']=subtract_intervals(b['seed_intervals'],a['seed_intervals'])
        if row['seed_intervals_added'] or row['seed_intervals_removed']:
            row['changes'].append('seed_sets_changed')
    if not valid:
        row["changes"].append("no_input_seeds" if not row['is_seed_region'] else "not_evaluable")
        return row
    da, db = _detected(a), _detected(b)
    if da or db:
        jac,details,reason=_class_footprint_comparison(a,b,scope)
        row.update(class_footprint_jaccard=jac,class_footprint_distance=1-jac if jac is not None else None,
                   class_footprint_details=details,class_footprint_reason=reason)
        if jac is None:
            row['changes'].append('class_footprint_information_unavailable')
        elif jac<1:
            row['changes'].append('class_footprints_changed')
    row.update(signal_positive=da or db, detection_difference=int(da)-int(db),
               detection_discordance=float(da != db), detection_agreement=float(da == db),
               extra_detection=float(da and not db), missed_detection=float(db and not da))
    # Unknown support must never be silently coerced into a negative category.
    comparable_a, comparable_b = (la & set(scope), lb & set(scope)) if scope is not None else (la, lb)
    row["label_discordance"] = float(comparable_a != comparable_b) if full_labels or scope else None
    row["class_change"] = float(comparable_a != comparable_b) if da and db and (full_labels or scope) else None
    if da != db:
        row["changes"].append("extra_detection" if da else "missed_detection")
    if comparable_a != comparable_b:
        row["changes"].append("classification_changed")
    ia, ib = a.get("intervals", []), b.get("intervals", [])
    intervals_valid = (not da or bool(ia)) and (not db or bool(ib))
    if intervals_valid and (da or db):
        row["interval_jaccard"] = interval_jaccard(ia, ib)
        row["interval_distance"] = 1 - row["interval_jaccard"] if row["interval_jaccard"] is not None else None
        row["intervals_added"], row["intervals_removed"] = subtract_intervals(ia, ib), subtract_intervals(ib, ia)
        if row["intervals_added"] or row["intervals_removed"]:
            row["changes"].append("intervals_changed")
    elif da or db:
        row["changes"].append("interval_information_missing")
    if not full_labels:
        row["changes"].append("category_support_differs_or_unknown")
    if not row["changes"]:
        row["changes"] = ["no_observed_change"]
    return row


def _sample_means(rows, metric):
    by_sample = defaultdict(list)
    for row in rows:
        if row.get(metric) is not None:
            by_sample[row["sample_id"]].append(float(row[metric]))
    return np.array([np.mean(by_sample[s]) for s in sorted(by_sample)], dtype=float), sum(map(len, by_sample.values()))


def _bootstrap(values, iterations, rng, alpha):
    if len(values) < 2:
        return None, None
    draws = []
    # Bound memory even when called with larger cohorts.
    for start in range(0, iterations, 500):
        inds = rng.integers(0, len(values), size=(min(500, iterations-start), len(values)))
        draws.extend(np.mean(values[inds], axis=1))
    lo, hi = np.quantile(draws, [alpha/2, 1-alpha/2])
    return float(lo), float(hi)


def _sign_flip(values, iterations, rng):
    """Paired whole-sample label-swap test; requires exchangeability under H0."""
    values = values[values != 0]
    if not len(values):
        return 1.0, "exact_cluster_sign_flip"
    observed = abs(float(values.sum()))
    if len(values) <= 16:
        extreme = sum(abs(sum(v*s for v, s in zip(values, signs))) >= observed-1e-12
                      for signs in itertools.product((-1, 1), repeat=len(values)))
        return extreme / (2 ** len(values)), "exact_cluster_sign_flip"
    extreme = 0
    for start in range(0, iterations, 500):
        signs = rng.choice([-1, 1], size=(min(500, iterations-start), len(values)))
        extreme += int(np.sum(np.abs(signs @ values) >= observed-1e-12))
    return (extreme+1)/(iterations+1), "monte_carlo_cluster_sign_flip"


def _holm(rows, source, target):
    usable = sorted((r for r in rows if r.get(source) is not None), key=lambda r: r[source])
    running = 0.0
    for i, row in enumerate(usable):
        running = max(running, min(1.0, row[source]*(len(usable)-i)))
        row[target] = running


def analyze(records, reference_version, options=None):
    """Return per-seed differences, metric estimates, and compatibility decisions.

    options: versions, margins, min_samples, alpha, bootstrap, permutations, seed.
    Margins are maximum tolerable mean losses in [0,1], chosen before analysis.
    All eligible seeds within each sample are averaged before inference.
    """
    options = options or {}
    analysis_scope=options.get('analysis_scope','fixed_seed')
    if analysis_scope not in ('fixed_seed','full_suite'):
        raise ValueError('analysis_scope must be fixed_seed or full_suite')
    if analysis_scope=='full_suite' and any(r['seed_id']!='__whole_sample__' for r in records):
        raise ValueError('full_suite requires exactly one __whole_sample__ record per sample/version, including failures')
    alpha = float(options.get("alpha", .05))
    min_samples = int(options.get("min_samples", 20))
    bootstrap = int(options.get("bootstrap", 2000))
    permutations = int(options.get("permutations", 10000))
    if not 0 < alpha < 1 or min_samples < 2 or bootstrap < 100 or permutations < 100:
        raise ValueError("alpha must be (0,1); min_samples >=2; bootstrap/permutations >=100")
    margins = options.get("margins", {})
    for name, value in margins.items():
        if name not in AGREEMENT_METRICS or value is not None and not 0 < float(value) < 1:
            raise ValueError(f"Unknown/invalid agreement margin: {name}={value}")
    versions = list(dict.fromkeys(options.get("versions", []) + [r["version"] for r in records]))
    if reference_version not in versions:
        raise ValueError(f"Reference version absent: {reference_version}")
    keys = sorted({(r["sample_id"], r["seed_id"]) for r in records})
    lookup = {}
    for r in records:
        key = (r["version"], r["sample_id"], r["seed_id"])
        if key in lookup:
            raise ValueError(f"Duplicate normalized record: {key}")
        if r["status"] not in VALID | {"failed", "missing", "unsupported", "invalid", "interrupted", "blocked", "no_seeds", "not_run"}:
            raise ValueError(f"Unknown status {r['status']}")
        lookup[key] = r
    no_seed_samples={r['sample_id'] for r in records if r['status']=='no_seeds'}
    if any(r['sample_id'] in no_seed_samples and r['seed_id']!='__no_seeds__' for r in records):
        raise ValueError('Sample has both no_seeds sentinel and actual seed records; seed universe inconsistent')
    def get(version, sample, seed):
        return lookup.get((version, sample, seed), dict(version=version, sample_id=sample, seed_id=seed,
                         status="missing", reason="No normalized record", stage="result_reading"))
    pairs = [(v, reference_version, "vs_reference") for v in versions if v != reference_version]
    historical = [v for v in versions if v != reference_version]
    pairs += [(a, b, "historical_pair") for a, b in itertools.combinations(historical, 2)]
    rng = np.random.default_rng(int(options.get("seed", 20261117)))
    pair_rows, stats_rows, compat_rows = [], [], []
    # Simultaneous one-sided limits across all version pairs and four metrics.
    decision_alpha = alpha/max(1, len(pairs)*len(AGREEMENT_METRICS))
    for version, comparator, kind in pairs:
        rows = [_pair(get(version,s,t), get(comparator,s,t), reference_version, kind) for s,t in keys]
        for row in rows:
            row['analysis_scope']=analysis_scope
            row['comparison_unit']='sample' if analysis_scope=='full_suite' else 'sample_seed'
        pair_rows.extend(rows)
        region_rows=[r for r in rows if r['is_seed_region']]
        n_no_seed_samples=len({r['sample_id'] for r in rows if not r['is_seed_region']})
        valid = [r for r in rows if r["evaluable"]]
        signal = [r for r in valid if r["signal_positive"]]
        local_stats = []
        main_metrics=("detection_difference", "detection_agreement", "detection_discordance", "extra_detection",
                      "missed_detection", "label_discordance", "class_change", "interval_jaccard", "interval_distance",
                      "class_footprint_jaccard","class_footprint_distance")
        selections=[("all_evaluable",valid,main_metrics),("signal_positive",signal,main_metrics)]
        if analysis_scope=='full_suite':
            selections.append(('seed_stage_evaluable',rows,('seed_count_difference','seed_interval_jaccard','seed_interval_distance')))
        for subset, selected,metrics in selections:
            for metric in metrics:
                values, n_regions = _sample_means(selected, metric)
                n = len(values)
                estimate = float(np.mean(values)) if n else None
                lo, hi = _bootstrap(values, bootstrap, rng, alpha) if n else (None, None)
                lower, upper = (-1, 1) if metric == "detection_difference" else (0, 1)
                radius = (upper-lower)*math.sqrt(math.log(2/alpha)/(2*n)) if n and metric!='seed_count_difference' else None
                p, test = _sign_flip(values, permutations, rng) if n >= 2 and metric in ("detection_difference","seed_count_difference") else (None, "")
                event_metrics = {"detection_agreement", "detection_discordance", "extra_detection", "missed_detection", "label_discordance", "class_change"}
                raw_events = int(sum(r[metric] for r in selected if r[metric] is not None)) if n_regions and metric in event_metrics else None
                stat = dict(reference_version=reference_version, version=version, comparator=comparator, pair_kind=kind,
                            analysis_scope=analysis_scope,comparison_unit='sample' if analysis_scope=='full_suite' else 'sample_seed',
                            subset=subset, metric=metric, estimate=estimate, bootstrap_ci_low=lo, bootstrap_ci_high=hi,
                            conservative_ci_low=max(lower, estimate-radius) if radius is not None else None,
                            conservative_ci_high=min(upper, estimate+radius) if radius is not None else None,
                            confidence_level=1-alpha, n_samples=n, n_regions=n_regions,n_comparison_units=n_regions,raw_region_events=raw_events,
                            total_regions=len(region_rows), excluded_regions=len(region_rows)-len(valid),no_seed_samples=n_no_seed_samples,
                            p_difference=p, p_difference_holm=None, difference_test=test,
                            inference="sample_equal_weight", bootstrap_degenerate=bool(n and lo == hi),
                            min_samples=min_samples, sufficient_samples=n >= min_samples,
                            margin=None, agreement_upper=None, p_agreement_bound=None, agreement_decision="not_assessed")
                if metric in AGREEMENT_METRICS and subset == "signal_positive":
                    margin = margins.get(metric)
                    u = min(1.0, estimate+math.sqrt(math.log(1/decision_alpha)/(2*n))) if n else None
                    p_agree = math.exp(-2*n*(float(margin)-estimate)**2) if n and margin is not None and estimate < float(margin) else 1.0 if n and margin is not None else None
                    decision = "insufficient_evidence"
                    if n >= min_samples and margin is not None:
                        if u < float(margin):
                            decision = "within_margin"
                        elif max(0., estimate-math.sqrt(math.log(1/decision_alpha)/(2*n))) > float(margin):
                            decision = "exceeds_margin"
                    stat.update(margin=margin, agreement_upper=u, p_agreement_bound=p_agree,
                                agreement_alpha=decision_alpha, agreement_decision=decision)
                stats_rows.append(stat)
                local_stats.append(stat)
        primary = [r for r in local_stats if r["subset"] == "signal_positive" and r["metric"] in AGREEMENT_METRICS]
        reasons = []
        if any(margins.get(m) is None for m in AGREEMENT_METRICS):
            reasons.append("agreement_margins_not_predefined")
        if any(r["n_samples"] < min_samples for r in primary):
            reasons.append("too_few_independent_samples")
        if len(valid) != len(region_rows):
            reasons.append("unevaluable_regions_present")
        if n_no_seed_samples:
            reasons.append('samples_without_input_seeds')
        if any(not r["full_label_comparable"] for r in valid):
            reasons.append("category_support_differs_or_unknown")
        if any(r["interval_distance"] is None for r in signal):
            reasons.append("interval_information_missing")
        if any(r['class_footprint_distance'] is None for r in signal):
            reasons.append('class_footprint_information_missing')
        scopes = sorted({r[key] for r in rows for key in ("reference_scope", "comparator_reference_scope")})
        if any(r["reference_scope"] != r["comparator_reference_scope"] for r in valid):
            reasons.append("reference_scopes_differ")
        decision = "insufficient_evidence"
        if not reasons and all(r["agreement_decision"] == "within_margin" for r in primary):
            decision = "within_margin"
        elif any(r["agreement_decision"] == "exceeds_margin" for r in primary):
            decision = "exceeds_margin"
        else:
            reasons.append("confidence_bounds_do_not_establish_agreement")
        compat_rows.append(dict(reference_version=reference_version, version=version, comparator=comparator,
                                analysis_scope=analysis_scope,comparison_unit='sample' if analysis_scope=='full_suite' else 'sample_seed',
                                pair_kind=kind, decision=decision,
                                conclusion={"within_margin":"在预设检出、分类及区段指标范围内支持初步兼容", "exceeds_margin":"至少一个指标超过预设差异范围", "insufficient_evidence":"证据不足"}[decision],
                                reasons=sorted(set(reasons)), total_regions=len(region_rows), evaluable_regions=len(valid),
                                signal_regions=len(signal), signal_samples=len({r['sample_id'] for r in signal}),
                                excluded_regions=len(region_rows)-len(valid), no_seed_samples=n_no_seed_samples,min_samples=min_samples,
                                reference_scopes=scopes,
                                family_alpha=alpha, per_metric_alpha=decision_alpha, margins=margins,
                                compatibility_metric_scope={'required_metrics':list(AGREEMENT_METRICS),
                                    'not_established':['same_cycle_order','same_breakpoint_topology','biological_truth','whole_version_interchangeability']},
                                scope=("tested samples; full Suite including CNVkit, seed processing, AA and AC; whole-sample label/feature unions, not amplicon-ID matching; " if analysis_scope=='full_suite' else "tested samples, fixed seeds, shared resources, AA+AC; CNVkit excluded; ")+"reference scope="+",".join(scopes)))
    _holm(stats_rows, "p_difference", "p_difference_holm")
    return pair_rows, stats_rows, compat_rows


def build_sample_similarity(records, reference_version, pair_rows, analysis_scope='fixed_seed'):
    """Descriptive per-input-sample similarity, with explicit missing/status scope.

    No confidence interval or equivalence decision is made for one sample.
    In fixed-seed mode the displayed metric is a mean across evaluable seeds;
    full-suite mode has one whole-sample record per version.
    """
    versions=list(dict.fromkeys([reference_version]+[r['version'] for r in records]+[r['version'] for r in pair_rows]))
    samples=sorted({r['sample_id'] for r in records}|{r['sample_id'] for r in pair_rows})
    by_sample=defaultdict(list)
    for r in records:by_sample[r['sample_id'],r['version']].append(r)
    result=[]
    def mean(rows,key):
        values=[r[key] for r in rows if r.get(key) is not None]
        return float(np.mean(values)) if values else None
    def one(rows,key):
        values={r.get(key) for r in rows if r.get(key) is not None}
        return next(iter(values)) if len(values)==1 else None
    for sample in samples:
        reference=by_sample[sample,reference_version]
        for version in versions:
            source=by_sample[sample,version]
            compared=([_pair(r,r,reference_version,'identity_anchor') for r in reference] if version==reference_version else
                      [r for r in pair_rows if r['sample_id']==sample and r['version']==version and r['comparator']==reference_version])
            regions=[r for r in compared if r['is_seed_region']]
            valid=[r for r in regions if r['evaluable']]
            statuses=sorted({r['status'] for r in source}) or ['missing']
            reference_statuses=sorted({r['status'] for r in reference}) or ['missing']
            if statuses==['no_seeds'] or any(not r['is_seed_region'] for r in compared):
                state='no_seeds'
            elif not valid:
                state='not_evaluable'
            elif len(valid)<len(regions):
                state='partially_evaluable'
            elif all(not r['signal_positive'] for r in valid):
                state='both_not_detected'
            else:
                state='evaluable'
            label_discordance=mean(valid,'label_discordance')
            result.append(dict(sample_id=sample,version=version,reference_version=reference_version,
                analysis_scope=analysis_scope,comparison_unit='sample' if analysis_scope=='full_suite' else 'sample_seed',
                comparison_state=state,observed_statuses=statuses,reference_statuses=reference_statuses,
                identity_anchor=version==reference_version,n_input_units=len(regions),n_evaluable_units=len(valid),
                n_signal_units=sum(r['signal_positive'] for r in valid),
                detection_agreement=mean(valid,'detection_agreement'),
                label_set_agreement=1-label_discordance if label_discordance is not None else None,
                interval_jaccard=mean(valid,'interval_jaccard'),
                class_footprint_jaccard=mean(valid,'class_footprint_jaccard'),
                class_footprint_distance=mean(valid,'class_footprint_distance'),
                full_label_comparable=all(r['full_label_comparable'] for r in valid) if valid else False,
                label_scopes=[r['label_scope'] for r in valid],
                labels=sorted({label for r in source if r['status'] in VALID for label in r.get('labels',[])}),
                reference_labels=sorted({label for r in reference if r['status'] in VALID for label in r.get('labels',[])}),
                seed_count=one(source,'seed_count'),reference_seed_count=one(reference,'seed_count'),
                seed_count_difference=mean(compared,'seed_count_difference'),seed_interval_jaccard=mean(compared,'seed_interval_jaccard'),
                detection_stages=sorted({r.get('detection_stage','') for r in source}),
                interpretation='descriptive_observed_similarity_not_equivalence',
                uncertainty='no_single_sample_confidence_interval'))
    return result
