"""Reports for cross-version AmpliconSuite concordance; no BAM execution or export.

Input-source comparisons, fixed-reference comparisons, and AC-only experiments
remain distinct. Inference uses a patient cluster when supplied, or a sample.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from .concordance_global_report import GLOBAL_METRIC, GLOBAL_SCOPE, LOCUS_RULE
from .reporting import write_tsv
from .statistics import (AGREEMENT_METRICS, VALID, _independent_means,
                         _is_primary_pair, _pair, _primary_pair_keys, analyze)


KINDS = {"aa_ac": "AA + AC", "ac_only": "AC only on fixed source AA output"}
DECISIONS = {"within_margin": "within margins", "exceeds_margin": "exceeds margins",
             "insufficient_evidence": "insufficient evidence"}
SIMILARITIES = (("detection_agreement", "Detection agreement"),
                ("label_agreement", "Label-set agreement"),
                ("interval_jaccard", "Interval Jaccard"),
                ("class_footprint_jaccard", "Class-footprint Jaccard"))


def _kind_label(kind, mode):
    if kind == "ac_only": return KINDS[kind]
    return "selected-seed AA + AC" if mode == "with_bam" else "existing AA + AC output metrics"


def _slug(value):
    import re
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")[:90] + "_" + hashlib.sha256(value.encode()).hexdigest()[:8]


def _json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _mean(values):
    return sum(values) / len(values) if values else None


def _bounds(values, alpha, two_sided=True, limits=(0., 1.)):
    if not values:
        return None, None
    estimate = _mean(values)
    lower, upper = limits
    radius = (upper-lower)*math.sqrt(math.log((2 if two_sided else 1) / alpha) / (2 * len(values)))
    return max(lower, estimate-radius), min(upper, estimate+radius)


def _label_scope(rows):
    scopes = []
    for row in rows:
        if row["label_scope"] not in scopes: scopes.append(row["label_scope"])
    known = bool(scopes) and all(isinstance(scope, list) for scope in scopes)
    common = sorted(set.intersection(*(set(scope) for scope in scopes))) if known else None
    return dict(full_label_comparable=bool(rows) and all(row["full_label_comparable"] for row in rows),
                common_supported_labels=common, label_scopes=scopes)


def _configuration(config):
    ref = config.get("reference_version", "suite-1.6.2")
    targets = list(dict.fromkeys(config.get("target_versions", [])))
    if not targets or any(not isinstance(v, str) or not v for v in targets):
        raise ValueError("target_versions must contain explicit Suite version IDs")
    if not isinstance(ref, str) or not ref:
        raise ValueError("reference_version must be fixed explicitly")
    if config.get("mode") not in {"with_bam", "results_only"}:
        raise ValueError("mode must be with_bam or results_only")
    samples = {}
    for sample in config.get("samples", []):
        sid, source = sample.get("sample_id"), sample.get("source_version")
        if not isinstance(sid, str) or not sid or not isinstance(source, str) or not source or sid in samples:
            raise ValueError("Each sample requires a unique sample_id and declared source_version")
        if "upstream_seed_ids" in sample:
            ids = sample["upstream_seed_ids"]
            if not isinstance(ids, list) or any(not isinstance(v, str) or not v for v in ids) or len(ids) != len(set(ids)):
                raise ValueError("upstream_seed_ids must be an explicit unique list, or [] for no input seeds")
        if "selected_seed_ids" in sample:
            ids = sample["selected_seed_ids"]
            if not isinstance(ids, list) or any(not isinstance(v, str) or not v for v in ids) or len(ids) != len(set(ids)):
                raise ValueError("selected_seed_ids must be a unique list of expected downstream seed IDs")
        if sample.get("source_bridge_status", "unverified") not in {"passed", "failed", "unverified"}:
            raise ValueError("source_bridge_status must be passed, failed or unverified")
        samples[sid] = dict(sample)
    if not samples:
        raise ValueError("The selected sample manifest is required, including samples with no results")
    if any(sample.get('patient_id') is not None for sample in samples.values()):
        patient_ids = {sid: sample.get('patient_id') for sid, sample in samples.items()}
        if any(not isinstance(value, str) or not value.strip() for value in patient_ids.values()):
            raise ValueError('Provide a nonempty patient_id for every sample when using patient clustering')
    else:
        patient_ids = None
    options = dict(config.get("statistics", {}))
    primary_keys = _primary_pair_keys(options.get("primary_pairs"))
    if isinstance(primary_keys, set):
        known_versions = {ref, *targets, *(sample['source_version'] for sample in samples.values())}
        if any(pair - known_versions for pair in primary_keys):
            raise ValueError('primary_pairs contains a version that is not in the selected comparison design')
    options['patient_ids'] = patient_ids
    alpha, minimum = float(options.get("alpha", .05)), int(options.get("min_samples", 20))
    margin = options.get("upstream_loss_margin")
    if not 0 < alpha < 1 or minimum < 2:
        raise ValueError("alpha must be (0,1) and min_samples >=2")
    if margin is not None and (not math.isfinite(float(margin)) or not 0 <= float(margin) <= 1):
        raise ValueError("upstream_loss_margin must be a predefined proportion in [0,1], or null")
    return ref, targets, samples, options, alpha, minimum, margin


def _prepare_downstream(records, samples, targets, ref, kinds, mode):
    """Preserve observed status separately when comparability has not been verified."""
    lookup, units = {}, defaultdict(set)
    for original in records:
        row = dict(original)
        sample, version, kind = row.get("sample_id"), row.get("version"), row.get("analysis_kind")
        if sample not in samples or kind not in KINDS or not version:
            raise ValueError("Unknown sample, analysis_kind or missing version in downstream record")
        unit = row.get("seed_id") or row.get("comparison_unit")
        if not isinstance(unit, str) or not unit:
            raise ValueError("Downstream records require a stable seed_id/comparison_unit identifier")
        source = samples[sample]["source_version"]
        if row.get("source_version", source) != source:
            raise ValueError("Record source_version differs from the sample manifest")
        key = (kind, sample, unit, version)
        if key in lookup:
            raise ValueError(f"Duplicate downstream comparison unit: {key}")
        row.update(seed_id=unit, source_version=source, observed_status=row.get("status", "missing"))
        row.setdefault("status", "missing")
        row.setdefault("labels", [])
        row.setdefault("intervals", [])
        if row.get("comparison_eligible") is not True and row["status"] in VALID:
            row.update(status="unsupported", stage="comparison_input", reason="Same-input comparison eligibility not explicitly verified; " + row.get("reason", ""))
        lookup[key] = row
        units[kind, sample].add(unit)
    groups = {}
    for kind in kinds:
        for source in sorted({s["source_version"] for s in samples.values()}):
            group = []
            versions = list(dict.fromkeys([ref, source] + targets))
            for sample in sorted(s for s in samples if samples[s]["source_version"] == source):
                if mode == "with_bam" and kind == "aa_ac":
                    units[kind, sample].update(samples[sample].get("selected_seed_ids", []))
                for unit in sorted(units[kind, sample] or {"__missing_input__"}):
                    local = []
                    for version in versions:
                        row = dict(lookup.get((kind, sample, unit, version), dict(
                            sample_id=sample, version=version, seed_id=unit, status="missing", observed_status="missing",
                            analysis_kind=kind, source_version=source, comparison_eligible=False,
                            labels=[], intervals=[], stage="comparison_input", reason="No result for selected sample/unit/version")))
                        local.append(row)
                    if kind == "ac_only":
                        present = [r for r in local if r["status"] in VALID]
                        aa_versions = {r.get("source_aa_version") for r in present}
                        graph_ids = {r.get("source_graph_identity") for r in present}
                        if None in aa_versions or "" in aa_versions or len(aa_versions) > 1:
                            for row in present:
                                row.update(status="unsupported", reason="AC-only fixed AA source version missing or inconsistent", stage="comparison_input")
                        elif None in graph_ids or "" in graph_ids or len(graph_ids) > 1:
                            for row in present:
                                row.update(status="unsupported", reason="AC-only fixed AA graph identity missing or inconsistent", stage="comparison_input")
                    group.extend(local)
            groups[kind, source] = group
    return groups


def _oriented_stat(stats, target, source, metric, subset='signal_positive'):
    direct = stats.get((target, source, subset, metric))
    if direct is not None:
        return direct
    reverse_metric = {'extra_detection': 'missed_detection', 'missed_detection': 'extra_detection'}.get(metric, metric)
    reverse = stats.get((source, target, subset, reverse_metric))
    if reverse is None or metric != 'detection_difference':
        return reverse
    # The stored pair order is relative to the fixed panel reference. Reverse
    # both the sign and interval endpoints for target versus input source.
    oriented = dict(reverse)
    if reverse['estimate'] is not None:
        oriented['estimate'] = -reverse['estimate']
    for interval in ('conservative', 'bootstrap'):
        low, high = (reverse[f'{interval}_ci_{bound}'] for bound in ('low', 'high'))
        oriented[f'{interval}_ci_low'] = -high if high is not None else None
        oriented[f'{interval}_ci_high'] = -low if low is not None else None
    return oriented


def _detection_change_fields(stats, target, source):
    delta = _oriented_stat(stats, target, source, 'detection_difference', 'all_evaluable') or {}
    extra = _oriented_stat(stats, target, source, 'extra_detection', 'all_evaluable') or {}
    missed = _oriented_stat(stats, target, source, 'missed_detection', 'all_evaluable') or {}
    classes = _oriented_stat(stats, target, source, 'class_change', 'all_evaluable') or {}
    return dict(detection_difference=delta.get('estimate'),
                detection_difference_ci_low=delta.get('conservative_ci_low'),
                detection_difference_ci_high=delta.get('conservative_ci_high'),
                detection_difference_n_samples=delta.get('n_samples', 0),
                detection_difference_n_independent_units=delta.get('n_independent_units', 0),
                detection_difference_n_regions=delta.get('n_regions', 0),
                expected_regions=delta.get('total_regions'),
                extra_detection_count=extra.get('raw_region_events'),
                missed_detection_count=missed.get('raw_region_events'),
                class_change_count=classes.get('raw_region_events'),
                class_change_n_regions=classes.get('n_regions', 0))


def _downstream(out, groups, targets, ref, options, group_alpha, mode, synthetic, global_mode=False):
    comparisons, sample_rows, links = [], [], []
    for (kind, source), records in groups.items():
        folder = out / "downstream" / kind / _slug(source)
        folder.mkdir(parents=True, exist_ok=True)
        settings = dict(options, alpha=group_alpha, analysis_scope="fixed_seed",
                        versions=list(dict.fromkeys([ref, source] + targets)))
        settings.pop("upstream_loss_margin", None)
        pairs, stats, compatibility = analyze(records, ref, settings)
        stats_index = {(s['version'], s['comparator'], s['subset'], s['metric']): s for s in stats}
        # The shared calculation is valid for these endpoints, but its AA+AC
        # wording is replaced for AC-only output, never reused as an AA claim.
        label = _kind_label(kind, mode)
        scope = label + "; conditional on selected input units; no genome-wide seed-discovery validation"
        endpoint = "ac_only_fixed_source_aa" if kind == "ac_only" else "selected_seed_aa_ac" if mode == "with_bam" else "existing_aa_ac_output_metrics"
        unit = "sample_seed" if mode == "with_bam" and kind == "aa_ac" else "sample_input_result"
        if global_mode:
            label = "AA + AC on version-native seeds at sampled loci"
            scope = "observed global-discovery pipeline outputs at sampled loci; native AA seed inputs vary by version"
            endpoint, unit = "native_seed_aa_ac", "sample_locus"
        for table in (pairs, stats, compatibility):
            for row in table:
                row.update(analysis_kind=kind, analysis_scope=endpoint, comparison_unit=unit,
                           source_version=source, benchmark_scope=scope, synthetic=synthetic)
        for row in compatibility:
            row["scope"] = scope
            row["conclusion"] = {"within_margin": "在本项已验证指标和预设范围内支持初步比较",
                                 "exceeds_margin": "至少一个已评价指标超过预设范围", "insufficient_evidence": "证据不足"}[row["decision"]]
        for filename, table in (("records.tsv", records), ("differences.tsv", pairs), ("statistics.tsv", stats), ("compatibility.tsv", compatibility)):
            write_tsv(folder / filename, table)
        _json(folder / "records.json", records)
        (folder / "report.md").write_text(f"# {label}\n\n固定参考：`{ref}`；输入来源：`{source}`。\n\n"
            + ("**SYNTHETIC TEST DATA / 人工合成测试数据：仅验证程序，不形成患者或生物学结论。**\n\n" if synthetic else "")
            + ("仅重做 AC 分类，不能称为 AA 重建兼容。\n\n" if kind == "ac_only" else "各版本以本版原生种子重建，按抽中的共同位点比较完整结构。\n\n" if global_mode else "相同输入种子的 AA＋AC 比较，未验证全基因组种子发现。\n\n" if mode == "with_bam" else "既有 AA＋AC 输出指标比较；未证明历史参数和资源相同，差异不能自动归因于算法版本。\n\n")
            + "[统计表](statistics.tsv) · [差异明细](differences.tsv) · [兼容表](compatibility.tsv)。先对样本内的种子取平均。提供患者标识时，再对患者内样本取平均并按患者等权。单样本只描述观察相似度。\n", encoding="utf-8")
        links.append(dict(analysis_kind=kind, source_version=source, path=str((folder / "report.md").relative_to(out))))
        lookup = {(r["sample_id"], r["seed_id"], r["version"]): r for r in records}
        units = sorted({(r["sample_id"], r["seed_id"]) for r in records})
        for target in list(dict.fromkeys([source] + targets)):
            target_pairs = [_pair(lookup[s, unit, target], lookup[s, unit, source], ref, "vs_input_source") for s, unit in units]
            native = next((r for r in compatibility if {r["version"], r["comparator"]} == {target, source}), None) if target != source else None
            row = dict(reference_version=ref, source_version=source, target_version=target, analysis_kind=kind,
                       comparison_scope=scope, identity_anchor=target == source,
                       confidence_level=1-group_alpha,
                       primary_pair=native['primary_pair'] if native else False,
                       independence_unit='patient' if options['patient_ids'] is not None else 'sample',
                       decision=native["decision"] if native else "insufficient_evidence",
                       reasons=native["reasons"] if native else ["identity_is_not_validation"],
                       n_selected_samples=len({s for s, _ in units}), signal_samples=native["signal_samples"] if native else None,
                       eligible_units=native["evaluable_regions"] if native else None,
                       excluded_units=native["excluded_regions"] if native else None,
                       signal_regions=native['signal_regions'] if native else None,
                       signal_independent_units=native['signal_independent_units'] if native else None,
                       compatibility_metric_scope=native["compatibility_metric_scope"] if native else {"required_metrics": list(AGREEMENT_METRICS)})
            row.update(_label_scope([p for p in target_pairs if p["evaluable"]]))
            row.update(_detection_change_fields(stats_index, target, source))
            row['observed_statuses'] = sorted({p[k] for p in target_pairs for k in ('status', 'comparator_status')})
            row['decision_metrics'] = [dict(stat) for metric in AGREEMENT_METRICS
                if target != source and (stat := _oriented_stat(stats_index, target, source, metric)) is not None]
            for name, _ in SIMILARITIES:
                metric = "label_discordance" if name == "label_agreement" else name
                stat = _oriented_stat(stats_index, target, source, metric) if target != source else None
                value, lo, hi = ((stat[k] for k in ("estimate", "conservative_ci_low", "conservative_ci_high")) if stat else (None, None, None))
                if name == "label_agreement" and value is not None:
                    value, lo, hi = 1-value, 1-hi, 1-lo
                row.update({name: value, name+"_ci_low": lo, name+"_ci_high": hi,
                            name+"_n_samples": stat["n_samples"] if stat else 0,
                            name+"_n_independent_units": stat['n_independent_units'] if stat else 0})
            if target in targets: comparisons.append(row)
            by_sample = defaultdict(list)
            for pair in target_pairs:
                by_sample[pair['sample_id']].append(pair)
            for sample, paired in sorted(by_sample.items()):
                valid = [p for p in paired if p["evaluable"] and p["is_seed_region"]]
                signal = [p for p in valid if p["signal_positive"]]
                state = "evaluable" if valid and len(valid) == len(paired) else "partial" if valid else "not_evaluable"
                item = dict(sample_id=sample, reference_version=ref, source_version=source, target_version=target,
                            analysis_kind=kind, comparison_scope=scope, identity_anchor=target == source,
                            status=state, n_units=sum(p["is_seed_region"] for p in paired), n_evaluable_units=len(valid), n_signal_units=len(signal),
                            target_statuses=sorted({p["status"] for p in paired}), source_statuses=sorted({p["comparator_status"] for p in paired}),
                            interpretation="observed_similarity_only; no_single_sample_equivalence")
                item.update(_label_scope(valid))
                for name, _ in SIMILARITIES:
                    metric = "label_discordance" if name == "label_agreement" else name
                    values = [p[metric] for p in signal if p[metric] is not None]
                    value = _mean(values)
                    item[name] = 1-value if name == "label_agreement" and value is not None else value
                for metric in ("extra_detection", "missed_detection", "class_change"):
                    item[metric] = sum(p[metric] for p in valid if p[metric] is not None) if valid else None
                sample_rows.append(item)
    return comparisons, sample_rows, links


def _upstream_evaluable(row):
    return (row.get('status') in {'retained', 'not_retained'}
            and row.get('comparison_eligible') is True
            and bool(row.get('matching_rule')) and bool(row.get('cnvkit_version')))


def _upstream_context(row):
    seed = row.get('seed')
    if not isinstance(seed, dict):
        return None
    chrom, start, end = (seed.get(k) for k in ('chrom', 'start', 'end'))
    threshold = row.get('min_seed_overlap')
    mode, rule = row.get('reference_mode'), row.get('matching_rule')
    if row.get('global_scope') and rule == LOCUS_RULE:
        if isinstance(chrom, str) and type(start) is int and type(end) is int and 0 <= start < end and mode:
            return chrom, start, end, mode, rule
        return None
    if (not isinstance(chrom, str) or not chrom or type(start) is not int or type(end) is not int
            or not 0 <= start < end or mode not in {'matched_normal', 'flat'}
            or not isinstance(rule, str) or not rule or type(threshold) not in (int, float)
            or not 0 < threshold <= 1):
        return None
    return chrom, start, end, mode, rule, threshold


def _control_status(row):
    # Earlier reports do not identify the calibration role. Do not promote
    # those rows into independently assessed source controls.
    if row.get('calibration_role') != 'source_control' or not _upstream_evaluable(row):
        return 'unverified'
    status = row.get('calibration_status', 'unverified')
    if not row.get('global_scope') and status == 'passed' and row['status'] != 'retained':
        return 'failed'
    return status


def _upstream(records, samples, targets, ref, alpha, decision_alpha, minimum, margin,
              patient_ids=None, primary_pairs=None, global_mode=False, seed_match_fraction=1.0):
    lookup, observed = {}, defaultdict(set)
    allowed = {"retained", "not_retained", "failed", "missing", "unsupported", "no_seeds", "invalid", "not_run"}
    for original in records:
        row = dict(original)
        sample, version, seed = row.get("sample_id"), row.get("version"), row.get("seed_id")
        if sample not in samples or not version or not seed or row.get("status") not in allowed:
            raise ValueError("Malformed upstream sample/version/seed/status")
        if row.get("source_version") != samples[sample]["source_version"]:
            raise ValueError("Upstream source_version differs from declared input source")
        if (sample, version, seed) in lookup:
            raise ValueError("Duplicate upstream sample/version/seed")
        fraction = row.get("retention_fraction")
        if fraction is not None and (not math.isfinite(float(fraction)) or not 0 <= float(fraction) <= 1):
            raise ValueError("retention_fraction must be [0,1] or null")
        if row.get("calibration_status", "unverified") not in {"passed", "unverified", "failed"}:
            raise ValueError("Unknown calibration_status")
        role = row.get('calibration_role')
        expected_role = 'source_control' if version == row['source_version'] else 'candidate_observation'
        if role is not None and role != expected_role:
            raise ValueError('Upstream calibration_role contradicts the source and target versions')
        lookup[sample, version, seed] = row
        if row["status"] != "no_seeds" and not (global_mode and seed == "__no_seeds__"):
            observed[sample].add(seed)
    universe = {sample: samples[sample].get("upstream_seed_ids", sorted(observed[sample])) for sample in samples}
    for sample in samples:
        if "upstream_seed_ids" in samples[sample] and observed[sample] - set(universe[sample]):
            raise ValueError("Upstream records contain seeds outside the declared candidate universe")
    seed_rows, sample_rows, statistics = [], [], []
    sources = sorted({s["source_version"] for s in samples.values()})
    primary_keys = _primary_pair_keys(primary_pairs)
    family_count = sum(_is_primary_pair(t, source, ref, primary_keys) for source in sources for t in targets)
    adjusted_alpha = decision_alpha / family_count if family_count else None
    unit = 'patient' if patient_ids is not None else 'sample'
    for source in sources:
        selected = sorted(s for s in samples if samples[s]["source_version"] == source)
        # Always preserve source controls even if the source is not a target.
        for target in dict.fromkeys([*targets, source]):
            local = []
            primary = _is_primary_pair(target, source, ref, primary_keys)
            for sample in selected:
                ids = universe[sample]
                rows = []
                for seed in ids or ["__no_seeds__" if "upstream_seed_ids" in samples[sample] else "__missing_seed_manifest__"]:
                    original = lookup.get((sample, target, seed))
                    control = lookup.get((sample, source, seed), {})
                    status = "no_seeds" if not ids and "upstream_seed_ids" in samples[sample] else "missing"
                    row = dict(original or dict(sample_id=sample, version=target, seed_id=seed, source_version=source,
                               status=status, reason="No declared input candidates" if status == "no_seeds" else "No upstream result for expected candidate",
                               comparison_eligible=False, calibration_status="unverified", retention_fraction=None))
                    row.update(reference_version=ref, target_version=target, candidate_universe_declared="upstream_seed_ids" in samples[sample])
                    row["evaluable"] = _upstream_evaluable(row)
                    row["seed_loss"] = float(row["status"] == "not_retained") if row["evaluable"] else None
                    if row["status"] in {"retained", "not_retained"} and not row["evaluable"]:
                        row["evaluation_reason"] = "Missing verified comparison eligibility, matching rule or actual CNVkit version"
                    context = _upstream_context(row)
                    paired = (target != source and row['evaluable'] and _upstream_evaluable(control)
                              and context is not None and context == _upstream_context(control))
                    row.update(source_control_status=_control_status(control),
                               source_control_recipe_match=control.get('native_recipe_match') is True,
                               source_control_cnvkit_version=control.get('cnvkit_version'),
                               retention_difference=(float(row['status'] == 'retained') - float(control['status'] == 'retained')
                                                     if paired else None))
                    if global_mode:
                        row['source_locus_present'] = control['status'] == 'retained' if _upstream_evaluable(control) else None
                        row['locus_discordance'] = None
                        if row['evaluable'] and _upstream_evaluable(control) and context and context == _upstream_context(control):
                            own, other = row['status'] == 'retained', row['source_locus_present']
                            if own or other:
                                row['locus_discordance'] = float(own != other)
                    rows.append(row)
                    seed_rows.append(row)
                valid = [r for r in rows if r["evaluable"]]
                fractions = [float(r["retention_fraction"]) for r in valid if r.get("retention_fraction") is not None]
                deltas = [r['retention_difference'] for r in rows if r['retention_difference'] is not None]
                control_statuses = sorted({r['source_control_status'] for r in rows})
                item = dict(sample_id=sample, reference_version=ref, source_version=source, target_version=target,
                            expected_seeds=len(ids), evaluable_seeds=len(valid), unavailable_seeds=len(ids)-len(valid),
                            no_seeds=not ids and "upstream_seed_ids" in samples[sample],
                            status_counts=dict(Counter(r["status"] for r in rows)),
                            seed_loss=_mean([r["seed_loss"] for r in valid]),
                            mean_retention_fraction=_mean(fractions),
                            retention_difference=_mean(deltas), paired_evaluable_seeds=len(deltas),
                            paired_candidate_observations=sum(r['retention_difference'] is not None
                                and r.get('calibration_role') == 'candidate_observation' for r in rows),
                            candidate_universe_declared="upstream_seed_ids" in samples[sample],
                            calibration_passed=bool(ids) and control_statuses == ['passed'],
                            calibration_statuses=control_statuses,
                            source_control_recipe_match=bool(ids) and all(r['source_control_recipe_match'] for r in rows),
                            source_control_cnvkit_versions=sorted({r['source_control_cnvkit_version'] for r in rows
                                                                  if r['source_control_cnvkit_version']}),
                            native_recipe_match=bool(valid) and len(valid) == len(ids) and all(r.get("native_recipe_match") is True for r in valid),
                            native_recipe_statuses=sorted({"matched" if r.get("native_recipe_match") is True else "controlled" if r.get("native_recipe_match") is False else "unverified" for r in rows}),
                            cnvkit_versions=sorted({str(r["cnvkit_version"]) for r in valid}),
                            interpretation="descriptive_known_candidate_retention; not_genome_wide_sensitivity")
                if global_mode:
                    from .concordance_upstream_similarity import native_seed_comparison
                    candidate_rows = {r['seed_id']: r for r in rows}
                    source_rows = {seed: lookup.get((sample, source, seed), {}) for seed in candidate_rows}
                    comparable = bool(ids) and all(
                        r['evaluable'] and _upstream_evaluable(source_rows[key])
                        and _upstream_context(r) is not None and _upstream_context(r) == _upstream_context(source_rows[key])
                        for key, r in candidate_rows.items())
                    item.update(native_seed_comparison(candidate_rows, source_rows, seed_match_fraction, comparable))
                    informative = [r for r in rows if r.get('locus_discordance') is not None]
                    item.update(global_scope=True,
                        no_seeds=(samples[sample].get('global_seed_sets', {}).get(target, {}).get('status') == 'ok'
                                  and not samples[sample]['global_seed_sets'][target]['seeds']),
                        interpretation='version_native_seed_presence_at_global_panel_union_loci',
                        native_seed_count=(len(samples[sample]['global_seed_sets'][target]['seeds'])
                            if samples[sample].get('global_seed_sets', {}).get(target, {}).get('status') == 'ok' else None),
                        informative_pair_loci=len(informative),
                        locus_discordance=_mean([r['locus_discordance'] for r in informative]),
                        extra_loci=sum(1 for r in informative if r['locus_discordance'] and not r['source_locus_present']),
                        missed_loci=sum(1 for r in informative if r['locus_discordance'] and r['source_locus_present']))
                local.append(item)
                sample_rows.append(item)
            sample_losses = {r['sample_id']: r['seed_loss'] for r in local if r['seed_loss'] is not None}
            if global_mode:
                sample_losses = {r['sample_id']: r[GLOBAL_METRIC] for r in local if r[GLOBAL_METRIC] is not None}
            values = _independent_means(sample_losses, patient_ids).tolist()
            sample_deltas = {r['sample_id']: r['retention_difference'] for r in local if r['retention_difference'] is not None}
            differences = _independent_means(sample_deltas, patient_ids).tolist()
            low, high = _bounds(values, alpha)
            delta_low, delta_high = _bounds(differences, alpha, limits=(-1., 1.))
            decision_low, decision_high = _bounds(values, adjusted_alpha, two_sided=False) if primary else (None, None)
            reasons = []
            if margin is None: reasons.append("upstream_loss_margin_not_predefined")
            if len(values) < minimum: reasons.append("too_few_independent_units")
            if any(r["unavailable_seeds"] or r["seed_loss"] is None for r in local): reasons.append("unavailable_or_no_candidate_samples")
            if any(not r["candidate_universe_declared"] for r in local): reasons.append("candidate_universe_not_declared")
            if global_mode and any(r[GLOBAL_METRIC] is None for r in local):
                reasons.append('missing_paired_native_seed_sets_or_both_empty')
            candidate_versions = {v for r in local for v in r['cnvkit_versions']}
            control_versions = {v for r in local for v in r['source_control_cnvkit_versions']}
            if any(len(r["cnvkit_versions"]) > 1 for r in local): reasons.append("multiple_cnvkit_versions_within_sample_target")
            if len(candidate_versions) > 1: reasons.append("cnvkit_version_differs_across_samples")
            decision = 'insufficient_evidence'
            if source == target:
                reasons = ['source_control_not_version_comparison']
            elif not primary:
                reasons = ['exploratory_pair_not_primary']
            elif not reasons and decision_high < float(margin):
                decision = 'within_margin'
            elif margin is not None and len(values) >= minimum and decision_low > float(margin):
                decision = 'exceeds_margin'
            else:
                reasons.append('confidence_bounds_do_not_establish_agreement')
            control_passed = all(r['calibration_passed'] and r['candidate_universe_declared'] for r in local)
            control_recipe = all(r['source_control_recipe_match'] for r in local) and len(control_versions) == 1
            native = all(r['native_recipe_match'] for r in local) and len(candidate_versions) == 1
            contexts_matched = (source != target and all(r['candidate_universe_declared'] and r['expected_seeds'] > 0
                                and r['paired_candidate_observations'] == r['expected_seeds'] for r in local))
            p_bound = None
            if primary and values and margin is not None:
                p_bound = math.exp(-2*len(values)*(float(margin)-_mean(values))**2) if _mean(values) < float(margin) else 1.
            statistics.append(dict(reference_version=ref, source_version=source, target_version=target,
                metric=GLOBAL_METRIC if global_mode else "known_candidate_not_retained_fraction", estimate=_mean(values),
                conservative_ci_low=low, conservative_ci_high=high, confidence_level=1-alpha,
                agreement_lower=decision_low, agreement_upper=decision_high,
                agreement_alpha=adjusted_alpha if primary else None, primary_pair=primary,
                n_samples=len(sample_losses), n_independent_units=len(values), independence_unit=unit,
                selected_samples=len(selected), expected_seeds=sum(r["expected_seeds"] for r in local),
                evaluable_seeds=sum(r["evaluable_seeds"] for r in local), no_seed_samples=sum(r["no_seeds"] for r in local),
                retention_difference=_mean(differences), retention_difference_ci_low=delta_low, retention_difference_ci_high=delta_high,
                retention_difference_n_samples=len(sample_deltas), retention_difference_n_independent_units=len(differences),
                retention_difference_n_seeds=sum(r['paired_evaluable_seeds'] for r in local),
                upstream_loss_margin=margin, min_samples=minimum, decision=decision, reasons=sorted(set(reasons)),
                calibration_passed=control_passed,
                observed_statuses=sorted({status for r in local for status in r['status_counts']}),
                calibration_statuses=sorted({c for r in local for c in r["calibration_statuses"]}),
                calibration_basis='original_source_full_seed_bed_reproduction' if global_mode else 'independent_same_version_source_control',
                source_control_recipe_match=control_recipe,
                source_control_cnvkit_versions=sorted(control_versions),
                paired_context_complete=contexts_matched,
                native_recipe_match=native,
                native_recipe_statuses=sorted({v for r in local for v in r["native_recipe_statuses"]}),
                p_agreement_bound=p_bound,
                comparison_scope=GLOBAL_SCOPE if global_mode else f"known_candidate_seed_retention_only; {unit}_equal_weight; no_genome_wide_detection_validation"))
            if global_mode:
                statistics[-1].update(seed_match_fraction=seed_match_fraction, informative_pair_loci=sum(r['informative_pair_loci'] for r in local),
                                      extra_loci=sum(r['extra_loci'] for r in local),
                                      missed_loci=sum(r['missed_loci'] for r in local))
    return seed_rows, sample_rows, statistics


def _combine(downstream, upstream, mode, samples, global_mode=False):
    lookup = {(r["source_version"], r["target_version"]): r for r in upstream}
    result = []
    for down in downstream:
        up = lookup.get((down["source_version"], down["target_version"]))
        includes_upstream = mode == "with_bam" and down["analysis_kind"] == "aa_ac"
        decision, reasons = down["decision"], list(down["reasons"])
        basis = ["downstream_metrics_exceed_margin"] if decision == "exceeds_margin" else ["downstream_metrics_within_margin"] if decision == "within_margin" else []
        scope = "conditional_local_upstream_and_fixed_seed_aa_ac" if includes_upstream else "ac_only_on_fixed_source_aa" if down["analysis_kind"] == "ac_only" else "provided_downstream_aa_ac_results_only"
        if includes_upstream:
            if global_mode:
                scope = GLOBAL_SCOPE
            upstream_validated = bool(up and up['calibration_passed'] and up['source_control_recipe_match']
                                      and up['native_recipe_match'] and up['paired_context_complete'])
            basis = ["downstream_metrics_exceed_margin"] if down["decision"] == "exceeds_margin" else []
            if up and up["decision"] == "exceeds_margin" and upstream_validated:
                basis.append("global_upstream_exceeds_margin_after_source_control" if global_mode else "local_upstream_exceeds_margin_after_source_control")
            if basis:
                decision = "exceeds_margin"
            elif down["decision"] == "within_margin" and up and up["decision"] == "within_margin" and upstream_validated:
                decision = "within_margin"
                basis = ["downstream_and_global_upstream_within_margin_after_source_control" if global_mode else "downstream_and_local_upstream_within_margin_after_source_control"]
            else:
                decision = "insufficient_evidence"
            reasons += up["reasons"] if up else ["upstream_results_missing"]
            if up and up["decision"] == "exceeds_margin" and not upstream_validated:
                reasons.append("upstream_exceedance_not_validated_by_independent_source_control")
            if not up or not up["calibration_passed"]: reasons.append("original_global_source_seed_set_not_reproduced" if global_mode else "local_upstream_calibration_not_passed")
            if not up or not up['source_control_recipe_match']: reasons.append('source_control_recipe_not_verified')
            if not up or not up['paired_context_complete']: reasons.append('upstream_paired_candidate_scope_incomplete_or_different')
            if not up or not up["native_recipe_match"]:
                reasons.append("controlled_CNVkit_not_historical_recipe" if up and "controlled" in up["native_recipe_statuses"] else "native_CNVkit_recipe_not_verified")
        bridge_statuses = sorted({s.get("source_bridge_status", "unverified") for s in samples.values()
                                  if s["source_version"] == down["source_version"]})
        bridge_passed = bool(bridge_statuses) and bridge_statuses == ["passed"]
        if mode == "with_bam" and not bridge_passed:
            decision, basis = "insufficient_evidence", []
            reasons.append("input_source_to_rerun_bridge_not_passed")
        limited = global_mode and any(
            s.get('global_seed_sets', {}).get(version, {}).get('historical_provenance_limited', False)
            for s in samples.values() if s['source_version'] == down['source_version']
            for version in (down['source_version'], down['target_version']))
        if limited:
            decision, basis = "insufficient_evidence", []
            reasons.append("historical_global_cnv_computation_provenance_limited")
        if not down['primary_pair'] and not down['identity_anchor']:
            decision, basis = 'insufficient_evidence', []
        result.append(dict(reference_version=down["reference_version"], source_version=down["source_version"],
            target_version=down["target_version"], analysis_kind=down["analysis_kind"], mode=mode,
            comparison_scope=scope, decision=decision, decision_basis=basis, reasons=sorted(set(reasons)),
            downstream_decision=down["decision"], upstream_decision=up["decision"] if includes_upstream and up else "not_assessed",
            upstream_calibration_passed=up["calibration_passed"] if includes_upstream and up else None,
            upstream_calibration_basis=up['calibration_basis'] if includes_upstream and up else None,
            upstream_source_control_statuses=up['calibration_statuses'] if includes_upstream and up else None,
            upstream_source_control_recipe_match=up['source_control_recipe_match'] if includes_upstream and up else None,
            upstream_paired_context_complete=up['paired_context_complete'] if includes_upstream and up else None,
            native_recipe_match=up["native_recipe_match"] if includes_upstream and up else None,
            source_bridge_passed=bridge_passed if mode == "with_bam" else None,
            source_bridge_statuses=bridge_statuses if mode == "with_bam" else [],
            n_selected_samples=down["n_selected_samples"], signal_samples=down["signal_samples"],
            signal_independent_units=down['signal_independent_units'], independence_unit=down['independence_unit'],
            primary_pair=down['primary_pair'],
            compatibility_metric_scope=down["compatibility_metric_scope"],
            not_established=["biological_seed_discovery_accuracy" if global_mode else "genome_wide_seed_discovery", "complete_suite_interchangeability", "biological_truth", "cycle_order_or_breakpoint_topology"]
                + (["aa_reconstruction_compatibility"] if down["analysis_kind"] == "ac_only" else [])))
    return result


def _sample_ladders(out, sample_rows, ref, mode, synthetic, upstream_samples=None):
    import matplotlib.pyplot as plt
    groups = defaultdict(list)
    for row in sample_rows: groups[row["sample_id"], row["analysis_kind"], row["source_version"]].append(row)
    folder = out / "sample_ladders"; folder.mkdir(exist_ok=True)
    indices = []
    upstream_lookup = {(r['sample_id'], r['source_version'], r['target_version']): r for r in (upstream_samples or [])}
    codes = {"failed": "F", "missing": "M", "unsupported": "U", "invalid": "I", "interrupted": "X",
             "blocked": "B", "no_seeds": "NS", "not_run": "NR"}
    for (sample, kind, source), rows in groups.items():
        has_upstream = mode == 'with_bam' and kind == 'aa_ac'
        global_mode = any(r.get('global_scope') for r in (upstream_samples or []))
        panels = (6 if global_mode else 7) if has_upstream else 4
        fig, axes = plt.subplots(1, panels,
                                 figsize=(23 if has_upstream else 14, max(4.6, .48*len(rows)+2.6)), sharey=True)
        labels = [r["target_version"] + (" [source]" if r["target_version"] == source else "")
                  + (" [ref]" if r["target_version"] == ref else "")
                  + f"   {r['n_evaluable_units']}/{r['n_units']} usable"
                  for r in rows]
        for ax, (metric, label) in zip(axes, SIMILARITIES):
            for y, row in enumerate(rows):
                value, status = row.get(metric), ""
                if row["status"] == "partial": value, status = None, "PART"
                elif row["status"] == "not_evaluable":
                    target_codes = "/".join(codes[s] for s in row["target_statuses"] if s in codes)
                    source_codes = "/".join(codes[s] for s in row["source_statuses"] if s in codes)
                    value, status = None, (target_codes + (" SRC:"+source_codes if source_codes else "")).strip() or "NA"
                elif not row["n_signal_units"]: value, status = None, "ND/ND"
                elif metric in {"label_agreement", "class_footprint_jaccard"} and not row["full_label_comparable"]:
                    value, status = None, "CLASS*"
                if value is None:
                    ax.axhspan(y-.33, y+.33, color="#eeeeee", zorder=0)
                    ax.text(.5, y, status or "NA", ha="center", va="center", color="#666666", fontsize=7)
                else:
                    identity = row["identity_anchor"]
                    ax.scatter(value, y, marker="D" if identity else "o", color="#333333" if identity else "#4178A9", s=25)
                    ax.text(value-.025 if value > .75 else value+.025, y-.16,
                            "identity" if identity else f"{value*100:.1f}%", ha="right" if value > .75 else "left", fontsize=7)
            ax.set_xlim(-.03, 1.03); ax.set_ylim(len(rows)-.5, -.5)
            ax.set_xticks([0, .5, 1], ["0", "50", "100"]); ax.set_xlabel(label+" (%)")
            ax.grid(axis="x", color="#eeeeee"); ax.set_axisbelow(True)
        # Counts retain the whole evaluable candidate scope, including regions
        # without signal. The four plotted similarities use signal-positive units.
        for y, row in enumerate(rows):
            axes[0].text(.02, y+.25,
                f"signal={row['n_signal_units']}  +{row.get('extra_detection', 0)} / -{row.get('missed_detection', 0)}  class={row.get('class_change', 0)}",
                fontsize=6, color='#666666')
        if has_upstream:
            sim_ax, ax = axes[-2:]
            source_ax = None if global_mode else axes[-3]
            baseline = None
            for y, row in enumerate(rows):
                up = upstream_lookup.get((sample, source, row['target_version']))
                matched = (up or {}).get('source_matched_cnv') or {}
                baseline = matched.get('source_cnvkit_version') or baseline
                matched_value = matched.get('interval_jaccard')
                if source_ax is not None:
                    if matched_value is None:
                        source_ax.axhspan(y-.33, y+.33, color='#eeeeee', zorder=0)
                        state = 'both empty' if matched.get('status') == 'both_empty' else 'NA'
                        source_ax.text(.5, y, state, ha='center', va='center', fontsize=7)
                    else:
                        source_ax.scatter(matched_value, y, marker='D' if matched.get('identity_anchor') else 'o',
                                          s=25, color='#315e77')
                        source_ax.text(.5, y+.24, f"CNVkit {matched['target_cnvkit_version']} · {matched_value:.1%}",
                                       ha='center', fontsize=6)
                sim = up.get('upstream_interval_jaccard') if up else None
                if sim is None:
                    sim_ax.axhspan(y-.33, y+.33, color='#eeeeee', zorder=0)
                    empty = up and up.get('upstream_similarity_status') == 'both_empty'
                    sim_ax.text(.5, y, 'both empty' if empty else 'NA', ha='center', va='center', fontsize=7)
                else:
                    sim_ax.scatter(sim, y, s=25, color='#4178A9')
                    sim_ax.text(sim-.025 if sim > .75 else sim+.025, y-.16, f'{100*sim:.1f}%',
                                ha='right' if sim > .75 else 'left', fontsize=7)
                value = None if not up or up.get('seed_loss') is None else 1-up['seed_loss']
                partial = bool(up and up['evaluable_seeds'] != up['expected_seeds'])
                if value is None or partial:
                    ax.axhspan(y-.33, y+.33, color='#eeeeee', zorder=0)
                    states = '/'.join(codes.get(s, s) for s in (up.get('status_counts', {}) if up else {}))
                    ax.text(.5, y, ('PART ' if partial else '') + (states or 'NA'), ha='center', va='center', fontsize=7, color='#666666')
                else:
                    ax.scatter(value, y, s=25, color='#4178A9')
                    ax.text(value-.025 if value > .75 else value+.025, y-.16, f'{100*value:.1f}%',
                            ha='right' if value > .75 else 'left', fontsize=7)
                if up:
                    flags = []
                    if not up.get('calibration_passed'): flags.append('SRC?')
                    if not up.get('native_recipe_match'): flags.append('NATIVE?')
                    ax.text(.02, y+.25, f"{up.get('upstream_retained_count', 'NA')}/{up['evaluable_seeds']} retained; expected {up['expected_seeds']} " + ' '.join(flags), fontsize=6, color='#666666')
            ax.set_xlim(-.03, 1.03); ax.set_ylim(len(rows)-.5, -.5)
            ax.set_xticks([0, .5, 1], ['0', '50', '100'])
            sim_ax.set_xlim(-.03, 1.03); sim_ax.set_xticks([0, .5, 1], ['0', '50', '100'])
            if source_ax is not None:
                source_ax.set_xlim(-.03, 1.03); source_ax.set_xticks([0, .5, 1], ['0', '50', '100'])
                source_ax.set_xlabel(f'Local CNVkit vs {baseline or "unknown"}\nseed Jaccard (%)')
                source_ax.grid(axis='x', color='#eeeeee')
            sim_ax.set_xlabel('Suite recipe seed Jaccard (%)'); sim_ax.grid(axis='x', color='#eeeeee')
            ax.set_xlabel('Union loci with native seeds (%)' if global_mode else 'Known candidates retained (%)'); ax.grid(axis='x', color='#eeeeee'); ax.set_axisbelow(True)
        axes[0].set_yticks(range(len(rows)), labels)
        fig.suptitle(f'Sample: {sample} | Input-source comparator: {source} | Fixed reference: {ref}\n'
                     + _kind_label(kind, mode) + " | observed signal-positive input units"
                     + (' | SYNTHETIC TEST DATA' if synthetic else ''), fontsize=10)
        fig.text(.5, .015, "One sample: observed similarity only; no confidence interval, compatibility probability or equivalence decision. Source diamond = identity.\nGray: F failed, M missing, U unsupported, I invalid, X interrupted, B blocked, NS no seeds, NR not run, PART incomplete, ND/ND both negative, NA unavailable; SRC = input source.\nCounts: usable/expected units; signal units; +extra/-missed detections; class changes. CLASS*: shared classes only. SRC?/NATIVE?: upstream control/recipe unverified.", ha="center", fontsize=6.5)
        fig.tight_layout(rect=[0, .14, 1, 1])
        index = dict(sample_id=sample, analysis_kind=kind, source_version=source, reference_version=ref,
                     comparison_scope=row["comparison_scope"], mode=mode, synthetic=synthetic)
        for ext in ("png", "pdf", "svg"):
            path = folder / (_slug(sample+"__"+kind+"__"+source) + "." + ext)
            fig.savefig(path, dpi=170, bbox_inches="tight"); index[ext] = str(path.relative_to(out))
        indices.append(index); plt.close(fig)
    return indices


def _upstream_ladder(out, upstream, ref, mode, synthetic):
    import matplotlib.pyplot as plt
    sources = list(dict.fromkeys(r["source_version"] for r in upstream))
    if not sources or mode != "with_bam": return
    longest = max(sum(r["source_version"] == s for r in upstream) for s in sources)
    fig, axes = plt.subplots(len(sources), 1, figsize=(9, max(4., len(sources)*(.37*longest+1.4)+1.2)), squeeze=False)
    global_mode = any(r.get("metric") == GLOBAL_METRIC for r in upstream)
    for ax, source in zip(axes.flat, sources):
        rows = [r for r in upstream if r["source_version"] == source]
        for y, row in enumerate(rows):
            if row["estimate"] is None:
                ax.axhspan(y-.33, y+.33, color="#eeeeee", zorder=0)
                ax.text(.5, y, "not evaluable", ha="center", va="center", color="#777777", fontsize=8)
                continue
            value, low, high = 1-row["estimate"], 1-row["conservative_ci_high"], 1-row["conservative_ci_low"]
            checked = (row['calibration_passed'] and row['source_control_recipe_match']
                       and row['native_recipe_match'] and (row['paired_context_complete'] or row['target_version'] == source))
            flags = []
            if not row["calibration_passed"]: flags.append("SRC!" if "failed" in row["calibration_statuses"] else "SRC?")
            if not row['source_control_recipe_match']: flags.append('SRC NATIVE?')
            if row['target_version'] != source and not row['paired_context_complete']: flags.append('PAIR?')
            if not row["native_recipe_match"]: flags.append("CTRL" if "controlled" in row["native_recipe_statuses"] else "NATIVE?")
            ax.errorbar(value, y, xerr=[[max(0, value-low)], [max(0, high-value)]], fmt="o" if checked else "^",
                        color="#4178A9" if checked else "#C5884E", markersize=5, capsize=3)
            ax.text(.02, y-.19, f'n={row["n_independent_units"]}' + ('; ' + ", ".join(flags) if flags else ''), fontsize=7, color="#666666")
        ax.set_yticks(range(len(rows)), [r["target_version"] + (" [source re-evaluation]" if r["target_version"] == source else "")
                                       + (" [ref]" if r["target_version"] == ref else "") for r in rows])
        ax.set_ylim(len(rows)-.5, -.5); ax.set_xlim(-.03, 1.03)
        ax.set_xticks([0, .25, .5, .75, 1], ["0", "25", "50", "75", "100"])
        ax.set_xlabel("Lower directional native-seed retention (%)" if global_mode else "Known candidates retained (%)"); ax.grid(axis="x", color="#eeeeee"); ax.set_axisbelow(True)
        ax.set_title(f"Source version: {source}", fontsize=9, loc="left")
    confidence = upstream[0]["confidence_level"]*100
    fig.suptitle((f'Global native-seed retention | Fixed reference: {ref}' if global_mode else f'Known-candidate retention | Fixed reference: {ref}\nBounded local upstream; not genome-wide detection sensitivity')
                 + (' | SYNTHETIC TEST DATA' if synthetic else ''), fontsize=11)
    unit = upstream[0]['independence_unit']
    fig.text(.5, .015, f'Equal-weight {unit} means and conservative {confidence:g}% bounds. Source control is an observed result.\nSRC?/SRC! unverified/failed source control; PAIR? incomplete matched candidates; CTRL/NATIVE? controlled/unverified recipe.', ha="center", fontsize=6.5)
    fig.tight_layout(rect=[0, .09, 1, 1])
    for ext in ("png", "pdf", "svg"): fig.savefig(out / f"upstream_retention_ladder.{ext}", dpi=170, bbox_inches="tight")
    plt.close(fig)


def _local_cnv_cohort_plot(out, rows, ref, synthetic):
    import matplotlib.pyplot as plt
    if not rows:
        return
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(11, max(3.3, 1.9 + .45 * len(rows))))
    labels = [f"{r['source_version']}, {r['reference_build'] or '?'} / {r['reference_mode']}  |  "
              f"{r['source_cnvkit_version'] or '?'} → {r['target_cnvkit_version'] or '?'}" for r in rows]
    for ax, metric, title, low_limit in zip(axes, ('interval_jaccard', 'retention_delta'),
            ('Final-seed interval Jaccard (%)', 'Candidate retention difference (percentage points)'), (0, -100)):
        ax.set_xlim(low_limit - 5, 105)
        ax.set_xticks([0, 50, 100] if low_limit == 0 else [-100, -50, 0, 50, 100])
        ax.set_xlabel(title)
        ax.grid(axis='x', color='#eee')
        ax.set_axisbelow(True)
        for index, row in enumerate(rows):
            value = row[metric]
            if value is None:
                ax.text((low_limit + 100) / 2, index, 'not evaluable', va='center', ha='center', color='#777777')
                continue
            low, high = row[metric + '_ci_low'], row[metric + '_ci_high']
            error = [[100*max(0, value-low)], [100*max(0, high-value)]] if low is not None and high is not None else None
            identity = row['identity_anchor']
            ax.errorbar(100*value, index, xerr=error, fmt='D' if identity else 'o',
                        markersize=4, capsize=3, color='#444444' if identity else '#4178A9')
            ax.text(low_limit, index-.18, 'identity' if identity else f"n={row[metric+'_n_independent_units']}",
                    fontsize=6, color='#555555')
    axes[0].set_yticks(range(len(rows)), labels)
    axes[0].set_ylim(len(rows)-.5, -.6)
    conf = 100*rows[0]['confidence_level']
    fig.suptitle(f"Local CNVkit controls: original CNVkit → target | Suite panel reference {ref}"+
                 (" | SYNTHETIC" if synthetic else ""), fontsize=10)
    fig.text(.5, .02, f"Same source AA filter within each contrast. {conf:g}% pointwise descriptive intervals. "
             "n = independent patients or samples.", ha='center', fontsize=7)
    fig.tight_layout(rect=[0, .08, 1, 1])
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(out / f'local_cnv_comparator_cohort.{ext}', dpi=170, bbox_inches='tight')
    plt.close(fig)


def _plots(out, downstream, combined, ref, mode, synthetic, sample_rows=None, upstream=None, upstream_samples=None, matched_cohort=None):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "ascompare-matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import numpy as np
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
    groups = list(dict.fromkeys((r["analysis_kind"], r["source_version"]) for r in downstream))
    targets = list(dict.fromkeys(r["target_version"] for r in downstream))
    height = max(4., (len(targets)*.34+1.2)*len(groups)+1.2)
    fig, axes = plt.subplots(len(groups), 4, figsize=(13, height), squeeze=False)
    for index, (kind, source) in enumerate(groups):
        rows = [r for r in downstream if r["analysis_kind"] == kind and r["source_version"] == source]
        for ax, (metric, title) in zip(axes[index], SIMILARITIES):
            for y, row in enumerate(rows):
                value = row.get(metric)
                if row["identity_anchor"]:
                    ax.text(.5, y, "input identity; not validation", color="#777777", ha="center", va="center", fontsize=6)
                elif value is None:
                    ax.text(.5, y, "not evaluable", color="#777777", ha="center", va="center", fontsize=7)
                elif metric in {"label_agreement", "class_footprint_jaccard"} and not row["full_label_comparable"]:
                    ax.axhspan(y-.33, y+.33, color="#eeeeee", zorder=0)
                    ax.text(.5, y, "shared classes only*", color="#666666", ha="center", va="center", fontsize=7)
                else:
                    low, high = row[metric+"_ci_low"], row[metric+"_ci_high"]
                    ax.errorbar(value, y, xerr=[[max(0, value-low)], [max(0, high-value)]], fmt="o", markersize=4,
                                color="#4178A9", capsize=3)
                    ax.text(.02, y-.18, f'n={row[metric+"_n_independent_units"]}', fontsize=6, color="#666666")
            ax.set_xlim(-.03, 1.03); ax.set_ylim(len(rows)-.5, -.5)
            ax.set_xticks([0, .5, 1], ["0", "50", "100"]); ax.set_xlabel(title+" (%)")
            ax.set_yticks(range(len(rows)), [r["target_version"] for r in rows] if ax is axes[index, 0] else [])
            ax.grid(axis="x", color="#eeeeee"); ax.set_axisbelow(True)
        axes[index, 0].set_title(f'{_kind_label(kind, mode)}\nInput-source version comparator: {source}', loc="left", fontsize=9)
    fig.suptitle(f'AmpliconSuite Cross-Version Concordance | Fixed reference: {ref}\nSource-to-target, signal-positive selected units | {mode}'
                 + (' | SYNTHETIC TEST DATA' if synthetic else ''), fontsize=11)
    confidence = 100*downstream[0]["confidence_level"] if downstream else 95
    unit = downstream[0]['independence_unit'] if downstream else 'sample'
    fig.text(.5, .01, f"Equal-weight {unit} means and conservative {confidence:g}% bounds. Signal-positive denominators can differ between pairs.\n* Class support differs or is unknown: gray label/class-footprint panels. Shared-class values are in the tables.\nAC-only panels assess classification on fixed AA output. Comparisons cover selected input units.", ha="center", fontsize=7)
    fig.tight_layout(rect=[0, .11, 1, 1])
    for extension in ("png", "pdf", "svg"): fig.savefig(out / f"similarity_ladder.{extension}", dpi=170, bbox_inches="tight")
    plt.close(fig)
    lookup = {(r["analysis_kind"], r["source_version"], r["target_version"]): r for r in combined}
    values = np.array([[{"insufficient_evidence": 0, "within_margin": 1, "exceeds_margin": 2}[lookup[k, s, t]["decision"]]
                        for t in targets] for k, s in groups])
    fig, ax = plt.subplots(figsize=(max(7, 1.45*len(targets)), max(3.4, .8*len(groups)+2.3)))
    ax.imshow(values, cmap=ListedColormap(["#e5e5e5", "#c8d9e8", "#ebc6a6"]), vmin=-.5, vmax=2.5, aspect="auto")
    for i, (kind, source) in enumerate(groups):
        for j, target in enumerate(targets):
            label = DECISIONS[lookup[kind, source, target]["decision"]]
            ax.text(j, i, label.replace(" ", "\n", 1), ha="center", va="center", fontsize=8)
    ax.set_xticks(range(len(targets)), targets, rotation=30, ha="right")
    ax.set_yticks(range(len(groups)), [f'{k}\nsource: {s}' for k, s in groups])
    ax.set_title(f'Conditional compatibility decisions | Fixed reference: {ref}\n{mode}; AC-only excludes AA; no genome-wide/full-Suite claim'+(' | SYNTHETIC TEST DATA' if synthetic else ''), fontsize=10)
    fig.text(.5, .01, "Missing margins, unavailable results or unverified source calibration prevent a combined compatibility claim. Single-sample similarity is not equivalence.", ha="center", fontsize=7)
    fig.tight_layout(rect=[0, .06, 1, 1])
    for extension in ("png", "pdf", "svg"): fig.savefig(out / f"compatibility_table.{extension}", dpi=170, bbox_inches="tight")
    plt.close(fig)
    if mode == "with_bam" and not any(r.get("metric") == GLOBAL_METRIC for r in upstream or []):
        _local_cnv_cohort_plot(out, matched_cohort or [], ref, synthetic)
    _upstream_ladder(out, upstream or [], ref, mode, synthetic)
    return _sample_ladders(out, sample_rows or [], ref, mode, synthetic, upstream_samples)


def _percentage_point_interval(row, name):
    value = row[name]
    if value is None:
        return '不可评价'
    low, high = row[name + '_ci_low'], row[name + '_ci_high']

    def percent_point(value):
        rounded = round(100 * value, 1)
        return f'{rounded:+.1f}' if rounded else '0.0'

    return f'{percent_point(value)} [{percent_point(low)}, {percent_point(high)}]'


def _report_effects(downstream, upstream, mode, independence_unit):
    unit = '患者' if independence_unit == 'patient' else '样本'
    header = ["## 配对变化概览", "",
              f"Δ 为目标版本减来源版本，单位为百分点。方括号给出保守置信区间。配对在同一样本、同一评价单元内进行，统计以{unit}为独立单位。",
              "", "| 分析 | 来源 → 目标 | 检出 Δ [区间] | 置信水平 | 新增 / 遗漏 | 分类改变 / 可比较 | 可评价 / 应评价 | 独立单位数 | 判定范围 |",
              "|---|---|---:|---:|---:|---:|---:|---:|---|"]
    for row in downstream:
        if row['identity_anchor']:
            continue
        extra, missed = (row[k] for k in ('extra_detection_count', 'missed_detection_count'))
        event_counts = f'{extra} / {missed}' if extra is not None else '—'
        classes = f"{row['class_change_count']} / {row['class_change_n_regions']}" if row['class_change_n_regions'] else '—'
        if row['class_change_n_regions'] and not row['full_label_comparable']:
            classes += '（共同类别）'
        level = f"{100*row['confidence_level']:g}%" if row['detection_difference'] is not None else '—'
        scope = '主要比较' if row['primary_pair'] else '探索性'
        header.append(f"| {KINDS[row['analysis_kind']]} | {row['source_version']} → {row['target_version']} | "
                      f"{_percentage_point_interval(row, 'detection_difference')} | {level} | {event_counts} | {classes} | "
                      f"{row['detection_difference_n_regions']} / {row['expected_regions']} | {row['detection_difference_n_independent_units']} | {scope} |")
    header += ["", "检出 Δ 使用该版本对的全部可评价单元。新增和遗漏是原始单元计数。分类改变限于双方检出且至少有共同可比类别的单元。各版本对的缺失或有信号分母可能不同，不能把这些结果直接当作共同队列的排名。"]
    if mode == 'with_bam':
        header += ["", "| 来源 → 目标 | 位点保留 Δ [区间] | 置信水平 | 双方可配对 / 应评价 | 独立单位数 |",
                   "|---|---:|---:|---:|---:|"]
        for row in upstream:
            if row['source_version'] == row['target_version']:
                continue
            level = f"{100*row['confidence_level']:g}%" if row['retention_difference'] is not None else '—'
            header.append(f"| {row['source_version']} → {row['target_version']} | "
                          f"{_percentage_point_interval(row, 'retention_difference')} | {level} | "
                          f"{row['retention_difference_n_seeds']} / {row['expected_seeds']} | {row['retention_difference_n_independent_units']} |")
        header += ["", "保留 Δ 仅使用双方已完成且候选坐标、参考策略和匹配规则相同的候选。缺失、失败或不匹配的候选保留在应评价分母中。"]
    return header


def write_concordance_reports(out, downstream_records, upstream_records, config, *, save_records=True, local_cnv_records=None):
    """Write local conditional reports, preserving unknown and unavailable states.

    ``comparison_eligible`` must be True. ``samples[].upstream_seed_ids`` declares
    the entire known upstream candidate universe, including [] for no candidates.
    Without it, observed rows are described but upstream agreement is not accepted.
    Combined with-BAM acceptance also requires passed source bridges and local
    independent same-version source control, comparable candidate records and
    native_recipe_match=True on both sides of the local paired comparison.
    AC-only records require a common nonempty source_graph_identity per unit.
    """
    ref, targets, samples, options, alpha, minimum, margin = _configuration(config)
    global_mode = config["mode"] == "with_bam" and config.get("cnv_scope") == "global"
    downstream_records, upstream_records = list(downstream_records), list(upstream_records)
    kinds = {r.get("analysis_kind") for r in downstream_records} or {"aa_ac" if config["mode"] == "with_bam" else "ac_only"}
    if any(k not in KINDS for k in kinds): raise ValueError("analysis_kind must be aa_ac or ac_only")
    kinds = sorted(kinds)
    groups = _prepare_downstream(downstream_records, samples, targets, ref, kinds, config["mode"])
    synthetic = bool(config.get("synthetic", config.get("_synthetic", False)))
    patient_ids = options['patient_ids']
    unit = 'patient' if patient_ids is not None else 'sample'
    # One family for each kind/source stratum plus the upstream family. The old
    # analyzer controls its version pairs and four primary metrics internally.
    group_alpha = alpha / (len(groups)+1)
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    downstream, sample_rows, links = _downstream(out, groups, targets, ref, options, group_alpha, config["mode"], synthetic, global_mode)
    seeds, upstream_samples, upstream = _upstream(upstream_records, samples, targets, ref, alpha, group_alpha, minimum, margin,
                                                  patient_ids, options.get('primary_pairs'), global_mode, config.get('seed_match_fraction', 1.0))
    from .concordance_downstream import write_downstream_similarity
    from .concordance_result import write_result
    from .concordance_upstream_similarity import upstream_similarity_rows, upstream_similarity_cohort_rows
    down_similarity = write_downstream_similarity(out, groups, targets, ref, patient_ids, synthetic, alpha)
    up_pairs = []
    if config['mode'] == 'with_bam':
        up_pairs = upstream_similarity_rows(seeds, upstream_samples, ref, config.get('seed_match_fraction', 1.0))
        for row in up_pairs:
            row['synthetic'] = synthetic
        up_cohort = upstream_similarity_cohort_rows(up_pairs, patient_ids, alpha)
        for row in up_cohort:
            row['synthetic'] = synthetic
        write_tsv(out / 'upstream_similarity_pairs.tsv', up_pairs)
        write_tsv(out / 'upstream_similarity_cohort.tsv', up_cohort)
        source_pairs = {(r['sample_id'], r['target_version']): r for r in up_pairs
                        if r['source_version'] == r['comparator_version']}
        for sample_row in upstream_samples:
            pair = source_pairs[sample_row['sample_id'], sample_row['target_version']]
            sample_row.update(upstream_interval_jaccard=pair['upstream_interval_jaccard'],
                              upstream_similarity_status=pair['upstream_similarity_status'])
    from .concordance_local_cnv_report import write_local_cnv_comparator
    matched = dict(samples=[], seeds=[], cohort=[], status="not_applicable_global_workflow") if global_mode else write_local_cnv_comparator(out, local_cnv_records, config, patient_ids, alpha, synthetic)
    matched_by_suite = {(row['sample_id'], suite): row for row in matched['samples']
                        for suite in row['target_suite_versions']}
    plot_upstream_samples = [dict(row, source_matched_cnv=matched_by_suite.get((row['sample_id'], row['target_version'])))
                             for row in upstream_samples]
    result = write_result(out, down_similarity['samples'], up_pairs, config['mode'], ref, synthetic, matched['samples'], global_mode=global_mode)
    combined = _combine(downstream, upstream, config["mode"], samples, global_mode)
    from .concordance_cohort_plot import cohort_metric_rows, plot_cohort_compatibility
    cohort_rows = cohort_metric_rows(downstream, combined, upstream)
    for row in cohort_rows: row['synthetic'] = synthetic
    write_tsv(out / 'cohort_compatibility_metrics.tsv', cohort_rows)
    for filename, table in (("downstream_summary.tsv", downstream), ("sample_similarity.tsv", sample_rows),
                            ("upstream_seeds.tsv", seeds), ("upstream_samples.tsv", upstream_samples),
                            ("upstream_statistics.tsv", upstream), ("compatibility_combined.tsv", combined), ("downstream_report_index.tsv", links)):
        for row in table: row["synthetic"] = synthetic
        write_tsv(out / filename, table)
    _json(out / "config.json", config)
    if save_records:
        _json(out / "downstream_records.json", downstream_records)
        _json(out / "upstream_records.json", upstream_records)
    plot_error, sample_plots = None, []
    try:
        figures = _plots(out, downstream, combined, ref, config["mode"], synthetic, sample_rows, upstream, plot_upstream_samples, matched["cohort"])
        if isinstance(figures, list): sample_plots = figures
        plot_cohort_compatibility(out, cohort_rows, ref, config['mode'], synthetic)
    except ImportError as exc: plot_error = str(exc)
    write_tsv(out / "sample_plot_index.tsv", sample_plots)
    sample_index = ['# Per-sample observed comparisons', '',
                    f'Fixed reference: {ref}. Single samples do not establish statistical equivalence.', '']
    sample_index += [f"- {r['sample_id']} / {r['analysis_kind']} / source {r['source_version']}: "
                     f"[PNG](<{r['png']}>) · [PDF](<{r['pdf']}>) · [SVG](<{r['svg']}>)" for r in sample_plots]
    (out / 'sample_reports.md').write_text('\n'.join(sample_index)+'\n', encoding='utf-8')
    summary = dict(reference_version=ref, mode=config["mode"], selected_samples=len(samples), target_versions=targets,
                   analysis_kinds=kinds, source_versions=sorted({s["source_version"] for s in samples.values()}),
                   synthetic=synthetic, compatibility_combined=combined, result=result,
                   downstream_similarity_files=down_similarity["files"],
                   local_cnv_comparator=matched,
                   downstream=downstream, upstream=upstream, family_alpha=alpha, allocated_group_alpha=group_alpha,
                   independence_unit=unit, selected_independent_units=len(set(patient_ids.values())) if patient_ids is not None else len(samples),
                   primary_pairs=options.get('primary_pairs'), plot_error=plot_error,
                   sample_plot_index=sample_plots, cohort_compatibility_metrics=cohort_rows,
                   scope=GLOBAL_SCOPE if global_mode else "conditional known-candidate comparison; never complete Suite or genome-wide compatibility",
                   report=str(out / "report.md"))
    _json(out / "summary.json", summary)
    lines = ["# AmpliconSuite Cross-Version Concordance", "", f"固定参考：**{ref}**；模式：`{config['mode']}`；输入样本：{len(samples)}。",
             "", "**人工合成测试数据：不形成患者或生物学结论。**" if synthetic else "本次比较覆盖各版全局种子和抽样位点的重建结果。" if global_mode else "这是所提供输入、已知候选区域和已评价指标的条件性比较。", ""]
    if options.get('primary_pairs') == []:
        lines += ["**探索性分析：未指定主要版本对，不作正式兼容判定。**", ""]
    else:
        outcomes = Counter(row['decision'] for row in combined if row['primary_pair'])
        if outcomes:
            lines += [f"**主要比较结论：** {outcomes['within_margin']} 项支持预设范围内一致性，"
                      f"{outcomes['exceeds_margin']} 项超出预设范围，"
                      f"{outcomes['insufficient_evidence']} 项证据不足。每项按分析类型、来源和目标计数。", ""]
    lines += ["[查看逐样本结果、上下游相似性和标签一致率](result.md) · [下载 result.tsv](result.tsv)", ""]
    if config['mode'] == 'with_bam' and not global_mode:
        lines += ["输入样本原用 CNVkit 的局部重跑是该样本的 CNVkit 对照。各版共用输入来源的 AA 筛种器。"
                  "[逐样本 CNVkit 对照](local_cnv_comparator_samples.tsv) · [按患者或样本汇总](local_cnv_comparator_cohort.tsv) · "
                  "[逐候选保留与区段](local_cnv_comparator_seeds.tsv)", ""]
    if config["mode"] == "with_bam" and not global_mode and (out / "local_cnv_comparator_cohort.png").is_file():
        lines += ["![局部 CNVkit 版本对照与描述性区间](local_cnv_comparator_cohort.png)", ""]
    lines += down_similarity["report_lines"] + [""]
    if config["mode"] == "with_bam":
        lines += ["", "### 来源对照", "",
                  "| 样本 | 来源 Suite | 原结果与同版重建 | 输入来源 CNVkit | 来源 Suite 配方实际 CNVkit | 来源种子校准 |",
                  "|---|---|---|---|---|---|"]
        controls = {(r["sample_id"], r["target_version"]): r for r in upstream_samples}
        for sid, spec in samples.items():
            control = controls[sid, spec["source_version"]]
            versions = ", ".join(control["cnvkit_versions"]) or "缺失"
            calibration = ", ".join(control["calibration_statuses"]) or "unverified"
            lines.append(f"| {sid} | {spec['source_version']} | {spec.get('source_bridge_status', 'unverified')} | "
                         f"{spec.get('source_cnvkit_version') or '未记录'} | {versions} | {calibration} |")
        if global_mode:
            lines += ["", "对照本样本原始种子 BED 与来源配方重新筛选的全套种子，另核对抽中位点的原始与重建 AA/AC 结果。来源版本或配方缺少证据时保留为未核实。", ""]
        else:
            lines += ["", "来源复现有差异时，先对照原始和重建的 AA 区段、graph、cycles 及 AC feature BED，再固定同一套图做同版本 AC 参数对照。"
                      "局部 CNV 校准先用已保存的 CNR、CNS 在相同可测坐标配对；再单独使用已核实的原 CNVkit 版本、参考策略和种子筛选规则做来源对照。"
                      "校准状态及原因可在 [逐种子记录](upstream_seeds.tsv) 中核对。", ""]
    lines += _report_effects(downstream, upstream, config['mode'], unit)
    if global_mode:
        lines += ["", "## 全局种子与下游的评价范围", "",
                  "各版按原生配方完成全基因组 CNV 分段和种子筛选。报告把成功版本至少有一版发现的种子按严格重叠连成共同位点；边界相邻的种子分属不同位点。共享位点可容纳种子的拆分或合并。报告另列逐版原生种子条数、另一版种子区段覆盖达到 seed_match_fraction 的条数与比例，以及完整区段覆盖碱基 Jaccard。",
                  "", "上游主指标取两个方向原生种子未保留率的较大值。每个方向逐条检查本版原生种子被另一版种子覆盖的比例是否达到 seed_match_fraction。仅一版有种子时失配为 100%，双方都空或任一方失败时留空。upstream_loss_margin 约束这个失配值，位点存在差和区段 Jaccard 单独作为描述指标。",
                  "", "下游从所有版本共同的位点并集中抽样，各版使用其落在抽中位点的原生种子运行 AA＋AC，保留延伸至其他位点的完整重建结构。下游差异可来自筛种、重建及分类。失败或缺失记录保留为未评价。",
                  "", "联合判断要求原数据来源复现及配方证据。观察范围是所分析样本、实际完成的全局筛种以及抽中位点的下游指标。基因组种子发现的生物学准确性需要独立真值。"]
    else:
        lines += [
                 "", "各表按输入来源 Suite 分层，比较该来源与每个目标版本；同时在分层子报告中保留固定参考和其他版本间比较。输入来源自身只作身份锚点，不能作为独立验证。",
                 "", "## 下游与上游分别解读", "", "with_bam 的 AA＋AC 使用同一组输入种子；results_only 的 AA＋AC 项只比较既有输出指标。仅 AC 使用相同来源 AA 输出，不能称为 AA 重建兼容。源数据与重跑数据的输入范围必须已核实，`comparison_eligible` 不明确为 true 的结果保留为不可评价。AC-only 还须记录一致的 source_aa_version 和 source_graph_identity；相同 AA 版本但不同图文件不构成固定输入。",
                 "", "上游同时报告原候选内两版最终种子区段的 Jaccard，以及目标版本的保留个数/可评价候选数。未保留比例是主失配指标；retention_fraction 仅解释覆盖和匹配。匹配规则不是兼容阈值。候选全集未声明时无法发现全部版本共同漏记的种子，不能给出上游兼容结论。",
                 "", "with_bam 的联合结论要求上游、固定种子 AA＋AC 均在预设范围内。上游还要求独立的同版本来源复现通过、来源和目标的 CNVkit 配方匹配，并对全部候选核对相同坐标、参考策略和匹配规则。目标版本的保留结果不作为校准依据，目标漏种仍可支持有条件的超差结论。",
                 "", "with_bam 还要求每例原导入结果与同一来源 Suite 的种子重跑结果完成来源桥接，source_bridge_status 必须为 passed；否则原数据的联合判定无论接受或超差均为证据不足，分项观察仍保留。桥接通过仅表示所检查的检出、分类与区段指标一致，不是生物学真实性或完整结构相同的证明。既有结果未证明参数、资源及运行条件相同，观察差异不能自动归因于 AA 算法版本。",
                 "", "上游配方受控或缺少可核查的独立来源控制时，仍报告描述性保留率和配对差。来源桥接通过后，可靠的 AA＋AC 指标超差可以独立支持该下游范围的超差。decision_basis 指明本次判定依据。联合结果覆盖局部候选保留和固定种子 AA＋AC 两个步骤。",]
    lines += [
             "", "## 统计与证据边界", "", "先对样本内种子取平均。有患者标识时，对患者内样本取平均并按患者等权。未提供患者标识时按样本等权，前提是样本独立。区间保留零观察失配时的不确定性。各分析分层与上游分配误差预算，主要版本对的判定另作 Bonferroni 校正。探索性比较报告观察差异和区间。",
             "", "容忍阈值未预设、独立样本不足、未知/失败或来源校准未验证时，联合接受结论为证据不足；P>0.05不表示兼容。单样本表只描述观察相似度，不给等效判断。全阴性区域不计为结构一致证据，有信号结果单列。",
             "", "[逐样本相似度](sample_similarity.tsv) · [下游汇总](downstream_summary.tsv) · [上游逐种子](upstream_seeds.tsv) · [上游逐样本](upstream_samples.tsv) · [上游统计](upstream_statistics.tsv) · [联合兼容表](compatibility_combined.tsv)", "",
             "![版本相似性天梯图](similarity_ladder.png)", "", "![整组统计兼容图](cohort_compatibility.png)", "",
             "整组统计兼容图显示差异大小、区间、预设阈值、独立单位数与联合判定；[绘图数据](cohort_compatibility_metrics.tsv)。主要比较使用与判定相同的多重校正单侧置信界，探索性比较仅显示描述性区间。", "",
             "![限定范围兼容表](compatibility_table.png)", "", "## 分层详细报告", ""]
    lines += ["逐样本天梯图见 [全部样本图](sample_reports.md) 和 [TSV索引](sample_plot_index.tsv)：每个样本×分析类型一张，标记输入来源与固定参考；灰色不可评价项不补为零。单样本只显示观察相似度，不给置信区间或等效判断。CLASS* 表示类别支持不同或未知，相关标签与类别区段指标灰显；表中保留 full_label_comparable、共同支持类别及限定类别范围的观察值，不能把它们解释为完整标签100%一致。", ""]
    lines += [f"- {r['sample_id']} / {r['analysis_kind']}：[PNG](<{r['png']}>) · [PDF](<{r['pdf']}>) · [SVG](<{r['svg']}>)" for r in sample_plots[:20]]
    if config["mode"] == "with_bam":
        if global_mode:
            lines += ["", "![两个方向较低的原生种子保留率](upstream_retention_ladder.png)", ""]
        else:
            lines += ["", "![已知候选种子保留率](upstream_retention_ladder.png)", "", "上游图显示候选保留率及独立单位等权区间。SRC?/SRC! 标明来源对照未知或失败，PAIR? 标明配对候选范围不完整，CTRL/NATIVE? 标明受控或未验证配方。", ""]
    lines += [f"- [{r['analysis_kind']} / 输入来源 {r['source_version']}](<{r['path']}>)" for r in links]
    if plot_error: lines += ["", "绘图依赖缺失；表格已保存：" + plot_error]
    (out / "report.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    return summary
