"""Reports for conditional AmpliconSuite benchmarks; no BAM execution or export.

Input-source comparisons, fixed-reference comparisons, and AC-only experiments
remain distinct. Every independent unit is a sample, never an individual seed.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from .reporting import write_tsv
from .statistics import AGREEMENT_METRICS, VALID, _pair, analyze


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


def _bounds(values, alpha, two_sided=True):
    if not values:
        return None, None
    estimate = _mean(values)
    radius = math.sqrt(math.log((2 if two_sided else 1) / alpha) / (2 * len(values)))
    return max(0., estimate-radius), min(1., estimate+radius)


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
    options = dict(config.get("statistics", {}))
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


def _oriented_stat(stats, target, source, metric):
    direct = next((s for s in stats if s["version"] == target and s["comparator"] == source
                   and s["subset"] == "signal_positive" and s["metric"] == metric), None)
    if direct is not None:
        return direct
    # The four displayed similarities are symmetric; all directed event details
    # below are generated explicitly as target versus the input source.
    return next((s for s in stats if s["version"] == source and s["comparator"] == target
                 and s["subset"] == "signal_positive" and s["metric"] == metric), None)


def _downstream(out, groups, targets, ref, options, group_alpha, mode, synthetic):
    comparisons, sample_rows, links = [], [], []
    for (kind, source), records in groups.items():
        folder = out / "downstream" / kind / _slug(source)
        folder.mkdir(parents=True, exist_ok=True)
        settings = dict(options, alpha=group_alpha, analysis_scope="fixed_seed",
                        versions=list(dict.fromkeys([ref, source] + targets)))
        settings.pop("upstream_loss_margin", None)
        pairs, stats, compatibility = analyze(records, ref, settings)
        # The shared calculation is valid for these endpoints, but its AA+AC
        # wording is replaced for AC-only output, never reused as an AA claim.
        label = _kind_label(kind, mode)
        scope = label + "; conditional on selected input units; no genome-wide seed-discovery validation"
        endpoint = "ac_only_fixed_source_aa" if kind == "ac_only" else "selected_seed_aa_ac" if mode == "with_bam" else "existing_aa_ac_output_metrics"
        unit = "sample_seed" if mode == "with_bam" and kind == "aa_ac" else "sample_input_result"
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
            + ("仅重做 AC 分类，不能称为 AA 重建兼容。\n\n" if kind == "ac_only" else "相同输入种子的 AA＋AC 比较，未验证全基因组种子发现。\n\n" if mode == "with_bam" else "既有 AA＋AC 输出指标比较；未证明历史参数和资源相同，差异不能自动归因于算法版本。\n\n")
            + "[统计表](statistics.tsv) · [差异明细](differences.tsv) · [兼容表](compatibility.tsv)。同一样本内先平均，再对样本等权；单样本只描述观察相似度。\n", encoding="utf-8")
        links.append(dict(analysis_kind=kind, source_version=source, path=str((folder / "report.md").relative_to(out))))
        lookup = {(r["sample_id"], r["seed_id"], r["version"]): r for r in records}
        units = sorted({(r["sample_id"], r["seed_id"]) for r in records})
        for target in list(dict.fromkeys([source] + targets)):
            target_pairs = [_pair(lookup[s, unit, target], lookup[s, unit, source], ref, "vs_input_source") for s, unit in units]
            native = next((r for r in compatibility if {r["version"], r["comparator"]} == {target, source}), None) if target != source else None
            row = dict(reference_version=ref, source_version=source, target_version=target, analysis_kind=kind,
                       comparison_scope=scope, identity_anchor=target == source,
                       confidence_level=1-group_alpha,
                       decision=native["decision"] if native else "insufficient_evidence",
                       reasons=native["reasons"] if native else ["identity_is_not_validation"],
                       n_selected_samples=len({s for s, _ in units}), signal_samples=native["signal_samples"] if native else None,
                       eligible_units=native["evaluable_regions"] if native else None,
                       excluded_units=native["excluded_regions"] if native else None,
                       compatibility_metric_scope=native["compatibility_metric_scope"] if native else {"required_metrics": list(AGREEMENT_METRICS)})
            row.update(_label_scope([p for p in target_pairs if p["evaluable"]]))
            for name, _ in SIMILARITIES:
                metric = "label_discordance" if name == "label_agreement" else name
                stat = _oriented_stat(stats, target, source, metric) if target != source else None
                value, lo, hi = ((stat[k] for k in ("estimate", "conservative_ci_low", "conservative_ci_high")) if stat else (None, None, None))
                if name == "label_agreement" and value is not None:
                    value, lo, hi = 1-value, 1-hi, 1-lo
                row.update({name: value, name+"_ci_low": lo, name+"_ci_high": hi, name+"_n_samples": stat["n_samples"] if stat else 0})
            if target in targets: comparisons.append(row)
            for sample in sorted({s for s, _ in units}):
                paired = [p for p in target_pairs if p["sample_id"] == sample]
                valid = [p for p in paired if p["evaluable"] and p["is_seed_region"]]
                signal = [p for p in valid if p["signal_positive"]]
                state = "evaluable" if valid and len(valid) == len(paired) else "partial" if valid else "not_evaluable"
                item = dict(sample_id=sample, reference_version=ref, source_version=source, target_version=target,
                            analysis_kind=kind, comparison_scope=scope, identity_anchor=target == source,
                            status=state, n_units=len(paired), n_evaluable_units=len(valid), n_signal_units=len(signal),
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


def _upstream(records, samples, targets, ref, alpha, decision_alpha, minimum, margin):
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
        lookup[sample, version, seed] = row
        if row["status"] != "no_seeds":
            observed[sample].add(seed)
    universe = {sample: samples[sample].get("upstream_seed_ids", sorted(observed[sample])) for sample in samples}
    for sample in samples:
        if "upstream_seed_ids" in samples[sample] and observed[sample] - set(universe[sample]):
            raise ValueError("Upstream records contain seeds outside the declared candidate universe")
    seed_rows, sample_rows, statistics = [], [], []
    sources = sorted({s["source_version"] for s in samples.values()})
    family_count = max(1, sum(t != source for source in sources for t in targets))
    adjusted_alpha = decision_alpha / family_count
    for source in sources:
        selected = sorted(s for s in samples if samples[s]["source_version"] == source)
        for target in targets:
            local = []
            for sample in selected:
                ids = universe[sample]
                rows = []
                for seed in ids or ["__no_seeds__" if "upstream_seed_ids" in samples[sample] else "__missing_seed_manifest__"]:
                    original = lookup.get((sample, target, seed))
                    status = "no_seeds" if not ids and "upstream_seed_ids" in samples[sample] else "missing"
                    row = dict(original or dict(sample_id=sample, version=target, seed_id=seed, source_version=source,
                               status=status, reason="No declared input candidates" if status == "no_seeds" else "No upstream result for expected candidate",
                               comparison_eligible=False, calibration_status="unverified", retention_fraction=None))
                    row.update(reference_version=ref, target_version=target, candidate_universe_declared="upstream_seed_ids" in samples[sample])
                    row["evaluable"] = (row["status"] in {"retained", "not_retained"} and row.get("comparison_eligible") is True
                                        and bool(row.get("matching_rule")) and bool(row.get("cnvkit_version")))
                    row["seed_loss"] = float(row["status"] == "not_retained") if row["evaluable"] else None
                    if row["status"] in {"retained", "not_retained"} and not row["evaluable"]:
                        row["evaluation_reason"] = "Missing verified comparison eligibility, matching rule or actual CNVkit version"
                    rows.append(row)
                    seed_rows.append(row)
                valid = [r for r in rows if r["evaluable"]]
                fractions = [float(r["retention_fraction"]) for r in valid if r.get("retention_fraction") is not None]
                item = dict(sample_id=sample, reference_version=ref, source_version=source, target_version=target,
                            expected_seeds=len(ids), evaluable_seeds=len(valid), unavailable_seeds=len(ids)-len(valid),
                            no_seeds=not ids and "upstream_seed_ids" in samples[sample],
                            status_counts=dict(Counter(r["status"] for r in rows)),
                            seed_loss=_mean([r["seed_loss"] for r in valid]),
                            mean_retention_fraction=_mean(fractions),
                            candidate_universe_declared="upstream_seed_ids" in samples[sample],
                            calibration_passed=bool(valid) and len(valid) == len(ids) and all(r.get("calibration_status") == "passed" for r in valid),
                            calibration_statuses=sorted({r.get("calibration_status", "unverified") for r in rows}),
                            native_recipe_match=bool(valid) and len(valid) == len(ids) and all(r.get("native_recipe_match") is True for r in valid),
                            native_recipe_statuses=sorted({"matched" if r.get("native_recipe_match") is True else "controlled" if r.get("native_recipe_match") is False else "unverified" for r in rows}),
                            cnvkit_versions=sorted({str(r["cnvkit_version"]) for r in valid}),
                            interpretation="descriptive_known_candidate_retention; not_genome_wide_sensitivity")
                local.append(item)
                sample_rows.append(item)
            values = [r["seed_loss"] for r in local if r["seed_loss"] is not None]
            low, high = _bounds(values, alpha)
            decision_low, decision_high = _bounds(values, adjusted_alpha, two_sided=False)
            reasons = []
            if margin is None: reasons.append("upstream_loss_margin_not_predefined")
            if len(values) < minimum: reasons.append("too_few_independent_samples")
            if any(r["unavailable_seeds"] or r["seed_loss"] is None for r in local): reasons.append("unavailable_or_no_candidate_samples")
            if any(not r["candidate_universe_declared"] for r in local): reasons.append("candidate_universe_not_declared")
            if any(len(r["cnvkit_versions"]) > 1 for r in local): reasons.append("multiple_cnvkit_versions_within_sample_target")
            if len({v for r in local for v in r["cnvkit_versions"]}) > 1: reasons.append("cnvkit_version_differs_across_samples")
            decision = "insufficient_evidence"
            if not reasons and decision_high < float(margin): decision = "within_margin"
            elif margin is not None and len(values) >= minimum and decision_low > float(margin): decision = "exceeds_margin"
            else: reasons.append("confidence_bounds_do_not_establish_agreement")
            if source == target:
                decision = "insufficient_evidence"
                reasons.append("identity_is_not_validation")
            statistics.append(dict(reference_version=ref, source_version=source, target_version=target,
                metric="known_candidate_not_retained_fraction", estimate=_mean(values),
                conservative_ci_low=low, conservative_ci_high=high, confidence_level=1-alpha,
                agreement_lower=decision_low, agreement_upper=decision_high, agreement_alpha=adjusted_alpha,
                n_samples=len(values), selected_samples=len(selected), expected_seeds=sum(r["expected_seeds"] for r in local),
                evaluable_seeds=sum(r["evaluable_seeds"] for r in local), no_seed_samples=sum(r["no_seeds"] for r in local),
                upstream_loss_margin=margin, min_samples=minimum, decision=decision, reasons=sorted(set(reasons)),
                calibration_passed=all(r["calibration_passed"] for r in local),
                calibration_statuses=sorted({c for r in local for c in r["calibration_statuses"]}),
                native_recipe_match=all(r["native_recipe_match"] for r in local),
                native_recipe_statuses=sorted({v for r in local for v in r["native_recipe_statuses"]}),
                p_agreement_bound=math.exp(-2*len(values)*(float(margin)-_mean(values))**2) if values and margin is not None and _mean(values) < float(margin) else 1. if values and margin is not None else None,
                comparison_scope="known_candidate_seed_retention_only; sample_equal_weight; no_genome_wide_detection_validation"))
    return seed_rows, sample_rows, statistics


def _combine(downstream, upstream, mode, samples):
    lookup = {(r["source_version"], r["target_version"]): r for r in upstream}
    result = []
    for down in downstream:
        up = lookup.get((down["source_version"], down["target_version"]))
        includes_upstream = mode == "with_bam" and down["analysis_kind"] == "aa_ac"
        decision, reasons = down["decision"], list(down["reasons"])
        basis = ["downstream_metrics_exceed_margin"] if decision == "exceeds_margin" else ["downstream_metrics_within_margin"] if decision == "within_margin" else []
        scope = "conditional_local_upstream_and_fixed_seed_aa_ac" if includes_upstream else "ac_only_on_fixed_source_aa" if down["analysis_kind"] == "ac_only" else "provided_downstream_aa_ac_results_only"
        if includes_upstream:
            upstream_validated = bool(up and up["calibration_passed"] and up["native_recipe_match"])
            basis = ["downstream_metrics_exceed_margin"] if down["decision"] == "exceeds_margin" else []
            if up and up["decision"] == "exceeds_margin" and upstream_validated:
                basis.append("calibrated_native_upstream_exceeds_margin")
            if basis:
                decision = "exceeds_margin"
            elif down["decision"] == "within_margin" and up and up["decision"] == "within_margin" and upstream_validated:
                decision = "within_margin"
                basis = ["downstream_and_calibrated_native_upstream_within_margin"]
            else:
                decision = "insufficient_evidence"
            reasons += up["reasons"] if up else ["upstream_results_missing"]
            if up and up["decision"] == "exceeds_margin" and not upstream_validated:
                reasons.append("upstream_exceedance_not_validated_for_historical_inference")
            if not up or not up["calibration_passed"]: reasons.append("local_upstream_calibration_not_passed")
            if not up or not up["native_recipe_match"]:
                reasons.append("controlled_CNVkit_not_historical_recipe" if up and "controlled" in up["native_recipe_statuses"] else "native_CNVkit_recipe_not_verified")
        bridge_statuses = sorted({s.get("source_bridge_status", "unverified") for s in samples.values()
                                  if s["source_version"] == down["source_version"]})
        bridge_passed = bool(bridge_statuses) and bridge_statuses == ["passed"]
        if mode == "with_bam" and not bridge_passed:
            decision, basis = "insufficient_evidence", []
            reasons.append("input_source_to_rerun_bridge_not_passed")
        result.append(dict(reference_version=down["reference_version"], source_version=down["source_version"],
            target_version=down["target_version"], analysis_kind=down["analysis_kind"], mode=mode,
            comparison_scope=scope, decision=decision, decision_basis=basis, reasons=sorted(set(reasons)),
            downstream_decision=down["decision"], upstream_decision=up["decision"] if includes_upstream and up else "not_assessed",
            upstream_calibration_passed=up["calibration_passed"] if includes_upstream and up else None,
            native_recipe_match=up["native_recipe_match"] if includes_upstream and up else None,
            source_bridge_passed=bridge_passed if mode == "with_bam" else None,
            source_bridge_statuses=bridge_statuses if mode == "with_bam" else [],
            n_selected_samples=down["n_selected_samples"], signal_samples=down["signal_samples"],
            compatibility_metric_scope=down["compatibility_metric_scope"],
            not_established=["genome_wide_seed_discovery", "complete_suite_interchangeability", "biological_truth", "cycle_order_or_breakpoint_topology"]
                + (["aa_reconstruction_compatibility"] if down["analysis_kind"] == "ac_only" else [])))
    return result


def _sample_ladders(out, sample_rows, ref, mode, synthetic):
    import matplotlib.pyplot as plt
    groups = defaultdict(list)
    for row in sample_rows: groups[row["sample_id"], row["analysis_kind"], row["source_version"]].append(row)
    folder = out / "sample_ladders"; folder.mkdir(exist_ok=True)
    indices = []
    codes = {"failed": "F", "missing": "M", "unsupported": "U", "invalid": "I", "interrupted": "X",
             "blocked": "B", "no_seeds": "NS", "not_run": "NR"}
    for (sample, kind, source), rows in groups.items():
        fig, axes = plt.subplots(1, 4, figsize=(13, max(4.0, .38*len(rows)+2.3)), sharey=True)
        labels = [r["target_version"] + (" [source]" if r["target_version"] == source else "")
                  + (" [ref]" if r["target_version"] == ref else "") for r in rows]
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
        axes[0].set_yticks(range(len(rows)), labels)
        fig.suptitle(f'Sample: {sample} | Input-source comparator: {source} | Fixed reference: {ref}\n'
                     + _kind_label(kind, mode) + " | observed signal-positive input units"
                     + (' | SYNTHETIC TEST DATA' if synthetic else ''), fontsize=10)
        fig.text(.5, .015, "One sample: observed similarity only; no confidence interval, compatibility probability or equivalence decision. Source diamond = identity.\nGray: F failed, M missing, U unsupported, I invalid, X interrupted, B blocked, NS no seeds, NR not run, PART incomplete, ND/ND both negative, NA unavailable; SRC = input source.\nCLASS*: class support differs or is unknown; shared-class values remain in the table but do not represent complete classification agreement.", ha="center", fontsize=6.5)
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
    for ax, source in zip(axes.flat, sources):
        rows = [r for r in upstream if r["source_version"] == source]
        for y, row in enumerate(rows):
            if row["estimate"] is None:
                ax.axhspan(y-.33, y+.33, color="#eeeeee", zorder=0)
                ax.text(.5, y, "not evaluable", ha="center", va="center", color="#777777", fontsize=8)
                continue
            value, low, high = 1-row["estimate"], 1-row["conservative_ci_high"], 1-row["conservative_ci_low"]
            checked = row["calibration_passed"] and row["native_recipe_match"]
            flags = []
            if not row["calibration_passed"]: flags.append("CAL!" if "failed" in row["calibration_statuses"] else "CAL?")
            if not row["native_recipe_match"]: flags.append("CTRL" if "controlled" in row["native_recipe_statuses"] else "NATIVE?")
            ax.errorbar(value, y, xerr=[[max(0, value-low)], [max(0, high-value)]], fmt="o" if checked else "^",
                        color="#4178A9" if checked else "#C5884E", markersize=5, capsize=3)
            ax.text(.02, y-.19, f'n={row["n_samples"]}' + ('; ' + ", ".join(flags) if flags else ''), fontsize=7, color="#666666")
        ax.set_yticks(range(len(rows)), [r["target_version"] + (" [source re-evaluation]" if r["target_version"] == source else "")
                                       + (" [ref]" if r["target_version"] == ref else "") for r in rows])
        ax.set_ylim(len(rows)-.5, -.5); ax.set_xlim(-.03, 1.03)
        ax.set_xticks([0, .25, .5, .75, 1], ["0", "25", "50", "75", "100"])
        ax.set_xlabel("Known candidates retained (%)"); ax.grid(axis="x", color="#eeeeee"); ax.set_axisbelow(True)
        ax.set_title(f"Input-source candidate set: {source}", fontsize=9, loc="left")
    confidence = upstream[0]["confidence_level"]*100
    fig.suptitle(f'Known-candidate retention | Fixed reference: {ref}\nBounded local upstream; not genome-wide detection sensitivity'
                 + (' | SYNTHETIC TEST DATA' if synthetic else ''), fontsize=11)
    fig.text(.5, .015, f'Equal-weight sample means and conservative {confidence:g}% bounds; a source-version re-evaluation is an observed local result, not an identity point.\nCAL? unverified calibration; CAL! failed calibration; CTRL controlled nonhistorical recipe; NATIVE? recipe unknown. Retention is descriptive, not compatibility.', ha="center", fontsize=6.5)
    fig.tight_layout(rect=[0, .09, 1, 1])
    for ext in ("png", "pdf", "svg"): fig.savefig(out / f"upstream_retention_ladder.{ext}", dpi=170, bbox_inches="tight")
    plt.close(fig)


def _plots(out, downstream, combined, ref, mode, synthetic, sample_rows=None, upstream=None):
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
                    ax.text(.02, y-.18, f'n={row[metric+"_n_samples"]}', fontsize=6, color="#666666")
            ax.set_xlim(-.03, 1.03); ax.set_ylim(len(rows)-.5, -.5)
            ax.set_xticks([0, .5, 1], ["0", "50", "100"]); ax.set_xlabel(title+" (%)")
            ax.set_yticks(range(len(rows)), [r["target_version"] for r in rows] if ax is axes[index, 0] else [])
            ax.grid(axis="x", color="#eeeeee"); ax.set_axisbelow(True)
        axes[index, 0].set_title(f'{_kind_label(kind, mode)}\nInput-source version comparator: {source}', loc="left", fontsize=9)
    fig.suptitle(f'AmpliconSuite Benchmarks | Fixed reference: {ref}\nSource-to-target, signal-positive selected units | {mode}'
                 + (' | SYNTHETIC TEST DATA' if synthetic else ''), fontsize=11)
    confidence = 100*downstream[0]["confidence_level"] if downstream else 95
    fig.text(.5, .01, f"Equal-weight sample means and conservative {confidence:g}% bounds; individual seeds are not independent samples.\n* Class support differs or is unknown: label/class-footprint panels are gray, not complete agreement; shared-class values are in the tables.\nAC-only panels do not assess AA. These conditional comparisons do not establish genome-wide/full-Suite compatibility.", ha="center", fontsize=7)
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
    ax.set_title(f'Conditional benchmark decisions | Fixed reference: {ref}\n{mode}; AC-only excludes AA; no genome-wide/full-Suite claim'+(' | SYNTHETIC TEST DATA' if synthetic else ''), fontsize=10)
    fig.text(.5, .01, "Missing margins, unavailable results or unverified local calibration prevent a combined compatibility claim. Single-sample similarity is not equivalence.", ha="center", fontsize=7)
    fig.tight_layout(rect=[0, .06, 1, 1])
    for extension in ("png", "pdf", "svg"): fig.savefig(out / f"compatibility_table.{extension}", dpi=170, bbox_inches="tight")
    plt.close(fig)
    _upstream_ladder(out, upstream or [], ref, mode, synthetic)
    return _sample_ladders(out, sample_rows or [], ref, mode, synthetic)


def write_benchmark_reports(out, downstream_records, upstream_records, config):
    """Write local conditional reports, preserving unknown and unavailable states.

    ``comparison_eligible`` must be True. ``samples[].upstream_seed_ids`` declares
    the entire known upstream candidate universe, including [] for no candidates.
    Without it, observed rows are described but upstream agreement is not accepted.
    Combined with-BAM acceptance also requires passed source bridges and local
    calibration, plus native_recipe_match=True on every upstream candidate.
    AC-only records require a common nonempty source_graph_identity per unit.
    """
    ref, targets, samples, options, alpha, minimum, margin = _configuration(config)
    downstream_records, upstream_records = list(downstream_records), list(upstream_records)
    kinds = {r.get("analysis_kind") for r in downstream_records} or {"aa_ac" if config["mode"] == "with_bam" else "ac_only"}
    if any(k not in KINDS for k in kinds): raise ValueError("analysis_kind must be aa_ac or ac_only")
    kinds = sorted(kinds)
    groups = _prepare_downstream(downstream_records, samples, targets, ref, kinds, config["mode"])
    synthetic = bool(config.get("synthetic", config.get("_synthetic", False)))
    # One family for each kind/source stratum plus the upstream family. The old
    # analyzer controls its version pairs and four primary metrics internally.
    group_alpha = alpha / (len(groups)+1)
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    downstream, sample_rows, links = _downstream(out, groups, targets, ref, options, group_alpha, config["mode"], synthetic)
    seeds, upstream_samples, upstream = _upstream(upstream_records, samples, targets, ref, alpha, group_alpha, minimum, margin)
    combined = _combine(downstream, upstream, config["mode"], samples)
    for filename, table in (("downstream_summary.tsv", downstream), ("sample_similarity.tsv", sample_rows),
                            ("upstream_seeds.tsv", seeds), ("upstream_samples.tsv", upstream_samples),
                            ("upstream_statistics.tsv", upstream), ("compatibility_combined.tsv", combined), ("downstream_report_index.tsv", links)):
        for row in table: row["synthetic"] = synthetic
        write_tsv(out / filename, table)
    _json(out / "config.json", config)
    _json(out / "downstream_records.json", downstream_records)
    _json(out / "upstream_records.json", upstream_records)
    plot_error, sample_plots = None, []
    try:
        figures = _plots(out, downstream, combined, ref, config["mode"], synthetic, sample_rows, upstream)
        if isinstance(figures, list): sample_plots = figures
    except ImportError as exc: plot_error = str(exc)
    write_tsv(out / "sample_plot_index.tsv", sample_plots)
    summary = dict(reference_version=ref, mode=config["mode"], selected_samples=len(samples), target_versions=targets,
                   analysis_kinds=kinds, source_versions=sorted({s["source_version"] for s in samples.values()}),
                   synthetic=synthetic, compatibility_combined=combined,
                   downstream=downstream, upstream=upstream, family_alpha=alpha, allocated_group_alpha=group_alpha,
                   independence_unit="sample", plot_error=plot_error,
                   sample_plot_index=sample_plots,
                   scope="conditional known-candidate benchmark; never complete Suite or genome-wide compatibility",
                   report=str(out / "report.md"))
    _json(out / "summary.json", summary)
    lines = ["# AmpliconSuite Benchmarks", "", f"固定参考：**{ref}**；模式：`{config['mode']}`；输入样本：{len(samples)}。",
             "", "**人工合成测试数据：不形成患者或生物学结论。**" if synthetic else "这是所提供输入、已知候选区域和已评价指标的条件性比较。",
             "", "各表按输入来源 Suite 分层，比较该来源与每个目标版本；同时在分层子报告中保留固定参考和其他版本间比较。输入来源自身只作身份锚点，不能作为独立验证。",
             "", "## 下游与上游分别解读", "", "with_bam 的 AA＋AC 使用同一组输入种子；results_only 的 AA＋AC 项只比较既有输出指标。仅 AC 使用相同来源 AA 输出，不能称为 AA 重建兼容。源数据与重跑数据的输入范围必须已核实，`comparison_eligible` 不明确为 true 的结果保留为不可评价。AC-only 还须记录一致的 source_aa_version 和 source_graph_identity；相同 AA 版本但不同图文件不构成固定输入。",
             "", "上游只回答已知候选种子是否仍入选，不代表全基因组检出率。未保留比例是主失配指标；retention_fraction 仅解释覆盖和匹配。匹配规则不是兼容阈值。候选全集未声明时无法发现全部版本共同漏记的种子，不能给出上游兼容结论。",
             "", "with_bam 的联合结论要求上游、固定种子 AA＋AC 均在预设范围内，且所有局部校准通过。两部分分别验证，不能据此宣称完整流程、未检查区域、路径顺序或断点拓扑可互换。results_only 仅报告可评价下游范围，不纳入上游兼容；AC-only 始终不评价 AA。",
             "", "with_bam 还要求每例原导入结果与同一来源 Suite 的种子重跑结果完成来源桥接，source_bridge_status 必须为 passed；否则原数据的联合判定无论接受或超差均为证据不足，分项观察仍保留。桥接通过仅表示所检查的检出、分类与区段指标一致，不是生物学真实性或完整结构相同的证明。既有结果未证明参数、资源及运行条件相同，观察差异不能自动归因于 AA 算法版本。",
             "", "上游若固定使用受控 CNVkit 而非该 Suite 的历史配置（native_recipe_match=false），仍保存描述性保留指标，但不能接受原生 Suite 的联合兼容；配置对应关系缺失时同样保留为未验证。",
             "", "局部上游超出允许范围，只有局部校准通过且使用匹配的历史 CNVkit 配置时，才可成为联合超差依据；未验证的局部近似偏差不能说明原历史数据不兼容。来源桥接通过后，可靠的 AA＋AC 指标超差可以独立支持该下游范围的超差，不要求未知上游也可评价。decision_basis 明确本次判定依据。联合结果是分步骤条件证据，并非将重新选出的种子传入 AA 后完成的一次端到端流程验证。",
             "", "## 统计与证据边界", "", "先对每例样本内种子取平均，再对样本等权。五例仍只有五个独立单位；多个种子不能扩大独立样本数。下游复用配对样本统计，上游使用有界样本均值的 Hoeffding 保守区间与允许失配上限检验，零观察失配也保留不确定性。各来源/分析类型与上游分配 Bonferroni 误差预算，下游另控制版本对及四项指标。",
             "", "容忍阈值未预设、独立样本不足、未知/失败或局部校准未验证时，联合接受结论为证据不足；P>0.05不表示兼容。单样本表只描述观察相似度，不给等效判断。全阴性区域不计为结构一致证据，有信号结果单列。",
             "", "[逐样本相似度](sample_similarity.tsv) · [下游汇总](downstream_summary.tsv) · [上游逐种子](upstream_seeds.tsv) · [上游逐样本](upstream_samples.tsv) · [上游统计](upstream_statistics.tsv) · [联合兼容表](compatibility_combined.tsv)", "",
             "![版本相似性天梯图](similarity_ladder.png)", "", "![限定范围兼容表](compatibility_table.png)", "", "## 分层详细报告", ""]
    lines += ["逐样本天梯图见 [完整图索引](sample_plot_index.tsv)：每个样本×分析类型一张，标记输入来源与固定参考；灰色不可评价项不补为零。单样本只显示观察相似度，不给置信区间或等效判断。CLASS* 表示类别支持不同或未知，相关标签与类别区段指标灰显；表中保留 full_label_comparable、共同支持类别及限定类别范围的观察值，不能把它们解释为完整标签100%一致。", ""]
    lines += [f"- {r['sample_id']} / {r['analysis_kind']}：[PNG](<{r['png']}>) · [PDF](<{r['pdf']}>) · [SVG](<{r['svg']}>)" for r in sample_plots[:20]]
    if config["mode"] == "with_bam":
        lines += ["", "![已知候选种子保留率](upstream_retention_ladder.png)", "", "上游图只显示候选保留率及样本等权区间；CAL?/CAL!/CTRL/NATIVE? 分别标明校准未知、校准失败、受控非历史配置、配置未知。来源版本的局部重评估是实测结果，不是恒等锚点；较高保留率不等于兼容。", ""]
    lines += [f"- [{r['analysis_kind']} / 输入来源 {r['source_version']}](<{r['path']}>)" for r in links]
    if plot_error: lines += ["", "绘图依赖缺失；表格已保存：" + plot_error]
    (out / "report.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    return summary
