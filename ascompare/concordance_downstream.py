"""Descriptive AA and AC similarity with explicit denominators and provenance.

These observations supplement the existing paired compatibility endpoints. AA
reconstruction and AC feature intervals measure different output stages. A raw
AA overlap here does not introduce a new compatibility acceptance criterion.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

from .reporting import write_tsv
from .statistics import VALID, _independent_means, _pair, interval_jaccard, merge_intervals


METRICS = (
    "detection_agreement", "full_label_agreement", "common_label_agreement",
    "ac_feature_interval_jaccard", "ac_class_footprint_jaccard", "raw_aa_interval_jaccard",
)


def _mean(values):
    return sum(values) / len(values) if values else None


def _raw_aa_eligibility(record):
    if record.get("stage") == "comparison_input":
        return "same_input_comparison_not_eligible"
    if record.get("comparison_eligible") is True and record.get("status") in VALID:
        return None
    # Complete AA may survive an AC-specific parsing failure. Recover its
    # overlap only when the version and input provenance are explicit.
    components = record.get("component_versions", {})
    traceable = (bool(record.get("run_dir")) and bool(components.get("aa"))
                 and bool(components.get("ac")))
    provenance = record.get("provenance_kind")
    imported = (provenance == "imported_declared"
                and record.get("declared_components_validated") is True
                and bool(record.get("source_graph_identity")))
    bounded = provenance == "bounded_reconstruction" and isinstance(record.get("seed"), dict)
    issues = record.get("issues", [])
    if (record.get("status") in {"failed", "missing", "invalid", "unsupported"}
            and record.get("stage") == "ac" and issues and all(issue.get("stage") == "ac" for issue in issues)
            and traceable and (imported or bounded)):
        return None
    return "same_input_comparison_not_eligible"


def _raw_aa(record):
    """Use normalized complete-AA evidence without reopening any source file."""
    unavailable = dict(status="unavailable", intervals=None, reason="complete_AA_evidence_unavailable")
    reason = _raw_aa_eligibility(record)
    if reason:
        return dict(unavailable, reason=reason)
    count = record.get("aa_amplicon_count")
    amplicons = record.get("amplicons")
    if (isinstance(count, bool) or not isinstance(count, int) or count < 0
            or not isinstance(amplicons, list) or "reconstruction_intervals" not in record):
        return unavailable
    if record.get("stage") in {"aa", "input_check", "seed_processing"} and record.get("status") not in VALID:
        return dict(unavailable, reason="AA_or_input_stage_incomplete")
    if any(issue.get("stage") in {"aa", "input_check", "seed_processing"}
           for issue in record.get("issues", [])):
        return dict(unavailable, reason="AA_or_input_issue_recorded")
    projected = record.get("projection") == "overlapping_AA_reconstruction_complete_amplicons_retained"
    if record.get("seed") is not None and not projected:
        return dict(unavailable, reason="AA_not_projected_to_selected_seed")
    if len(amplicons) > count or (not projected and len(amplicons) != count):
        return dict(unavailable, reason="AA_summary_amplicon_count_mismatch")
    if not amplicons:
        if record.get("coverage_status") != "complete":
            return dict(unavailable, reason="empty_AA_completion_unverified")
        if record.get("reconstruction_intervals"):
            return dict(unavailable, reason="empty_AA_has_reconstruction_intervals")
        return dict(status="empty", intervals=[], reason="complete_AA_has_no_reconstructed_interval_in_this_unit")
    try:
        ids = [amp["amplicon_id"] for amp in amplicons]
        if len(set(ids)) != len(ids):
            return dict(unavailable, reason="duplicate_AA_amplicon_id")
        for amp in amplicons:
            if (not amp.get("reconstruction_intervals")
                    or amp.get("stage") == "aa" and amp.get("status") not in VALID):
                return dict(unavailable, reason="AA_amplicon_incomplete")
        expected = merge_intervals([interval for amp in amplicons for interval in amp["reconstruction_intervals"]])
        intervals = merge_intervals(record["reconstruction_intervals"])
    except (ValueError, TypeError, KeyError):
        return dict(unavailable, reason="malformed_AA_reconstruction_intervals")
    if intervals != expected:
        return dict(unavailable, reason="AA_interval_union_differs_from_complete_amplicons")
    return dict(status="complete", intervals=intervals,
                reason="parsed_complete_AA_amplicons_seed_projection" if projected else "parsed_complete_AA_amplicons")


def _same_raw_aa_input(candidate, comparator):
    if (candidate["sample_id"] != comparator["sample_id"]
            or candidate["seed_id"] != comparator["seed_id"]
            or candidate.get("analysis_kind") != comparator.get("analysis_kind")):
        return False
    if (candidate.get("reference_build") and comparator.get("reference_build")
            and candidate["reference_build"] != comparator["reference_build"]):
        return False
    first, second = candidate.get("seed"), comparator.get("seed")
    if first is None and second is None:
        return True
    if not isinstance(first, dict) or not isinstance(second, dict):
        return False
    first_coordinates = tuple(first.get(key) for key in ("chrom", "start", "end"))
    second_coordinates = tuple(second.get(key) for key in ("chrom", "start", "end"))
    return (isinstance(first_coordinates[0], str) and bool(first_coordinates[0])
            and all(isinstance(value, int) and not isinstance(value, bool) for value in first_coordinates[1:])
            and 0 <= first_coordinates[1] < first_coordinates[2]
            and first_coordinates == second_coordinates
            and first.get("seed_id") == candidate["seed_id"]
            and second.get("seed_id") == comparator["seed_id"])


def _trace(record):
    fields = ("run_dir", "source_record_dir", "source_graph_identity", "component_versions",
              "provenance_kind", "completion_evidence", "comparison_limitations", "projection")
    trace = {field: record[field] for field in fields if field in record}
    trace["evidence"] = list(record.get("evidence", []))
    trace["input_artifacts"] = record.get("input_artifacts", [])
    # Keep only recorded AA hashes. Full job inventories remain in the saved
    # execution records and are not repeated for every projected seed.
    trace["aa_artifact_sha256"] = record.get("aa_artifact_sha256", {})
    return trace


def downstream_unit_similarity(candidate, comparator, reference_version, source_version, roles):
    """Compare one fixed input unit, retaining full and shared label scopes."""
    paired = _pair(candidate, comparator, reference_version, "descriptive_stage_similarity")
    signal = paired["evaluable"] and paired["signal_positive"] and paired["is_seed_region"]
    scope = paired["label_scope"] if isinstance(paired["label_scope"], list) else None
    candidate_labels = set(paired["labels"])
    comparator_labels = set(paired["comparator_labels"])
    shared_signal = bool(scope and (candidate_labels | comparator_labels) & set(scope))
    full = bool(paired["full_label_comparable"] and scope)
    full_agreement = float(candidate_labels == comparator_labels) if signal and full else None
    common_agreement = 1 - paired["label_discordance"] if signal and shared_signal and paired["label_discordance"] is not None else None
    if not signal:
        label_reason = "no_signal_positive_evaluable_pair"
    elif scope is None:
        label_reason = "category_support_unknown"
    elif not shared_signal:
        label_reason = "no_common_supported_positive_labels"
    elif not full:
        label_reason = "full_category_support_differs_shared_classes_only"
    else:
        label_reason = "full_category_support_matches"
    aa, aa_comparator = _raw_aa(candidate), _raw_aa(comparator)
    aa_jaccard = None
    aa_reason = "complete_AA_evidence_unavailable"
    if candidate.get("analysis_kind") == "ac_only":
        aa_reason = "not_applicable_fixed_source_AA_not_rerun"
    elif not _same_raw_aa_input(candidate, comparator):
        aa_reason = "raw_AA_comparison_input_mismatch"
    elif paired["is_seed_region"] and all(row["status"] in {"complete", "empty"} for row in (aa, aa_comparator)):
        aa_jaccard = interval_jaccard(aa["intervals"], aa_comparator["intervals"])
        aa_reason = "complete_AA_interval_union_overlap" if aa_jaccard is not None else "both_AA_interval_sets_empty"
    row = dict(
        reference_version=reference_version, source_version=source_version,
        target_version=candidate["version"], comparator_version=comparator["version"],
        comparator_roles=list(roles), analysis_kind=candidate["analysis_kind"],
        sample_id=candidate["sample_id"], seed_id=candidate["seed_id"],
        identity_anchor=candidate["version"] == comparator["version"],
        target_status=paired["status"], comparator_status=paired["comparator_status"],
        target_observed_status=candidate.get("observed_status", paired["status"]),
        comparator_observed_status=comparator.get("observed_status", paired["comparator_status"]),
        target_stage=candidate.get("stage"), comparator_stage=comparator.get("stage"),
        target_reason=candidate.get("reason"), comparator_reason=comparator.get("reason"),
        evaluable=paired["evaluable"] and paired["is_seed_region"], signal_positive=bool(signal),
        full_label_comparable=full, common_supported_labels=scope,
        target_labels=paired["labels"], comparator_labels=paired["comparator_labels"],
        target_supported_labels=candidate.get("supported_labels"),
        comparator_supported_labels=comparator.get("supported_labels"),
        full_label_agreement=full_agreement, common_label_agreement=common_agreement,
        label_comparison_reason=label_reason,
        detection_agreement=paired["detection_agreement"] if signal else None,
        detection_agreement_all_evaluable=paired["detection_agreement"] if paired["is_seed_region"] else None,
        ac_feature_interval_jaccard=paired["interval_jaccard"] if signal else None,
        ac_class_footprint_jaccard=paired["class_footprint_jaccard"] if signal else None,
        ac_class_footprint_reason=paired["class_footprint_reason"],
        raw_aa_interval_jaccard=aa_jaccard, raw_aa_comparison_reason=aa_reason,
        target_raw_aa_status=aa["status"], comparator_raw_aa_status=aa_comparator["status"],
        target_raw_aa_reason=aa["reason"], comparator_raw_aa_reason=aa_comparator["reason"],
        target_seed=candidate.get("seed"), comparator_seed=comparator.get("seed"),
        target_reconstruction_intervals=aa["intervals"], comparator_reconstruction_intervals=aa_comparator["intervals"],
        target_traceability=_trace(candidate), comparator_traceability=_trace(comparator),
        interpretation="descriptive_output_similarity; raw_AA_is_separate_from_AC_feature_intervals",
    )
    return row


def _sample_similarity(rows):
    first = rows[0]
    keys = ("reference_version", "source_version", "target_version", "comparator_version",
            "comparator_roles", "analysis_kind", "sample_id", "identity_anchor")
    result = {key: first[key] for key in keys}
    valid = [row for row in rows if row["evaluable"]]
    result.update(n_expected_units=sum(row["seed_id"] != "__no_seeds__" for row in rows),
                  n_evaluable_units=len(valid),
                  n_signal_units=sum(row["signal_positive"] for row in rows),
                  n_both_negative_units=sum(not row["signal_positive"] for row in valid),
                  target_statuses=sorted({row["target_status"] for row in rows}),
                  comparator_statuses=sorted({row["comparator_status"] for row in rows}),
                  full_label_comparable=bool(valid) and all(row["full_label_comparable"] for row in valid),
                  common_label_scopes=list(dict.fromkeys(tuple(row["common_supported_labels"] or []) for row in rows)),
                  label_comparison_reasons=sorted({row["label_comparison_reason"] for row in rows}),
                  raw_aa_comparison_reasons=sorted({row["raw_aa_comparison_reason"] for row in rows}),
                  interpretation="per_sample_observed_similarity; no_single_sample_equivalence")
    for metric in (*METRICS, "detection_agreement_all_evaluable"):
        values = [row[metric] for row in rows if row[metric] is not None]
        result[metric] = _mean(values)
        result[metric + "_n_units"] = len(values)
    for prefix in ("full_label", "common_label"):
        metric = prefix + "_agreement"
        values = [row[metric] for row in rows if row[metric] is not None]
        count = sum(value == 1. for value in values)
        result[prefix + "_agreement_count"] = count
        result[prefix + "_compared_count"] = len(values)
        result[prefix + "_agreement_ratio"] = f"{count}/{len(values)}" if values else None
    return result


def _cohort_similarity(rows, patient_ids, alpha):
    first = rows[0]
    keys = ("reference_version", "source_version", "target_version", "comparator_version",
            "comparator_roles", "analysis_kind", "identity_anchor")
    result = {key: first[key] for key in keys}
    result.update(n_selected_samples=len(rows), independence_unit="patient" if patient_ids is not None else "sample",
                  n_expected_units=sum(row["n_expected_units"] for row in rows),
                  n_evaluable_units=sum(row["n_evaluable_units"] for row in rows),
                  n_signal_units=sum(row["n_signal_units"] for row in rows),
                  full_label_comparable=all(row["full_label_comparable"] for row in rows),
                  confidence_level=1-alpha, interval_method="descriptive_two_sided_bounded_mean_Hoeffding",
                  interpretation="patient_equal_weight_mean_of_sample_means; supplemental_descriptive_metrics")
    for metric in METRICS:
        sample_values = {row["sample_id"]: row[metric] for row in rows if row[metric] is not None}
        values = list(_independent_means(sample_values, patient_ids))
        estimate = _mean(values)
        low, high = None, None
        if values and not first["identity_anchor"]:
            radius = math.sqrt(math.log(2 / alpha) / (2 * len(values)))
            low, high = max(0., estimate-radius), min(1., estimate+radius)
        result.update({metric: estimate, metric+"_ci_low": low, metric+"_ci_high": high,
                       metric+"_n_samples": len(sample_values), metric+"_n_independent_units": len(values),
                       metric+"_n_units": sum(row[metric+"_n_units"] for row in rows)})
    for prefix in ("full_label", "common_label"):
        agreed = sum(row[prefix+"_agreement_count"] for row in rows)
        total = sum(row[prefix+"_compared_count"] for row in rows)
        result[prefix+"_agreement_count"] = agreed
        result[prefix+"_compared_count"] = total
        result[prefix+"_agreement_ratio"] = f"{agreed}/{total}" if total else None
    return result


def _percent(value):
    return f"{100*value:.1f}%" if value is not None else "不可评价"


def write_downstream_similarity(out, groups, targets, reference_version, patient_ids=None, synthetic=False, alpha=.05):
    """Write stage-specific observed comparisons against source and fixed reference.

    ``groups`` uses the report's normalized (analysis_kind, source_version)
    records. The returned summaries can be included in summary.json. Detailed
    per-unit evidence stays in its TSV, keeping the JSON summary compact.
    """
    if not 0 < alpha < 1:
        raise ValueError("Descriptive similarity alpha must be between zero and one")
    units, samples, cohort, provenance = [], [], [], {}
    provenance_file = "downstream_similarity_provenance.json"
    for (kind, source), records in groups.items():
        lookup = {(row["sample_id"], row["seed_id"], row["version"]): row for row in records}
        comparators = defaultdict(list)
        comparators[reference_version].append("fixed_reference")
        comparators[source].append("input_source")
        selected = list(dict.fromkeys([source, reference_version, *targets]))
        by_sample = defaultdict(list)
        for sample, seed in sorted({(row["sample_id"], row["seed_id"]) for row in records}):
            for target in selected:
                for comparator, roles in comparators.items():
                    row = downstream_unit_similarity(lookup[sample, seed, target], lookup[sample, seed, comparator],
                                                     reference_version, source, roles)
                    for side in ("target", "comparator"):
                        trace = row[side+"_traceability"]
                        identifier = hashlib.sha256(json.dumps(trace, sort_keys=True, ensure_ascii=False).encode()).hexdigest()
                        provenance[identifier] = trace
                        row[side+"_traceability"] = dict(id=identifier, file=provenance_file)
                    row["synthetic"] = synthetic
                    units.append(row)
                    by_sample[sample, target, comparator].append(row)
        local_samples = [_sample_similarity(rows) for rows in by_sample.values()]
        for row in local_samples:
            row["synthetic"] = synthetic
        samples.extend(local_samples)
        by_comparison = defaultdict(list)
        for row in local_samples:
            by_comparison[row["target_version"], row["comparator_version"]].append(row)
        for rows in by_comparison.values():
            summary = _cohort_similarity(rows, patient_ids, alpha)
            summary["synthetic"] = synthetic
            cohort.append(summary)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    files = dict(units="downstream_similarity_units.tsv", samples="downstream_similarity_samples.tsv",
                 cohort="downstream_similarity_cohort.tsv", report="downstream_similarity.md", provenance=provenance_file)
    for key, rows in (("units", units), ("samples", samples), ("cohort", cohort)):
        write_tsv(out / files[key], rows)
    (out / provenance_file).write_text(json.dumps(provenance, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    lines = ["# 下游 AA 与 AC 相似性", "", f"固定参考 `{reference_version}`。同时对照每例声明的输入来源版本。", ""]
    if synthetic:
        lines.extend(["**人工合成测试数据，用于验证程序。**", ""])
    lines.extend([
        "原始 AA 区段取自完整解析的 graph/cycles 重建区段并集。AC 区段取自分类后的 feature BED 并集，类别区段指标同时检查类别和位置。AC 分类失败时，只有明确完整的 AA 证据才保留原始 AA 区段比较。AC-only 固定来源 AA，不评价 AA 重建版本相似度。", "",
        "完整标签一致率要求两版支持完全相同的类别集合。共同类别一致率仅比较双方明确支持的类别。全阴性单位和共同类别内无信号的单位不抬高标签一致率，缺失、失败、空区段与类别未知保留为不可评价。标签计数的分母仅包含相应可评价单位。", "",
        "整组比例先平均每例内的单位，再在患者内平均样本并按独立患者等权。未提供患者编号时按样本独立处理。标签计数比是原始单位计数，等权均值与该计数比可能不同。区间为描述性双侧 Hoeffding 界，身份锚点不提供区间。本表补充观察相似度，兼容判定继续使用预先指定的端点与阈值。", "",
        "| 样本 | 分析 | 目标 | 对照 | 原始 AA 区段 | AC feature 区段 | AC 类别区段 | 完整标签一致率 | 共同类别一致率 |",
        "|---|---|---|---|---|---|---|---|---|",
    ])
    for row in samples:
        full = _percent(row["full_label_agreement"]) + f" ({row['full_label_agreement_ratio'] or '0/0'})"
        common = _percent(row["common_label_agreement"]) + f" ({row['common_label_agreement_ratio'] or '0/0'})"
        lines.append(f"| {row['sample_id']} | {row['analysis_kind']} | {row['target_version']} | {row['comparator_version']} | "
                     f"{_percent(row['raw_aa_interval_jaccard'])} | {_percent(row['ac_feature_interval_jaccard'])} | "
                     f"{_percent(row['ac_class_footprint_jaccard'])} | {full} | {common} |")
    lines.extend(["", f"[逐单位区段、标签、状态与证据]({files['units']}) · [去重证据目录]({files['provenance']}) · [逐样本精确数值]({files['samples']}) · [整组均值与区间]({files['cohort']})", ""])
    (out / files["report"]).write_text("\n".join(lines), encoding="utf-8")
    return dict(samples=samples, cohort=cohort, files=files,
                report_lines=[f"[下游 AA/AC 相似性及完整与共同类别标签一致率]({files['report']}) · "
                              f"[逐样本数值]({files['samples']}) · [逐单位证据]({files['units']})"])
