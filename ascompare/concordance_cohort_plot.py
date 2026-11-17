"""Cohort plots of the endpoints actually used in compatibility decisions.

The figure is a view of existing statistics, never a second inference engine.
In particular it does not turn a descriptive interval, a source identity, or
all-negative seeds into evidence of compatibility. Source versions and AA+AC
versus fixed-AA AC-only comparisons remain separate strata.
"""
from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path

from .reporting import write_tsv
from .statistics import AGREEMENT_METRICS


METRIC_LABELS = {
    "detection_discordance": "Detection discordance",
    "label_discordance": "Class-set discordance",
    "interval_distance": "Interval distance (1 - Jaccard)",
    "class_footprint_distance": "Class footprint distance",
    "known_candidate_not_retained_fraction": "Known candidate seed loss",
}
UPSTREAM_METRIC = "known_candidate_not_retained_fraction"
DECISION_LABELS = {
    "within_margin": "Within margins",
    "exceeds_margin": "Exceeds margin",
    "insufficient_evidence": "Insufficient evidence",
}
COLORS = {"within_margin": "#32785a", "exceeds_margin": "#ac3c36", "insufficient_evidence": "#6b7280"}


def _number(value):
    if isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _metric_row(down, stat, metric, final, stage):
    stat = stat or {}
    primary = bool(down.get("primary_pair")) and bool(stat.get("primary_pair"))
    identity = bool(down.get("identity_anchor"))
    estimate = _number(stat.get("estimate"))
    lower, upper = (_number(stat.get("agreement_" + side)) for side in ("lower", "upper"))
    alpha = _number(stat.get("agreement_alpha"))
    margin = _number(stat.get("margin" if stage == "downstream" else "upstream_loss_margin"))
    formal_limits = primary and not identity and alpha is not None and lower is not None and upper is not None
    # These are separate one-sided limits from the existing decision, not an
    # ordinary two-sided confidence interval. Exploratory plots use only the
    # descriptive two-sided interval and remain unassessed.
    low = lower if formal_limits else _number(stat.get("conservative_ci_low"))
    high = upper if formal_limits else _number(stat.get("conservative_ci_high"))
    interval_kind = "adjusted_one_sided_decision_limits" if formal_limits else "descriptive_two_sided_ci"
    confidence = 1-alpha if formal_limits else _number(stat.get("confidence_level"))
    n = stat.get("n_independent_units", 0) or 0
    if identity:
        status = "identity_anchor"
        estimate, low, high, confidence, interval_kind = None, None, None, None, "not_assessed"
    elif not n or estimate is None:
        status = ("no_amplification_signal" if stage == "downstream" and down.get("eligible_units", 0)
                  and down.get("signal_regions") == 0 else "unavailable")
    elif not primary:
        status = "exploratory"
    elif margin is None:
        status = "margin_not_defined"
    else:
        status = "estimated"
    metric_decision = stat.get("agreement_decision", stat.get("decision", "not_assessed"))
    class_scope_limited = (not identity and stage == "downstream"
                           and metric in {"label_discordance", "class_footprint_distance"}
                           and down.get("full_label_comparable") is not True)
    annotation = "class support differs/unknown; shared classes only" if class_scope_limited else ""
    # The final decision also includes label/resource completeness, upstream
    # controls, source bridges, and the predefined primary-pair restriction.
    final_decision = final.get("decision", "insufficient_evidence")
    if not primary or identity:
        final_decision = "insufficient_evidence"
    return dict(
        reference_version=down["reference_version"], source_version=down["source_version"],
        target_version=down["target_version"], analysis_kind=down["analysis_kind"],
        comparison_scope=final.get("comparison_scope", down.get("comparison_scope", "")),
        stage=stage, metric=metric, subset="signal_positive" if stage == "downstream" else "known_candidates",
        estimate=estimate, interval_low=low, interval_high=high, interval_kind=interval_kind,
        displayed_confidence_level=confidence, agreement_lower=lower, agreement_upper=upper,
        agreement_alpha=alpha, margin=margin,
        conservative_ci_low=_number(stat.get("conservative_ci_low")),
        conservative_ci_high=_number(stat.get("conservative_ci_high")),
        descriptive_confidence_level=_number(stat.get("confidence_level")),
        n_samples=stat.get("n_samples", 0), n_independent_units=n,
        independence_unit=stat.get("independence_unit", down.get("independence_unit", "sample")),
        n_selected_samples=down.get("n_selected_samples", 0),
        n_evaluable_units=stat.get("n_regions", stat.get("evaluable_seeds", 0)),
        n_expected_units=stat.get("total_regions", stat.get("expected_seeds", 0)),
        n_excluded_units=stat.get("excluded_regions", down.get("excluded_units", 0)),
        min_independent_units=stat.get("min_samples"), primary_pair=primary,
        identity_anchor=identity, metric_decision=metric_decision, final_decision=final_decision,
        status=status, observed_statuses=(stat.get("observed_statuses", []) if stage == "upstream"
                                        else down.get("observed_statuses", [])),
        reasons=final.get("reasons", down.get("reasons", [])),
        component_reasons=stat.get("reasons", []),
        full_label_comparable=down.get("full_label_comparable"),
        class_scope_limited=class_scope_limited, scope_annotation=annotation,
        common_supported_labels=down.get("common_supported_labels"),
        label_scopes=down.get("label_scopes"),
        upstream_calibration_passed=final.get("upstream_calibration_passed"),
        source_bridge_passed=final.get("source_bridge_passed"),
        synthetic=bool(down.get("synthetic", False)),
    )


def cohort_metric_rows(downstream, combined, upstream=None):
    """Return exact plotted endpoints, retaining source/kind and decision gates.

    Each downstream summary must carry ``decision_metrics``, the four oriented
    signal-positive statistics generated by ``statistics.analyze``. No metric
    is inferred from a rounded similarity value. Missing statistics stay blank.
    """
    final_index = {(r["analysis_kind"], r["source_version"], r["target_version"]): r for r in combined}
    up_index = {(r["source_version"], r["target_version"]): r for r in (upstream or [])}
    rows = []
    for down in downstream:
        key = down["analysis_kind"], down["source_version"], down["target_version"]
        final = final_index.get(key, {})
        stats = {r["metric"]: r for r in down.get("decision_metrics", []) if r is not None}
        for metric in AGREEMENT_METRICS:
            rows.append(_metric_row(down, stats.get(metric), metric, final, "downstream"))
        if final.get("comparison_scope") == "conditional_local_upstream_and_fixed_seed_aa_ac":
            rows.append(_metric_row(down, up_index.get(key[1:]), UPSTREAM_METRIC, final, "upstream"))
    return rows


def _missing_text(row):
    if row["status"] == "identity_anchor":
        return "identity; not validation"
    if row["status"] == "no_amplification_signal":
        return "no positive-signal units"
    statuses = set(row.get("observed_statuses", []))
    known = [s for s in ("failed", "missing", "unsupported", "invalid", "interrupted", "blocked", "no_seeds", "not_run") if s in statuses]
    return ", ".join(known) if known else "not evaluable"


def plot_cohort_compatibility(out, rows, ref, mode, synthetic):
    """Write a cohort forest figure and the exact table behind it.

    Lower values mean smaller differences. Each endpoint shows its own
    independent patient/sample count. Empty and all-negative rows are explicit;
    only signal-positive downstream observations enter these endpoints.
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = [dict(row, synthetic=bool(synthetic)) for row in rows]
    write_tsv(out / "cohort_compatibility_metrics.tsv", rows)
    if not rows:
        return []
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    groups = defaultdict(list)
    for row in rows:
        groups[row["analysis_kind"], row["source_version"]].append(row)
    has_upstream = any(r["stage"] == "upstream" for r in rows)
    metrics = list(AGREEMENT_METRICS) + ([UPSTREAM_METRIC] if has_upstream else [])
    targets = {key: list(dict.fromkeys(r["target_version"] for r in local)) for key, local in groups.items()}
    heights = [max(2.4, .39 * len(targets[key]) + .9) for key in groups]
    with plt.rc_context({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.spines.left": False, "axes.linewidth": .6, "savefig.dpi": 180}):
        figure_height = sum(heights) + 2.70
        fig = plt.figure(figsize=(3.0 * len(metrics) + 3.8, figure_height))
        grid = fig.add_gridspec(len(groups), len(metrics) + 1, height_ratios=heights,
                               width_ratios=[1] * len(metrics) + [1.35],
                               left=.10, right=.99, bottom=1.05/figure_height, top=1-1.65/figure_height,
                               wspace=.12, hspace=.7)
        for group_no, ((kind, source), local) in enumerate(groups.items()):
            versions = targets[kind, source]
            index = {(r["target_version"], r["metric"]): r for r in local}
            for col, metric in enumerate(metrics):
                ax = fig.add_subplot(grid[group_no, col])
                ax.set_xlim(-.03, 1.03)
                ax.set_ylim(len(versions) - .4, -.8)
                ax.set_xticks([0, .25, .5, .75, 1])
                ax.set_xticklabels(["0", "25", "50", "75", "100"])
                ax.set_xlabel("Difference (%)", fontsize=8)
                ax.grid(axis="x", color="#ececec", linewidth=.6)
                ax.set_axisbelow(True)
                ax.set_yticks(range(len(versions)))
                ax.set_yticklabels(versions if col == 0 else [""] * len(versions))
                ax.tick_params(axis="y", length=0, pad=6)
                ax.set_title(METRIC_LABELS[metric], fontsize=9, pad=10)
                for y, target in enumerate(versions):
                    row = index.get((target, metric))
                    if row is None:
                        ax.text(.04, y, "not included in this analysis", color="#777777", fontsize=6.4, va="center")
                        continue
                    if row["estimate"] is None:
                        ax.text(.04, y, _missing_text(row), color="#777777", fontsize=6.4, va="center")
                        if row.get("scope_annotation"):
                            ax.text(.04, y+.28, row["scope_annotation"], color="#777777", fontsize=5.5, va="center")
                        continue
                    formal = row["interval_kind"] == "adjusted_one_sided_decision_limits"
                    color = COLORS.get(row["metric_decision"], COLORS["insufficient_evidence"])
                    if not formal or row.get("class_scope_limited"):
                        color = COLORS["insufficient_evidence"]
                    if row["interval_low"] is not None and row["interval_high"] is not None:
                        ax.plot([row["interval_low"], row["interval_high"]], [y, y],
                                color=color, linewidth=1.3, linestyle="-" if formal else ":", zorder=3)
                    ax.plot(row["estimate"], y, marker="o" if formal else "o", markersize=4,
                            color=color, markerfacecolor=color if formal else "white", zorder=4)
                    if row["margin"] is not None:
                        ax.plot([row["margin"], row["margin"]], [y-.21, y+.21], color="#222222", linewidth=1.2, zorder=2)
                    suffix = "" if row["status"] == "estimated" else " | margin unset" if row["status"] == "margin_not_defined" else " | exploratory"
                    counts = f"n={row['n_independent_units']}; units={row['n_evaluable_units']}/{row['n_expected_units']}"
                    ax.text(1.0, y+.20 if row.get("scope_annotation") else y+.25, counts + suffix,
                            color="#555555", fontsize=5.9, ha="right", va="center")
                    if row.get("scope_annotation"):
                        ax.text(1.0, y+.43, row["scope_annotation"], color="#777777", fontsize=5.5, ha="right", va="center")
                if col == 0:
                    label = "AA+AC" if kind == "aa_ac" else "AC only (fixed source AA)"
                    ax.annotate(f"{label} | input source: {source}", xy=(0, 1), xycoords="axes fraction",
                                xytext=(0, 32), textcoords="offset points", fontsize=9,
                                fontweight="bold", va="bottom")
            decision_ax = fig.add_subplot(grid[group_no, -1])
            decision_ax.set_ylim(len(versions) - .4, -.8)
            decision_ax.set_xlim(0, 1)
            decision_ax.axis("off")
            decision_ax.set_title("Combined cohort conclusion", fontsize=9, pad=10)
            for y, target in enumerate(versions):
                row = index[target, AGREEMENT_METRICS[0]]
                conclusion = "Identity anchor" if row["identity_anchor"] else DECISION_LABELS.get(row["final_decision"], "Insufficient evidence")
                decision_ax.text(.02, y-.05, conclusion, color=COLORS.get(row["final_decision"], "#6b7280"), fontsize=8, va="center")
                selected = row["n_selected_samples"]
                unit = "patients" if row["independence_unit"] == "patient" else "samples"
                label = f"{selected} input samples; n = independent {unit}"
                decision_ax.text(.02, y+.22, label, color="#777777", fontsize=5.8, va="center")
        title = f"Cohort cross-version compatibility | fixed reference: {ref}"
        if synthetic:
            title += "\nSYNTHETIC TEST DATA - no patient or biological conclusion"
        fig.suptitle(title, y=1-.10/figure_height, fontsize=12)
        legend = [Line2D([0], [0], color="#6b7280", marker="o", markersize=4, label="Adjusted one-sided decision limits"),
                  Line2D([0], [0], color="#6b7280", marker="o", markerfacecolor="white", markersize=4, linestyle=":", label="Descriptive two-sided CI (exploratory)"),
                  Line2D([0], [0], color="#222222", marker="|", markersize=9, linestyle="none", label="Predefined maximum difference")]
        fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(.5, 1-.75/figure_height), ncol=3, frameon=False, fontsize=7)
        fig.text(.5, .10/figure_height,
                 "Downstream endpoints use positive-signal units; units = metric-evaluable / expected. Seeds are averaged within samples; linked samples within patients.\n"
                 "Each limit is one-sided at 1 - endpoint alpha (see plotted-data TSV); it is not an ordinary two-sided CI.\n"
                 "Final decisions retain completeness, source-control and bridge requirements. Missing margins or inadequate data remain insufficient evidence.",
                 ha="center", va="bottom", fontsize=6.8)
        paths = []
        for suffix in ("png", "pdf", "svg"):
            path = out / f"cohort_compatibility.{suffix}"
            fig.savefig(path, bbox_inches="tight")
            paths.append(str(path.name))
        plt.close(fig)
    return [dict(figure="cohort_compatibility", reference_version=ref, mode=mode,
                 synthetic=bool(synthetic), metrics_table="cohort_compatibility_metrics.tsv", files=paths)]
