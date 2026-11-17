import csv
import math
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from ascompare.concordance_cohort_plot import (UPSTREAM_METRIC, cohort_metric_rows,
                                              plot_cohort_compatibility)
from ascompare.statistics import AGREEMENT_METRICS, analyze


def record(sample, seed, version, detected=True):
    intervals = [["chr1", 10, 20]] if detected else []
    return dict(sample_id=sample, seed_id=seed, version=version,
                status="ok" if detected else "not_detected", labels=["ecDNA"] if detected else [],
                intervals=intervals, features=[dict(label="ecDNA", intervals=intervals)] if detected else [],
                supported_labels=["ecDNA", "BFB"])


def inputs(n=30, primary=True, positive=True, margins=True):
    records = [record(str(i), "seed", v, positive) for i in range(n) for v in ("old", "ref")]
    options = dict(min_samples=2, bootstrap=100, permutations=100,
                   primary_pairs="reference" if primary else [],
                   margins={metric: .8 for metric in AGREEMENT_METRICS} if margins else {})
    _, stats, compatible = analyze(records, "ref", options)
    down = dict(reference_version="ref", source_version="ref", target_version="old", analysis_kind="aa_ac",
                primary_pair=primary, identity_anchor=False, decision=compatible[0]["decision"],
                n_selected_samples=n, eligible_units=n, excluded_units=0,
                signal_regions=n if positive else 0, independence_unit="sample",
                full_label_comparable=True, common_supported_labels=["ecDNA", "BFB"],
                label_scopes=[["ecDNA", "BFB"]],
                observed_statuses=["ok" if positive else "not_detected"],
                decision_metrics=[row for row in stats if row["subset"] == "signal_positive" and row["metric"] in AGREEMENT_METRICS])
    combined = dict(down, comparison_scope="provided_downstream_aa_ac_results_only", reasons=compatible[0]["reasons"])
    return down, combined, stats


class CohortMetricTests(unittest.TestCase):
    def test_plots_exact_adjusted_decision_limits_and_margins(self):
        down, final, stats = inputs()
        rows = cohort_metric_rows([down], [final])
        self.assertEqual(len(rows), 4)
        for row in rows:
            original = next(s for s in stats if s["metric"] == row["metric"] and s["subset"] == "signal_positive")
            self.assertEqual(row["interval_low"], original["agreement_lower"])
            self.assertEqual(row["interval_high"], original["agreement_upper"])
            self.assertEqual(row["interval_kind"], "adjusted_one_sided_decision_limits")
            self.assertEqual(row["agreement_alpha"], .05 / 4)
            self.assertEqual(row["margin"], .8)
            self.assertEqual(row["final_decision"], "within_margin")
            self.assertNotEqual(row["interval_high"], original["conservative_ci_high"])

    def test_added_lower_bound_matches_existing_formula_and_patient_count(self):
        records = [record(sample, str(seed), version, detected=(version == "ref" or sample == "C"))
                   for sample, seeds in (("A", 8), ("B", 1), ("C", 1))
                   for seed in range(seeds) for version in ("old", "ref")]
        _, stats, _ = analyze(records, "ref", dict(min_samples=2, bootstrap=100, permutations=100,
            patient_ids={"A": "P1", "B": "P1", "C": "P2"}, margins={m: .8 for m in AGREEMENT_METRICS}))
        row = next(s for s in stats if s["metric"] == "detection_discordance" and s["subset"] == "signal_positive")
        self.assertEqual((row["n_samples"], row["n_independent_units"], row["n_regions"]), (3, 2, 10))
        self.assertEqual(row["estimate"], .5)
        radius = math.sqrt(math.log(1 / row["agreement_alpha"]) / (2 * 2))
        self.assertEqual(row["agreement_lower"], max(0., .5 - radius))
        self.assertEqual(row["agreement_upper"], min(1., .5 + radius))
        ordinary = next(s for s in stats if s["metric"] == "detection_discordance" and s["subset"] == "all_evaluable")
        self.assertIsNone(ordinary["agreement_lower"])

    def test_exploratory_is_descriptive_even_when_results_identical(self):
        down, final, _ = inputs(primary=False)
        final["decision"] = "within_margin"  # Defence against an inconsistent caller.
        rows = cohort_metric_rows([down], [final])
        for row in rows:
            self.assertEqual(row["interval_kind"], "descriptive_two_sided_ci")
            self.assertEqual(row["status"], "exploratory")
            self.assertEqual(row["final_decision"], "insufficient_evidence")
            self.assertFalse(row["primary_pair"])

    def test_no_margins_preserves_estimate_and_insufficient_conclusion(self):
        down, final, _ = inputs(margins=False)
        for row in cohort_metric_rows([down], [final]):
            self.assertEqual(row["estimate"], 0.)
            self.assertIsNone(row["margin"])
            self.assertEqual(row["status"], "margin_not_defined")
            self.assertEqual(row["final_decision"], "insufficient_evidence")

    def test_all_negative_and_unavailable_do_not_become_zero_loss(self):
        down, final, _ = inputs(positive=False)
        for row in cohort_metric_rows([down], [final]):
            self.assertIsNone(row["estimate"])
            self.assertEqual(row["n_independent_units"], 0)
            self.assertEqual(row["status"], "no_amplification_signal")
            self.assertEqual(row["final_decision"], "insufficient_evidence")
        down.update(eligible_units=0, excluded_units=30, observed_statuses=["missing", "unsupported", "failed"])
        down["decision_metrics"] = []
        row = cohort_metric_rows([down], [final])[0]
        self.assertEqual(row["status"], "unavailable")
        self.assertEqual(row["observed_statuses"], ["missing", "unsupported", "failed"])

    def test_identity_does_not_supply_compatibility_evidence(self):
        down, final, _ = inputs()
        down["identity_anchor"] = True
        row = cohort_metric_rows([down], [final])[0]
        self.assertIsNone(row["estimate"])
        self.assertEqual(row["interval_kind"], "not_assessed")
        self.assertEqual(row["status"], "identity_anchor")
        self.assertEqual(row["final_decision"], "insufficient_evidence")

    def test_local_upstream_uses_actual_bounds_and_keeps_combined_gates(self):
        down, final, _ = inputs()
        final.update(comparison_scope="conditional_local_upstream_and_fixed_seed_aa_ac",
                     decision="insufficient_evidence", upstream_calibration_passed=False,
                     reasons=["local_upstream_calibration_not_passed"])
        up = dict(source_version="ref", target_version="old", estimate=.01,
                  agreement_lower=0., agreement_upper=.17, agreement_alpha=.008,
                  conservative_ci_low=0., conservative_ci_high=.12, confidence_level=.95,
                  upstream_loss_margin=.3, primary_pair=True, decision="within_margin",
                  n_samples=30, n_independent_units=24, independence_unit="patient", min_samples=20,
                  expected_seeds=90, evaluable_seeds=90)
        rows = cohort_metric_rows([down], [final], [up])
        row = rows[-1]
        self.assertEqual(row["metric"], UPSTREAM_METRIC)
        self.assertEqual((row["interval_low"], row["interval_high"], row["margin"]), (0., .17, .3))
        self.assertEqual((row["n_independent_units"], row["n_evaluable_units"]), (24, 90))
        self.assertEqual(row["metric_decision"], "within_margin")
        self.assertEqual(row["final_decision"], "insufficient_evidence")
        self.assertFalse(row["upstream_calibration_passed"])

    def test_shared_class_scope_is_annotated_without_discarding_values(self):
        down, final, _ = inputs()
        down.update(full_label_comparable=False, common_supported_labels=["ecDNA"],
                    label_scopes=[["ecDNA"]])
        final.update(decision="insufficient_evidence", reasons=["category_support_differs_or_unknown"])
        rows = cohort_metric_rows([down], [final])
        for row in rows:
            limited = row["metric"] in {"label_discordance", "class_footprint_distance"}
            self.assertEqual(row["class_scope_limited"], limited)
            self.assertEqual(row["estimate"], 0.)
            self.assertEqual(row["metric_decision"], "within_margin")
            self.assertEqual(row["final_decision"], "insufficient_evidence")
            if limited:
                self.assertEqual(row["scope_annotation"], "class support differs/unknown; shared classes only")
                self.assertEqual(row["common_supported_labels"], ["ecDNA"])
                self.assertEqual(row["label_scopes"], [["ecDNA"]])
        with tempfile.TemporaryDirectory() as folder:
            plot_cohort_compatibility(folder, rows, "ref", "results_only", True)
            svg = (Path(folder) / "cohort_compatibility.svg").read_text()
            self.assertIn("class support differs/unknown; shared classes only", svg)
            image = ET.fromstring(svg)
            for axis_id in ("axes_2", "axes_4"):
                panel = ET.tostring(image.find(f".//*[@id='{axis_id}']"), encoding="unicode")
                self.assertIn("#6b7280", panel)
                self.assertNotIn("#32785a", panel)
            detection = ET.tostring(image.find(".//*[@id='axes_1']"), encoding="unicode")
            self.assertIn("#32785a", detection)

    def test_upstream_unavailable_uses_upstream_statuses_only(self):
        down, final, _ = inputs()
        down["observed_statuses"] = ["failed", "missing"]
        final["comparison_scope"] = "conditional_local_upstream_and_fixed_seed_aa_ac"
        up = dict(source_version="ref", target_version="old", estimate=None,
                  primary_pair=True, observed_statuses=["unsupported"], n_independent_units=0)
        row = cohort_metric_rows([down], [final], [up])[-1]
        self.assertEqual(row["observed_statuses"], ["unsupported"])
        self.assertEqual(row["status"], "unavailable")
        row = cohort_metric_rows([down], [final], [])[-1]
        self.assertEqual(row["observed_statuses"], [])

    def test_sources_and_ac_only_remain_separate_plot_strata(self):
        down, final, _ = inputs()
        ac_down = dict(down, analysis_kind="ac_only", source_version="older")
        ac_final = dict(final, analysis_kind="ac_only", source_version="older",
                        comparison_scope="ac_only_on_fixed_source_aa")
        rows = cohort_metric_rows([down, ac_down], [final, ac_final])
        self.assertEqual(len(rows), 8)
        self.assertEqual({(r["analysis_kind"], r["source_version"]) for r in rows},
                         {("aa_ac", "ref"), ("ac_only", "older")})
        with tempfile.TemporaryDirectory() as folder:
            result = plot_cohort_compatibility(folder, rows, "ref", "results_only", True)
            self.assertEqual(len(result), 1)
            for suffix in ("png", "pdf", "svg"):
                self.assertGreater((Path(folder) / f"cohort_compatibility.{suffix}").stat().st_size, 100)
            with (Path(folder) / "cohort_compatibility_metrics.tsv").open() as handle:
                saved = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(saved), 8)
            self.assertEqual(saved[0]["interval_high"], str(rows[0]["interval_high"]))
            self.assertEqual(saved[0]["synthetic"], "True")


if __name__ == "__main__":
    unittest.main()
