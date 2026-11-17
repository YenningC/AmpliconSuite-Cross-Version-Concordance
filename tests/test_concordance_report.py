import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ascompare.concordance_report import _upstream, write_concordance_reports
from ascompare.statistics import AGREEMENT_METRICS


REF = "suite-1.6.2"
SOURCE = "suite-1.4.0"
OTHER = "suite-1.6.0"


def fixture(n=3, kind="aa_ac", mode="with_bam"):
    samples = [dict(sample_id=f"SYNTHETIC_{i}", source_version=SOURCE, upstream_seed_ids=["a"], source_bridge_status="passed") for i in range(n)]
    config = dict(reference_version=REF, target_versions=[SOURCE, OTHER, REF], samples=samples,
                  mode=mode, synthetic=True, statistics=dict(min_samples=2, bootstrap=100, permutations=100,
                      margins={}, upstream_loss_margin=None, seed=91))
    downstream, upstream = [], []
    for sample in samples:
        for version in config["target_versions"]:
            downstream.append(dict(sample_id=sample["sample_id"], version=version, seed_id="a", analysis_kind=kind,
                source_version=SOURCE, source_aa_version="AA_TEST_SOURCE", comparison_eligible=True,
                source_graph_identity="SYNTHETIC_GRAPH_"+sample["sample_id"],
                status="ok", labels=["ecDNA"], supported_labels=["ecDNA", "BFB"], intervals=[["chr1", 10, 20]],
                features=[dict(label="ecDNA", intervals=[["chr1", 10, 20]])]))
            upstream.append(dict(sample_id=sample["sample_id"], version=version, seed_id="a", source_version=SOURCE,
                cnvkit_version="SYNTHETIC_OLD" if version == SOURCE else "SYNTHETIC_NEW",
                status="retained", comparison_eligible=True,
                native_recipe_match=True,
                seed=dict(chrom='chr1', start=10, end=20), reference_mode='matched_normal', min_seed_overlap=1.,
                retention_fraction=1., matching_rule="SYNTHETIC overlap rule, not an acceptance cutoff",
                calibration_role='source_control' if version == SOURCE else 'candidate_observation',
                calibration_status='passed' if version == SOURCE else 'unverified'))
    return config, downstream, upstream


def accepted_options(config):
    # Deliberately broad synthetic-only tolerances test the decision machinery.
    config["statistics"].update(margins={m: .8 for m in AGREEMENT_METRICS}, upstream_loss_margin=.8)


def run(config, downstream, upstream):
    with tempfile.TemporaryDirectory() as folder, patch("ascompare.concordance_report._plots"):
        summary = write_concordance_reports(folder, downstream, upstream, config)
        with (Path(folder)/"sample_similarity.tsv").open() as handle:
            summary["sample_rows_for_test"] = list(csv.DictReader(handle, delimiter="\t"))
        summary["report_for_test"] = (Path(folder)/"report.md").read_text()
        return summary


class ConcordanceReportTests(unittest.TestCase):
    def test_source_to_target_is_not_silently_replaced_by_reference_comparison(self):
        config, downstream, upstream = fixture()
        for row in downstream:
            if row["version"] == SOURCE:
                row.update(labels=["BFB"], features=[dict(label="BFB", intervals=row["intervals"])])
        result = run(config, downstream, upstream)
        compared = next(r for r in result["downstream"] if r["target_version"] == OTHER)
        self.assertEqual(compared["source_version"], SOURCE)
        self.assertEqual(compared["reference_version"], REF)
        self.assertEqual(compared["label_agreement"], 0.)
        self.assertEqual(compared["interval_jaccard"], 1.)
        self.assertEqual(compared["class_footprint_jaccard"], 0.)

    def test_no_margins_or_five_samples_cannot_establish_equivalence(self):
        config, downstream, upstream = fixture(5)
        config["statistics"]["min_samples"] = 20
        result = run(config, downstream, upstream)
        self.assertTrue(all(r["decision"] == "insufficient_evidence" for r in result["compatibility_combined"]))
        self.assertIn("upstream_loss_margin_not_predefined", result["upstream"][-1]["reasons"])

    def test_unverified_calibration_blocks_combined_but_keeps_retention_statistics(self):
        config, downstream, upstream = fixture(30)
        accepted_options(config)
        for row in upstream: row["calibration_status"] = "unverified"
        result = run(config, downstream, upstream)
        combined = next(r for r in result["compatibility_combined"] if r["target_version"] == OTHER)
        self.assertEqual(combined["downstream_decision"], "within_margin")
        self.assertEqual(combined["upstream_decision"], "within_margin")
        self.assertEqual(combined["decision"], "insufficient_evidence")
        self.assertIn("local_upstream_calibration_not_passed", combined["reasons"])
        self.assertEqual(next(r for r in result["upstream"] if r["target_version"] == OTHER)["estimate"], 0.)

    def test_validated_components_allow_only_narrow_conditional_acceptance(self):
        config, downstream, upstream = fixture(30)
        accepted_options(config)
        result = run(config, downstream, upstream)
        combined = next(r for r in result["compatibility_combined"] if r["target_version"] == OTHER)
        self.assertEqual(combined["decision"], "within_margin")
        self.assertEqual(combined["comparison_scope"], "conditional_local_upstream_and_fixed_seed_aa_ac")
        self.assertIn("complete_suite_interchangeability", combined["not_established"])
        self.assertEqual(next(r for r in result["compatibility_combined"] if r["target_version"] == SOURCE)["decision"], "insufficient_evidence")

    def test_missing_original_source_bridge_blocks_combined_acceptance(self):
        config, downstream, upstream = fixture(30)
        accepted_options(config)
        config["samples"][0].pop("source_bridge_status")
        result = run(config, downstream, upstream)
        combined = next(r for r in result["compatibility_combined"] if r["target_version"] == OTHER)
        self.assertEqual(combined["downstream_decision"], "within_margin")
        self.assertEqual(combined["upstream_decision"], "within_margin")
        self.assertEqual(combined["decision"], "insufficient_evidence")
        self.assertIn("input_source_to_rerun_bridge_not_passed", combined["reasons"])

    def test_controlled_cnvkit_cannot_establish_native_suite_compatibility(self):
        config, downstream, upstream = fixture(30)
        accepted_options(config)
        for row in upstream: row["native_recipe_match"] = False
        result = run(config, downstream, upstream)
        combined = next(r for r in result["compatibility_combined"] if r["target_version"] == OTHER)
        self.assertEqual(combined["upstream_decision"], "within_margin")
        self.assertEqual(combined["decision"], "insufficient_evidence")
        self.assertIn("controlled_CNVkit_not_historical_recipe", combined["reasons"])

    def test_unvalidated_upstream_exceedance_is_not_historical_incompatibility(self):
        for change in ({"calibration_status": "unverified"}, {"calibration_status": "failed"},
                       {"native_recipe_match": False}, {"native_recipe_match": None}):
            config, downstream, upstream = fixture(30)
            accepted_options(config)
            config["statistics"]["upstream_loss_margin"] = .2
            for row in upstream:
                row.update(status="not_retained", retention_fraction=0., **change)
            result = run(config, downstream, upstream)
            combined = next(r for r in result["compatibility_combined"] if r["target_version"] == OTHER)
            self.assertEqual(combined["upstream_decision"], "exceeds_margin")
            self.assertEqual(combined["downstream_decision"], "within_margin")
            self.assertEqual(combined["decision"], "insufficient_evidence")
            self.assertEqual(combined["decision_basis"], [])
            self.assertIn("upstream_exceedance_not_validated_by_independent_source_control", combined["reasons"])

    def test_validated_upstream_exceedance_can_support_conditional_rejection(self):
        config, downstream, upstream = fixture(30)
        accepted_options(config)
        config["statistics"]["upstream_loss_margin"] = .2
        for row in upstream:
            if row['version'] == OTHER:
                row.update(status="not_retained", retention_fraction=0.)
        result = run(config, downstream, upstream)
        combined = next(r for r in result["compatibility_combined"] if r["target_version"] == OTHER)
        self.assertEqual(combined["decision"], "exceeds_margin")
        self.assertEqual(combined["decision_basis"], ["local_upstream_exceeds_margin_after_source_control"])
        summary = next(r for r in result['upstream'] if r['target_version'] == OTHER)
        self.assertEqual(summary['retention_difference'], -1.)
        self.assertEqual(summary['retention_difference_n_seeds'], 30)
        self.assertEqual(summary['source_control_cnvkit_versions'], ['SYNTHETIC_OLD'])

    def test_nonzero_losses_can_pass_after_independent_control(self):
        config, downstream, upstream = fixture(30)
        accepted_options(config)
        up = next(r for r in upstream if r['version'] == OTHER)
        up.update(status='not_retained', retention_fraction=0.)
        result = run(config, downstream, upstream)
        compared = next(r for r in result['compatibility_combined'] if r['target_version'] == OTHER)
        summary = next(r for r in result['upstream'] if r['target_version'] == OTHER)
        self.assertEqual(compared['decision'], 'within_margin')
        self.assertAlmostEqual(summary['retention_difference'], -1/30)

    def test_source_control_requires_native_recipe_same_scope_and_complete_roles(self):
        for edit in ('old_rows', 'coordinates', 'threshold', 'reference_mode', 'missing_control', 'mixed_source_cnvkit'):
            with self.subTest(edit=edit):
                config, downstream, upstream = fixture(30)
                accepted_options(config)
                if edit == 'old_rows':
                    upstream[0].pop('calibration_role')
                elif edit == 'coordinates':
                    upstream[1]['seed'] = dict(chrom='chr2', start=10, end=20)
                elif edit == 'threshold':
                    upstream[1]['min_seed_overlap'] = .9
                elif edit == 'reference_mode':
                    upstream[1]['reference_mode'] = 'flat'
                elif edit == 'mixed_source_cnvkit':
                    upstream[0]['cnvkit_version'] = 'SYNTHETIC_OTHER'
                else:
                    upstream.pop(0)
                result = run(config, downstream, upstream)
                compared = next(r for r in result['compatibility_combined'] if r['target_version'] == OTHER)
                self.assertEqual(compared['decision'], 'insufficient_evidence')

    def test_first_page_preserves_directed_delta_and_event_denominators(self):
        config, downstream, upstream = fixture(4)
        for row in downstream:
            if row['version'] == REF and row['sample_id'] in ('SYNTHETIC_0', 'SYNTHETIC_1'):
                row.update(status='not_detected', labels=[], intervals=[], features=[])
            elif row['version'] == REF and row['sample_id'] == 'SYNTHETIC_2':
                row.update(labels=['BFB'], features=[dict(label='BFB', intervals=row['intervals'])])
            elif row['version'] == REF and row['sample_id'] == 'SYNTHETIC_3':
                row.update(status='missing')
        result = run(config, downstream, upstream)
        # Stored reference stats use SOURCE minus REF. Report uses REF minus SOURCE.
        compared = next(r for r in result['downstream'] if r['target_version'] == REF)
        self.assertAlmostEqual(compared['detection_difference'], -2/3)
        self.assertEqual(compared['detection_difference_ci_low'], -1.)
        self.assertEqual(compared['detection_difference_ci_high'], 1.)
        self.assertEqual((compared['extra_detection_count'], compared['missed_detection_count']), (0, 2))
        self.assertEqual((compared['class_change_count'], compared['class_change_n_regions']), (1, 1))
        self.assertEqual((compared['detection_difference_n_regions'], compared['expected_regions']), (3, 4))
        self.assertIn('-66.7', result['report_for_test'])

    def test_patient_cluster_counts_and_exploratory_mode(self):
        config, downstream, upstream = fixture(3)
        for sid, sample in enumerate(config['samples']):
            sample['patient_id'] = 'P1' if sid < 2 else 'P2'
        for row in upstream:
            if row['version'] == OTHER and row['sample_id'] != 'SYNTHETIC_2':
                row.update(status='not_retained', retention_fraction=0.)
        config['statistics'].update(primary_pairs=[])
        accepted_options(config)
        result = run(config, downstream, upstream)
        row = next(r for r in result['upstream'] if r['target_version'] == OTHER)
        self.assertEqual(result['independence_unit'], 'patient')
        self.assertEqual(result['selected_independent_units'], 2)
        self.assertEqual((row['n_samples'], row['n_independent_units']), (3, 2))
        self.assertEqual(row['retention_difference'], -.5)
        self.assertIsNone(row['agreement_alpha'])
        self.assertTrue(all(r['decision'] == 'insufficient_evidence' for r in result['compatibility_combined']))
        self.assertIn('探索性分析', result['report_for_test'])

    def test_reliable_downstream_exceedance_survives_unknown_upstream_but_not_failed_source_bridge(self):
        for bridge in ("passed", "failed", "unverified"):
            config, downstream, upstream = fixture(30)
            accepted_options(config)
            config["statistics"]["margins"] = {m: .2 for m in AGREEMENT_METRICS}
            for row in downstream:
                if row["version"] == OTHER:
                    row.update(labels=["BFB"], features=[dict(label="BFB", intervals=row["intervals"])])
            for row in upstream: row.update(calibration_status="unverified", native_recipe_match=False)
            config["samples"][0]["source_bridge_status"] = bridge
            result = run(config, downstream, upstream)
            combined = next(r for r in result["compatibility_combined"] if r["target_version"] == OTHER)
            self.assertEqual(combined["downstream_decision"], "exceeds_margin")
            self.assertEqual(combined["decision"], "exceeds_margin" if bridge == "passed" else "insufficient_evidence")
            self.assertEqual(combined["decision_basis"], ["downstream_metrics_exceed_margin"] if bridge == "passed" else [])
            if bridge != "passed": self.assertIn("input_source_to_rerun_bridge_not_passed", combined["reasons"])

    def test_ac_only_and_results_only_never_claim_aa_or_upstream_compatibility(self):
        config, downstream, upstream = fixture(30, "ac_only", "results_only")
        accepted_options(config)
        result = run(config, downstream, [])
        combined = next(r for r in result["compatibility_combined"] if r["target_version"] == OTHER)
        self.assertEqual(combined["decision"], "within_margin")
        self.assertEqual(combined["upstream_decision"], "not_assessed")
        self.assertEqual(combined["comparison_scope"], "ac_only_on_fixed_source_aa")
        self.assertIn("aa_reconstruction_compatibility", combined["not_established"])

    def test_missing_or_different_ac_source_blocks_fixed_input_inference(self):
        config, downstream, upstream = fixture(30, "ac_only", "results_only")
        accepted_options(config)
        downstream[0]["source_aa_version"] = "OTHER_SOURCE"
        result = run(config, downstream, [])
        self.assertTrue(all(r["decision"] == "insufficient_evidence" for r in result["compatibility_combined"]))

    def test_same_aa_version_with_different_or_unknown_graph_is_not_fixed_input(self):
        for graph_id in ("ANOTHER_SYNTHETIC_GRAPH", None):
            config, downstream, upstream = fixture(30, "ac_only", "results_only")
            accepted_options(config)
            downstream[0]["source_graph_identity"] = graph_id
            result = run(config, downstream, [])
            self.assertTrue(all(r["decision"] == "insufficient_evidence" for r in result["compatibility_combined"]))

    def test_eligibility_requires_explicit_boolean_true(self):
        config, downstream, upstream = fixture(30)
        accepted_options(config)
        downstream[0]["comparison_eligible"] = "true"
        upstream[1].pop("comparison_eligible")
        result = run(config, downstream, upstream)
        row = next(r for r in result["compatibility_combined"] if r["target_version"] == OTHER)
        self.assertEqual(row["decision"], "insufficient_evidence")
        self.assertIn("unevaluable_regions_present", row["reasons"])

    def test_sample_equal_upstream_weight_and_missing_are_not_zero_loss(self):
        samples = {"A": dict(source_version=SOURCE, upstream_seed_ids=[str(i) for i in range(100)]),
                   "B": dict(source_version=SOURCE, upstream_seed_ids=["one"]),
                   "C": dict(source_version=SOURCE, upstream_seed_ids=["missing"])}
        records = [dict(sample_id="A", version=OTHER, source_version=SOURCE, seed_id=str(i), status="retained",
                    retention_fraction=1., matching_rule="rule", cnvkit_version="TEST", comparison_eligible=True, calibration_status="passed") for i in range(100)]
        records.append(dict(sample_id="B", version=OTHER, source_version=SOURCE, seed_id="one", status="not_retained",
                    retention_fraction=0., matching_rule="rule", cnvkit_version="TEST", comparison_eligible=True, calibration_status="passed"))
        seeds, cases, stats = _upstream(records, samples, [OTHER], REF, .05, .025, 2, .8)
        self.assertEqual(stats[0]["n_samples"], 2)
        self.assertEqual(stats[0]["estimate"], .5)
        self.assertEqual(stats[0]["evaluable_seeds"], 101)
        self.assertEqual(next(r for r in seeds if r["sample_id"] == "C")["status"], "missing")
        self.assertIsNone(next(r for r in cases if r["sample_id"] == "C")["seed_loss"])
        self.assertEqual(stats[0]["decision"], "insufficient_evidence")

    def test_no_seeds_and_undeclared_universe_do_not_count_as_successes(self):
        config, downstream, upstream = fixture(30)
        accepted_options(config)
        config["samples"][0]["upstream_seed_ids"] = []
        upstream = [r for r in upstream if r["sample_id"] != config["samples"][0]["sample_id"]]
        config["samples"][1].pop("upstream_seed_ids")
        result = run(config, downstream, upstream)
        up = next(r for r in result["upstream"] if r["target_version"] == OTHER)
        self.assertEqual(up["n_samples"], 29)
        self.assertEqual(up["no_seed_samples"], 1)
        self.assertIn("candidate_universe_not_declared", up["reasons"])
        self.assertEqual(up["decision"], "insufficient_evidence")

    def test_missing_entire_input_sample_remains_in_comparison_scope(self):
        config, downstream, upstream = fixture(30)
        accepted_options(config)
        downstream = [r for r in downstream if r["sample_id"] != config["samples"][0]["sample_id"]]
        result = run(config, downstream, upstream)
        combined = next(r for r in result["compatibility_combined"] if r["target_version"] == OTHER)
        self.assertEqual(combined["n_selected_samples"], 30)
        self.assertEqual(combined["decision"], "insufficient_evidence")

    def test_seed_missing_from_every_version_still_enters_declared_denominator(self):
        config, downstream, upstream = fixture(30)
        accepted_options(config)
        for sample in config["samples"]: sample["selected_seed_ids"] = ["a", "absent_in_every_version"]
        result = run(config, downstream, upstream)
        compared = next(r for r in result["downstream"] if r["target_version"] == OTHER)
        self.assertEqual(compared["excluded_units"], 30)
        self.assertEqual(compared["decision"], "insufficient_evidence")
        sample = next(r for r in result["sample_rows_for_test"] if r["target_version"] == OTHER)
        self.assertEqual(sample["n_units"], "2")
        self.assertEqual(sample["n_evaluable_units"], "1")
        self.assertEqual(sample["status"], "partial")

    def test_results_only_whole_outputs_do_not_expand_into_unrun_selected_seeds(self):
        config, downstream, upstream = fixture(30, mode="results_only")
        accepted_options(config)
        for sample in config["samples"]: sample["selected_seed_ids"] = ["seed_not_a_whole_result_unit"]
        result = run(config, downstream, [])
        compared = next(r for r in result["downstream"] if r["target_version"] == OTHER)
        self.assertEqual(compared["excluded_units"], 0)
        self.assertEqual(compared["decision"], "within_margin")

    def test_kinds_and_source_versions_are_separate_strata(self):
        config, downstream, upstream = fixture(2)
        config["samples"][1]["source_version"] = OTHER
        for row in downstream:
            if row["sample_id"] == config["samples"][1]["sample_id"]: row["source_version"] = OTHER
        for row in upstream:
            if row["sample_id"] == config["samples"][1]["sample_id"]:
                role = 'source_control' if row['version'] == OTHER else 'candidate_observation'
                row.update(source_version=OTHER, calibration_role=role)
        downstream += [{**r, "analysis_kind": "ac_only"} for r in downstream]
        result = run(config, downstream, upstream)
        self.assertEqual(len(result["compatibility_combined"]), 12)
        self.assertEqual({r["n_selected_samples"] for r in result["compatibility_combined"]}, {1})
        self.assertEqual(set(result["analysis_kinds"]), {"ac_only", "aa_ac"})

    def test_sample_ladder_keeps_source_anchor_even_if_source_is_not_a_target(self):
        config, downstream, upstream = fixture(2)
        config["target_versions"] = [OTHER, REF]
        result = run(config, downstream, upstream)
        self.assertEqual(len(result["compatibility_combined"]), 2)
        anchors = [r for r in result["sample_rows_for_test"] if r["target_version"] == SOURCE]
        self.assertEqual(len(anchors), 2)
        self.assertTrue(all(r["identity_anchor"] == "True" for r in anchors))

    def test_report_and_scoped_figures_are_written_from_synthetic_data(self):
        config, downstream, upstream = fixture(2)
        downstream[1]["status"] = "missing"
        downstream[2]["supported_labels"] = ["ecDNA"]
        upstream[0].update(calibration_status="unverified", native_recipe_match=False)
        with tempfile.TemporaryDirectory() as folder:
            result = write_concordance_reports(folder, downstream, upstream, config)
            out = Path(folder)
            self.assertEqual(result["reference_version"], REF)
            self.assertTrue(all((out/f"similarity_ladder.{ext}").is_file() and (out/f"compatibility_table.{ext}").is_file()
                                for ext in ("png", "pdf", "svg")))
            svg = (out/"similarity_ladder.svg").read_text()
            self.assertIn(REF, svg)
            self.assertIn(SOURCE, svg)
            self.assertIn("SYNTHETIC TEST DATA", svg)
            self.assertIn("shared classes only*", svg)
            self.assertTrue(all((out/f"upstream_retention_ladder.{ext}").is_file() for ext in ("png", "pdf", "svg")))
            self.assertIn("SRC?", (out/"upstream_retention_ladder.svg").read_text())
            self.assertIn("CTRL", (out/"upstream_retention_ladder.svg").read_text())
            self.assertEqual(len(result["sample_plot_index"]), 2)
            for item in result["sample_plot_index"]:
                self.assertTrue(all((out/item[ext]).is_file() for ext in ("png", "pdf", "svg")))
                text = (out/item["svg"]).read_text()
                self.assertIn("SYNTHETIC TEST DATA", text)
                self.assertIn("no confidence interval", text)
                self.assertIn(REF, text)
                self.assertIn(SOURCE, text)
            self.assertIn('>M</text>', (out/result["sample_plot_index"][0]["svg"]).read_text())
            self.assertIn('>CLASS*</text>', (out/result["sample_plot_index"][0]["svg"]).read_text())
            with (out/"sample_similarity.tsv").open() as handle:
                row = next(r for r in csv.DictReader(handle, delimiter="\t") if r["sample_id"] == "SYNTHETIC_0" and r["target_version"] == REF)
            self.assertEqual(row["full_label_comparable"], "False")
            self.assertEqual(json.loads(row["common_supported_labels"]), ["ecDNA"])
            self.assertEqual(float(row["label_agreement"]), 1.)
            self.assertIn('sample_plot_index.tsv', (out/"report.md").read_text())
            self.assertIn("same_cycle_order", json.dumps(result))
            self.assertIn("人工合成测试数据", (out/"report.md").read_text())
            reports = list((out/"downstream").glob("*/*/report.md"))
            self.assertTrue(reports)
            self.assertTrue(all("SYNTHETIC TEST DATA" in p.read_text() for p in reports))
            with (out/"upstream_statistics.tsv").open() as handle:
                self.assertTrue(all(r["synthetic"] == "True" for r in csv.DictReader(handle, delimiter="\t")))


if __name__ == "__main__":
    unittest.main()
