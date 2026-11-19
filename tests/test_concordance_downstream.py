"""Stage similarity regressions using invented normalized and native artifacts."""
import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest

from ascompare.concordance import project_record
from ascompare.concordance_demo import _write_result
from ascompare.concordance_downstream import downstream_unit_similarity, write_downstream_similarity
from ascompare.parsing import parse_run
from ascompare.versions import load_lock


SOURCE = "suite-1.6.0"
REF = "suite-1.6.2"


def record(version=SOURCE, sample="SYNTHETIC_01", labels=None, support=None):
    labels = ["ecDNA"] if labels is None else labels
    support = ["ecDNA", "BFB", "FAN"] if support is None else support
    intervals = [["chr1", 10, 20]] if labels else []
    reconstructed = [["chr1", 0, 100]] if labels else []
    amplicons = [dict(amplicon_id="amplicon1", reconstruction_intervals=reconstructed,
                     intervals=intervals, labels=labels, status="ok", stage="result_reading")] if labels else []
    return dict(sample_id=sample, seed_id="candidate1", analysis_kind="aa_ac", version=version,
                source_version=SOURCE, comparison_eligible=True, status="ok" if labels else "not_detected",
                stage="result_reading", detected=bool(labels), labels=labels, supported_labels=support,
                intervals=intervals, features=[dict(amplicon_id="amplicon1", label=label, intervals=intervals) for label in labels],
                reconstruction_intervals=reconstructed, aa_amplicon_count=len(amplicons),
                amplicons=amplicons, coverage_status="complete", issues=[],
                run_dir="/SYNTHETIC_ONLY/run", evidence=["/SYNTHETIC_ONLY/run/stages.json"],
                component_versions={"aa": "SYNTHETIC", "ac": "SYNTHETIC"})


def compare(candidate, comparator):
    return downstream_unit_similarity(candidate, comparator, REF, SOURCE, ["input_source"])


class DownstreamSimilarityTests(unittest.TestCase):
    def test_raw_AA_AC_features_and_multilabels_are_independent_observations(self):
        source = record()
        target = record(REF, labels=["ecDNA", "BFB"])
        target["intervals"] = [["chr1", 10, 30]]
        target["features"] = [dict(label=label, intervals=target["intervals"]) for label in target["labels"]]
        result = compare(target, source)
        self.assertEqual(result["raw_aa_interval_jaccard"], 1.)
        self.assertEqual(result["ac_feature_interval_jaccard"], .5)
        self.assertEqual(result["full_label_agreement"], 0.)
        self.assertEqual(result["common_label_agreement"], 0.)
        self.assertEqual(result["target_labels"], ["BFB", "ecDNA"])
        self.assertEqual(result["target_traceability"]["evidence"], target["evidence"])

    def test_legacy_category_support_exposes_shared_scope_without_full_agreement(self):
        source = record(support=["ecDNA", "BFB"])
        target = record(REF, labels=["ecDNA", "FAN"])
        result = compare(target, source)
        self.assertFalse(result["full_label_comparable"])
        self.assertIsNone(result["full_label_agreement"])
        self.assertEqual(result["common_label_agreement"], 1.)
        self.assertEqual(result["common_supported_labels"], ["BFB", "ecDNA"])
        source = record(labels=[], support=["ecDNA", "BFB"])
        target = record(REF, labels=["FAN"])
        result = compare(target, source)
        self.assertIsNone(result["common_label_agreement"])
        self.assertEqual(result["label_comparison_reason"], "no_common_supported_positive_labels")

    def test_both_negative_and_unknown_support_do_not_create_label_agreement(self):
        result = compare(record(REF, labels=[]), record(labels=[]))
        self.assertIsNone(result["full_label_agreement"])
        self.assertIsNone(result["common_label_agreement"])
        self.assertIsNone(result["raw_aa_interval_jaccard"])
        self.assertEqual(result["raw_aa_comparison_reason"], "both_AA_interval_sets_empty")
        result = compare(record(REF, support=[]), record(support=[]))
        self.assertIsNone(result["full_label_agreement"])
        self.assertIsNone(result["common_label_agreement"])

    def test_raw_AA_empty_one_side_is_zero_overlap(self):
        result = compare(record(REF), record(labels=[]))
        self.assertEqual(result["raw_aa_interval_jaccard"], 0.)
        self.assertEqual(result["full_label_agreement"], 0.)

    def test_AC_failure_keeps_only_explicitly_complete_AA(self):
        candidate = record(REF)
        candidate.update(status="missing", stage="ac", comparison_eligible=False,
                         provenance_kind="imported_declared", declared_components_validated=True,
                         source_graph_identity="SYNTHETIC_AA_IDENTITY",
                         issues=[dict(stage="ac", status="missing", reason="SYNTHETIC missing profile")])
        candidate["amplicons"][0].update(status="missing", stage="ac")
        result = compare(candidate, record())
        self.assertFalse(result["evaluable"])
        self.assertIsNone(result["ac_feature_interval_jaccard"])
        self.assertIsNone(result["full_label_agreement"])
        self.assertEqual(result["raw_aa_interval_jaccard"], 1.)
        candidate["issues"].append(dict(stage="aa", status="missing"))
        self.assertIsNone(compare(candidate, record())["raw_aa_interval_jaccard"])
        candidate = record(REF)
        candidate["amplicons"] = []
        self.assertIsNone(compare(candidate, record())["raw_aa_interval_jaccard"])
        candidate = record(REF)
        candidate.update(status="failed", stage="aa")
        self.assertIsNone(compare(candidate, record())["raw_aa_interval_jaccard"])

    def test_complete_AA_does_not_override_comparison_ineligibility(self):
        candidate = record(REF)
        candidate["comparison_eligible"] = False
        self.assertIsNone(compare(candidate, record())["raw_aa_interval_jaccard"])
        candidate.update(status="unsupported", observed_status="ok", stage="comparison_input")
        result = compare(candidate, record())
        self.assertIsNone(result["raw_aa_interval_jaccard"])
        self.assertEqual(result["target_raw_aa_reason"], "same_input_comparison_not_eligible")
        candidate["comparison_eligible"] = True
        self.assertIsNone(compare(candidate, record())["raw_aa_interval_jaccard"])

    def test_AC_failure_needs_traceable_version_and_input(self):
        candidate = record(REF)
        candidate.update(status="missing", stage="ac", comparison_eligible=False,
                         issues=[dict(stage="ac", status="missing")])
        self.assertIsNone(compare(candidate, record())["raw_aa_interval_jaccard"])
        candidate.update(provenance_kind="imported_declared", declared_components_validated=True,
                         source_graph_identity="SYNTHETIC_VERIFIED_GRAPH")
        self.assertEqual(compare(candidate, record())["raw_aa_interval_jaccard"], 1.)
        candidate["declared_components_validated"] = False
        self.assertIsNone(compare(candidate, record())["raw_aa_interval_jaccard"])

    def test_fixed_seed_requires_actual_projection_for_every_provenance(self):
        seed = dict(chrom="chr1", start=10, end=20, seed_id="candidate1")
        source = project_record(record(), seed)
        for provenance in ("bounded_reconstruction", "imported_declared", "unknown"):
            with self.subTest(provenance=provenance):
                candidate = record(REF)
                candidate.update(provenance_kind=provenance, seed=seed)
                result = compare(candidate, source)
                self.assertIsNone(result["raw_aa_interval_jaccard"])
                self.assertEqual(result["target_raw_aa_reason"], "AA_not_projected_to_selected_seed")
        candidate = record(REF)
        candidate.update(status="missing", stage="ac", comparison_eligible=False,
                         provenance_kind="bounded_reconstruction", issues=[dict(stage="ac", status="missing")])
        projected = project_record(candidate, seed)
        self.assertNotIn("projection", projected)
        self.assertIsNone(compare(projected, source)["raw_aa_interval_jaccard"])
        # An explicit complete-AA projection may be supplied by a stage-aware
        # producer, but its coordinates must match the compared source seed.
        projected["projection"] = "overlapping_AA_reconstruction_complete_amplicons_retained"
        self.assertEqual(compare(projected, source)["raw_aa_interval_jaccard"], 1.)
        projected["seed"] = dict(seed, start=11)
        result = compare(projected, source)
        self.assertIsNone(result["raw_aa_interval_jaccard"])
        self.assertEqual(result["raw_aa_comparison_reason"], "raw_AA_comparison_input_mismatch")

    def test_raw_AA_same_sample_and_unit_are_required(self):
        for key in ("sample_id", "seed_id", "analysis_kind"):
            with self.subTest(key=key):
                candidate = record(REF)
                candidate[key] = "SYNTHETIC_DIFFERENT"
                self.assertIsNone(compare(candidate, record())["raw_aa_interval_jaccard"])

    def test_seed_projection_keeps_whole_amplicons_and_joint_summary_count(self):
        source = record()
        source["aa_amplicon_count"] = 5
        source["projection"] = "overlapping_AA_reconstruction_complete_amplicons_retained"
        self.assertEqual(compare(record(REF), source)["raw_aa_interval_jaccard"], 1.)
        source.update(amplicons=[], intervals=[], features=[], reconstruction_intervals=[], labels=[],
                      status="not_detected", detected=False)
        self.assertEqual(compare(record(REF), source)["raw_aa_interval_jaccard"], 0.)

    def test_AC_only_does_not_claim_AA_version_similarity(self):
        source, target = record(), record(REF)
        source["analysis_kind"] = target["analysis_kind"] = "ac_only"
        result = compare(target, source)
        self.assertIsNone(result["raw_aa_interval_jaccard"])
        self.assertEqual(result["raw_aa_comparison_reason"], "not_applicable_fixed_source_AA_not_rerun")
        self.assertEqual(result["full_label_agreement"], 1.)

    def test_cohort_is_patient_weighted_and_keeps_label_count_denominators(self):
        rows = []
        for sample, labels in (("SYNTHETIC_1", ["ecDNA"]), ("SYNTHETIC_2", ["BFB"]), ("SYNTHETIC_3", ["BFB"])):
            rows.extend([record(SOURCE, sample), record(REF, sample, labels)])
        patient_ids = {"SYNTHETIC_1": "SYNTHETIC_P1", "SYNTHETIC_2": "SYNTHETIC_P1", "SYNTHETIC_3": "SYNTHETIC_P2"}
        with tempfile.TemporaryDirectory() as folder:
            result = write_downstream_similarity(folder, {("aa_ac", SOURCE): rows}, [REF], REF, patient_ids, synthetic=True)
            cohort = next(row for row in result["cohort"] if row["target_version"] == REF and row["comparator_version"] == SOURCE)
            self.assertEqual(cohort["full_label_agreement"], .25)
            self.assertEqual(cohort["full_label_agreement_ratio"], "1/3")
            self.assertEqual(cohort["full_label_agreement_n_independent_units"], 2)
            self.assertIsNotNone(cohort["full_label_agreement_ci_high"])
            self.assertEqual(len(result["samples"]), 12)
            self.assertEqual(set(result), {"samples", "cohort", "files", "report_lines"})
            json.dumps(result, allow_nan=False)
            path = Path(folder)
            with (path/result["files"]["units"]).open() as handle:
                units = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(units), 12)
            self.assertIn("target_traceability", units[0])
            traces = json.loads((path/result["files"]["provenance"]).read_text())
            self.assertLess(len(traces), len(units))
            for unit in units:
                pointer = json.loads(unit["target_traceability"])
                self.assertEqual(pointer["file"], result["files"]["provenance"])
                self.assertIn(pointer["id"], traces)
            self.assertIn("完整标签一致率", (path/result["files"]["report"]).read_text())

    def test_fixed_reference_and_source_comparator_are_deduplicated(self):
        source = record(REF)
        source["source_version"] = REF
        with tempfile.TemporaryDirectory() as folder:
            result = write_downstream_similarity(folder, {("aa_ac", REF): [source]}, [REF], REF)
        self.assertEqual(len(result["samples"]), 1)
        row = result["samples"][0]
        self.assertEqual(row["comparator_roles"], ["fixed_reference", "input_source"])
        self.assertTrue(row["identity_anchor"])
        self.assertIsNone(result["cohort"][0]["full_label_agreement_ci_low"])

    def test_native_parsed_graph_coordinates_and_projected_empty_remain_traceable(self):
        entries = {entry["id"]: entry for entry in load_lock()["versions"]}
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            records = []
            for version, intervals in ((SOURCE, [["chr1", 10, 100], ["chr7", 200, 300]]),
                                       (REF, [["chr1", 10, 150], ["chr7", 200, 300]])):
                run = root/version
                _write_result(run, "SYNTHETIC_01", [dict(labels=["ecDNA"], intervals=intervals)], entries[version])
                supported = entries[version].get("supported_labels", entries[version].get("supported_classes"))
                parsed = parse_run(run, "SYNTHETIC_01", "candidate1", version, supported)
                parsed.update(analysis_kind="aa_ac", comparison_eligible=True, run_dir=str(run))
                records.append(parsed)
            result = compare(records[1], records[0])
            self.assertAlmostEqual(result["raw_aa_interval_jaccard"], 190/240)
            self.assertEqual(result["target_reconstruction_intervals"], [["chr1", 10, 150], ["chr7", 200, 300]])
            self.assertTrue(any(path.endswith("_graph.txt") for path in result["target_traceability"]["evidence"]))
            empty_seed = dict(chrom="chr2", start=10, end=30, seed_id="empty_candidate")
            left, right = (project_record(copy.deepcopy(row), empty_seed) for row in records)
            result = compare(right, left)
            self.assertIsNone(result["raw_aa_interval_jaccard"])
            self.assertEqual(result["raw_aa_comparison_reason"], "both_AA_interval_sets_empty")


if __name__ == "__main__":
    unittest.main()
