"""Synthetic fixtures only: schema and error-path tests, no patient data."""
import json
import tempfile
import unittest
from pathlib import Path

from ascompare.parsing import LABELS, aa_amplicon_count, merge_intervals, parse_run


def write_fixture_run(path, labels=("ecDNA",), kind="Cyclic", count=1, modern=True):
    """Write minimal source-format artifacts without pretending to execute AA/AC."""
    path = Path(path)
    (path / "aa").mkdir(parents=True, exist_ok=True)
    (path / "ac").mkdir(exist_ok=True)
    (path / "stages.json").write_text(json.dumps({
        "aa": {"status": "complete", "returncode": 0},
        "ac": {"status": "complete" if count else "skipped_zero_amplicons", "returncode": 0}}))
    (path / "aa/reconstruction_summary.txt").write_text(f"#Amplicons = {count}\n")
    if not count:
        return path
    graph = "sequence\tchr1:100-\tchr1:199+\t8\t20\t100\t10\nsequence\tchr2:1000-\tchr2:1099+\t8\t20\t100\t10\n"
    cycles = "Segment\t1\tchr1\t100\t199\nSegment\t2\tchr2\t1000\t1099\nCycle=1;Copy_count=8;Segments=1+,2+\n"
    (path / "aa/reconstruction_amplicon1_graph.txt").write_text(graph)
    (path / "aa/reconstruction_amplicon1_cycles.txt").write_text(cycles)
    flags = ["ecDNA", "BFB"] + (["FAN"] if modern else [])
    fields = ["sample_name", "amplicon_number", "amplicon_decomposition_class"] + [x + "+" for x in flags]
    values = ["reconstruction", "amplicon1", kind] + ["Positive" if x in labels else "None detected" for x in flags]
    (path / "ac/classification_amplicon_classification_profiles.tsv").write_text("\t".join(fields) + "\n" + "\t".join(values) + "\n")
    bd = path / "ac/classification_classification_bed_files"
    bd.mkdir(exist_ok=True)
    for label in labels:
        (bd / f"reconstruction_amplicon1_{label}_1_intervals.bed").write_text("chr1\t100\t200\nchr2\t1000\t1100\n")
    return path


def add_synthetic_amplicon(path, kind="No amp/Invalid"):
    """Append an AA graph and AC row that may be negative, invalid, or lack its positive BED."""
    (path / "aa/reconstruction_summary.txt").write_text("#Amplicons = 2\n")
    for suffix in ("graph", "cycles"):
        text = (path / f"aa/reconstruction_amplicon1_{suffix}.txt").read_text()
        (path / f"aa/reconstruction_amplicon2_{suffix}.txt").write_text(text.replace("chr1", "chr3"))
    profile = next((path / "ac").glob("*profiles.tsv"))
    flags = ["Positive" if kind == "Cyclic" else "None detected", "None detected"]
    if "FAN+" in profile.read_text().splitlines()[0]:
        flags.append("None detected")
    with profile.open("a") as handle:
        handle.write("\t".join(["reconstruction", "amplicon2", kind, *flags]) + "\n")


class ParserTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def parse(self, supported=LABELS):
        return parse_run(self.root, "TEST_SAMPLE", "TEST_SEED", "TEST_VERSION", supported)

    def test_multilabel_and_connected_regions(self):
        write_fixture_run(self.root, ("ecDNA", "BFB", "FAN"))
        record = self.parse()
        self.assertEqual(record["status"], "ok")
        self.assertEqual(set(record["labels"]), {"ecDNA", "BFB", "FAN"})
        self.assertEqual(record["reconstruction_intervals"], [["chr1", 100, 200], ["chr2", 1000, 1100]])
        self.assertEqual(record["intervals"], record["reconstruction_intervals"])
        self.assertEqual(len(record["features"]), 3)

    def test_zero_amplicons_explicit_negative(self):
        write_fixture_run(self.root, count=0)
        self.assertEqual(aa_amplicon_count(self.root), 0)
        self.assertEqual(self.parse()["status"], "not_detected")

    def test_zero_without_success_is_not_negative(self):
        write_fixture_run(self.root, count=0)
        (self.root / "stages.json").write_text('{"aa":{"status":"interrupted"}}')
        self.assertEqual(self.parse()["status"], "failed")

    def test_missing_summary_is_missing(self):
        write_fixture_run(self.root)
        (self.root / "aa/reconstruction_summary.txt").unlink()
        self.assertIsNone(aa_amplicon_count(self.root))
        self.assertEqual(self.parse()["status"], "missing")

    def test_profile_rows_missing(self):
        write_fixture_run(self.root)
        (self.root / "ac/classification_amplicon_classification_profiles.tsv").unlink()
        self.assertEqual(self.parse()["status"], "missing")

    def test_missing_cycles_not_negative(self):
        write_fixture_run(self.root)
        (self.root / "aa/reconstruction_amplicon1_cycles.txt").unlink()
        record = self.parse()
        self.assertEqual((record["status"], record["stage"]), ("missing", "aa"))

    def test_legacy_invalid_cannot_be_called_negative(self):
        write_fixture_run(self.root, (), "No amp/Invalid", modern=False)
        record = self.parse([x for x in LABELS if x != "FAN"])
        self.assertEqual(record["status"], "invalid")
        self.assertIsNone(record["detected"])
        self.assertEqual(record["category_status"]["FAN"], "unsupported")

    def test_legacy_invalid_amplicon_keeps_other_valid_calls_as_partial_observations(self):
        write_fixture_run(self.root, modern=False)
        add_synthetic_amplicon(self.root)
        record = self.parse([x for x in LABELS if x != "FAN"])
        self.assertEqual((record['status'], record['coverage_status']), ('invalid', 'partial'))
        self.assertIsNone(record['detected'])
        self.assertEqual((record['evaluable_amplicon_count'], record['unevaluable_amplicon_count']), (1, 1))
        self.assertEqual([(a['amplicon_id'], a['status']) for a in record['amplicons']],
                         [('amplicon1', 'ok'), ('amplicon2', 'invalid')])
        self.assertEqual(record['features'][0]['label'], 'ecDNA')
        self.assertEqual(record['category_status']['BFB'], 'unavailable')

    def test_missing_cycle_or_positive_feature_preserves_other_amplicon(self):
        write_fixture_run(self.root)
        add_synthetic_amplicon(self.root, 'Cyclic')
        record = self.parse()
        self.assertEqual((record['status'], record['stage'], record['coverage_status']), ('missing', 'ac', 'partial'))
        self.assertEqual(record['amplicons'][0]['status'], 'ok')
        self.assertIn('feature BED', record['amplicons'][1]['reason'])
        (self.root / 'aa/reconstruction_amplicon2_cycles.txt').unlink()
        record = self.parse()
        self.assertEqual((record['status'], record['stage'], record['coverage_status']), ('missing', 'aa', 'partial'))
        self.assertEqual(record['amplicons'][0]['status'], 'ok')
        self.assertIn('cycles', record['amplicons'][1]['reason'])

    def test_modern_negative(self):
        write_fixture_run(self.root, (), "No-FSCNA")
        self.assertEqual(self.parse()["status"], "not_detected")

    def test_modern_invalid(self):
        write_fixture_run(self.root, (), "Invalid")
        self.assertEqual(self.parse()["status"], "invalid")

    def test_missing_feature_bed(self):
        write_fixture_run(self.root)
        next((self.root / "ac/classification_classification_bed_files").glob("*.bed")).unlink()
        self.assertEqual(self.parse()["status"], "missing")

    def test_unknown_category(self):
        write_fixture_run(self.root, (), "Future-new-category")
        self.assertEqual(self.parse()["status"], "unsupported")

    def test_invalid_bed_coordinates(self):
        write_fixture_run(self.root)
        next((self.root / "ac/classification_classification_bed_files").glob("*.bed")).write_text("chr1\t200\t100\n")
        self.assertEqual(self.parse()["status"], "invalid")

    def test_stale_summary_count(self):
        write_fixture_run(self.root)
        (self.root / "aa/reconstruction_summary.txt").write_text("#Amplicons = 0\n")
        self.assertEqual(self.parse()["status"], "missing")

    def test_missing_completion_evidence(self):
        write_fixture_run(self.root)
        (self.root / "stages.json").unlink()
        self.assertEqual(self.parse()["status"], "missing")

    def test_unknown_fan_flag_not_negative(self):
        write_fixture_run(self.root)
        p = self.root / "ac/classification_amplicon_classification_profiles.tsv"
        p.write_text(p.read_text().replace("None detected\tNone detected", "None detected\tNA"))
        self.assertEqual(self.parse()["status"], "unsupported")

    def test_overlaps_union(self):
        self.assertEqual(merge_intervals([["chr1", 0, 9], ["chr1", 7, 12], ["chr1", 12, 15]]), [["chr1", 0, 15]])


if __name__ == "__main__":
    unittest.main()
