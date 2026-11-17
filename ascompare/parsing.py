"""Conservative adapters for AA graphs and AC 0.4/1.x/2.x profile outputs.

All normalized coordinates are BED 0-based, half-open. AA graph/Segment
endpoints are inclusive (AA ref_util.interval.size == end-start+1), hence
their right endpoint is increased by one. AC feature BEDs already use
half-open intervals. A missing file and the legacy 'No amp/Invalid' category
are never interpreted as a validated negative call.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path


LABELS = ("ecDNA", "BFB", "FAN", "Linear", "Complex-non-cyclic", "Virus")
LEGACY_LABELS = tuple(x for x in LABELS if x != "FAN")


class OutputError(ValueError):
    def __init__(self, message, status="invalid", stage="result_reading"):
        super().__init__(message)
        self.status = status
        self.stage = stage


def merge_intervals(intervals):
    """Union genomic intervals, preserving all chromosomes and remote regions."""
    merged = []
    for chrom, start, end in sorted(set(tuple(x) for x in intervals)):
        if start < 0 or end <= start:
            raise OutputError("Invalid genomic interval")
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    return merged


def aa_amplicon_count(run_dir):
    """Return the explicit AA summary count, or None for absent/invalid summary.

    The caller must additionally check successful AA exit before skipping AC.
    No counts are inferred from missing graphs or empty BEDs.
    """
    files = sorted((Path(run_dir) / "aa").glob("*_summary.txt"))
    if len(files) != 1:
        return None
    hits = re.findall(r"^\s*#Amplicons\s*=\s*(\d+)\s*$", files[0].read_text(errors="replace"), re.M)
    return int(hits[0]) if len(hits) == 1 else None


def _graph_intervals(path):
    intervals = []
    for line in path.read_text().splitlines():
        if not line.startswith("sequence"):
            continue
        fields = line.split()
        if len(fields) < 4:
            raise OutputError("Malformed sequence row: " + str(path))
        left = re.fullmatch(r"(.+):(\d+)[+-]", fields[1])
        right = re.fullmatch(r"(.+):(\d+)[+-]", fields[2])
        if not left or not right or left[1] != right[1] or int(right[2]) < int(left[2]):
            raise OutputError("Malformed graph coordinates: " + str(path))
        intervals.append([left[1], int(left[2]), int(right[2]) + 1])
    if not intervals:
        raise OutputError("Graph has no sequence intervals: " + str(path))
    return merge_intervals(intervals)


def _bed_intervals(path):
    intervals = []
    for n, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split()
        try:
            intervals.append([fields[0], int(fields[1]), int(fields[2])])
        except (IndexError, ValueError) as exc:
            raise OutputError(f"Malformed feature BED {path}:{n}") from exc
    return merge_intervals(intervals)


def _flag(value, column):
    if value == "Positive":
        return True
    if value == "None detected":
        return False
    # NA in an advertised class is unavailable evidence, not a negative.
    raise OutputError(f"Unsupported/unavailable {column} value: {value!r}", "unsupported")


def _profile_labels(row, fields, supported):
    flags = {label for label in ("ecDNA", "BFB", "FAN")
             if label + "+" in fields and _flag(row[label + "+"], label + "+")}
    kind = row["amplicon_decomposition_class"]
    if kind in ("No amp/Invalid", "Invalid"):
        # 1.x combines nonamplified and invalid; do not invent a split.
        raise OutputError("AC legacy No amp/Invalid is ambiguous" if kind == "No amp/Invalid"
                          else "AC classified reconstruction as Invalid", "invalid", "ac")
    if kind == "No-FSCNA":
        if flags:
            raise OutputError("Contradictory No-FSCNA and positive mechanism flags")
        return set()
    if kind in ("Linear", "Complex-non-cyclic", "Virus") and not flags:
        flags.add(kind)
    elif kind not in ("Linear", "Complex-non-cyclic", "Virus", "Cyclic"):
        raise OutputError("Unrecognized AC decomposition class: " + kind, "unsupported")
    if not flags:
        raise OutputError("Cyclic class has no resolved ecDNA/BFB/FAN mechanism", "unsupported")
    if flags - supported:
        raise OutputError("Output contains a class not declared in version capabilities: "
                          + ",".join(sorted(flags - supported)), "unsupported")
    return flags


def parse_run(run_dir, sample_id, seed_id, version, supported_labels=None, reference_lengths=None,
              allowed_missing_contigs=()):
    """Normalize one seed-anchored run, including all connected AA outputs.

    Expected directories: aa/ and ac/; stages.json is written by the runner.
    `intervals` are the union of positive classified feature BEDs;
    `reconstruction_intervals` preserve all reconstructed AA segments.
    `amplicons`/`features` retain IDs and labels rather than choosing one winner.
    Unknown category support is explicit per category, never coded as absence.
    """
    root = Path(run_dir)
    supported = set(LEGACY_LABELS if supported_labels is None else supported_labels)
    result = dict(sample_id=sample_id, seed_id=seed_id, version=version,
                  status="missing", labels=[], intervals=[], reconstruction_intervals=[],
                  amplicons=[], features=[], evidence=[], stage="result_reading", reason="",
                  detected=None, supported_labels=sorted(supported),
                  unsupported_labels=sorted(set(LABELS) - supported),
                  category_status={x: "unavailable" if x in supported else "unsupported" for x in LABELS})

    def evidence(path):
        result["evidence"].append(str(path.resolve()))

    try:
        stage_path = root / "stages.json"
        if not stage_path.is_file():
            raise OutputError("Stage completion record is missing; output presence cannot establish success", "missing")
        evidence(stage_path)
        stages = json.loads(stage_path.read_text())
        for stage in ("input_check", "seed_processing", "aa", "ac"):
            state = stages.get(stage, {})
            if state.get("status") in ("failed", "interrupted", "unsupported"):
                raise OutputError(state.get("reason") or state.get("error") or "Execution did not complete",
                                  "unsupported" if state["status"] == "unsupported" else "failed", stage)
        aa = stages.get("aa", {})
        if aa.get("status") != "complete" or aa.get("returncode", 0) != 0:
            raise OutputError("Successful AA completion is not recorded", "missing", "aa")
        count = aa_amplicon_count(root)
        if count is None:
            raise OutputError("AA summary missing, duplicated, or count unrecognized", "missing", "aa")
        for p in (root / "aa").glob("*_summary.txt"):
            evidence(p)
        result["aa_amplicon_count"] = count
        graphs = sorted((root / "aa").glob("*_amplicon*_graph.txt"))
        cycles = sorted((root / "aa").glob("*_amplicon*_cycles.txt"))
        if count != len(graphs) or count != len(cycles):
            raise OutputError("AA summary count does not match complete graph/cycles files", "missing", "aa")
        if count == 0:
            result.update(status="not_detected", detected=False, stage="aa",
                          reason="AA exited successfully and explicitly reported zero amplicons")
            result["category_status"] = {x: "not_detected" if x in supported else "unsupported" for x in LABELS}
            return result
        graph_by_amp = {}
        for p in graphs:
            match = re.fullmatch(r"(.+)_(amplicon\d+)_graph\.txt", p.name)
            if not match or match[2] in graph_by_amp:
                raise OutputError("Ambiguous AA graph naming")
            cp = p.with_name(p.name.replace("_graph.txt", "_cycles.txt"))
            if not cp.is_file() or cp.stat().st_size == 0:
                raise OutputError("Missing/empty paired cycles file", "missing", "aa")
            # Segment/cycle evidence is needed, not merely an empty placeholder.
            ctext = cp.read_text()
            if not re.search(r"^Segment\s", ctext, re.M) or not re.search(r"^Cycle=", ctext, re.M):
                raise OutputError("Cycles file lacks segment or cycle records", "invalid", "aa")
            ints = _graph_intervals(p)
            graph_by_amp[match[2]] = (match[1], ints)
            result["reconstruction_intervals"].extend(ints)
            evidence(p)
            evidence(cp)
        result["reconstruction_intervals"] = merge_intervals(result["reconstruction_intervals"])
        if reference_lengths is not None:
            for chrom, start, end in result["reconstruction_intervals"]:
                if chrom not in reference_lengths or not 0 <= start < end <= reference_lengths[chrom]:
                    raise OutputError("AA extension reaches a missing or out-of-bounds reference contig",
                                      "unsupported" if chrom in allowed_missing_contigs else "invalid", "aa")
        ac = stages.get("ac", {})
        if ac.get("status") != "complete" or ac.get("returncode", 0) != 0:
            raise OutputError("Successful AC completion is not recorded", "missing", "ac")
        profile_files = sorted((root / "ac").glob("*_amplicon_classification_profiles.tsv"))
        if len(profile_files) != 1:
            raise OutputError("Expected one AC classification profile", "missing", "ac")
        profile = profile_files[0]
        evidence(profile)
        with profile.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = set(reader.fieldnames or [])
            required = {"sample_name", "amplicon_number", "amplicon_decomposition_class", "ecDNA+", "BFB+"}
            if not required <= fields:
                raise OutputError("AC profile header not supported: " + ",".join(sorted(required - fields)), "unsupported")
            if "FAN" in supported and "FAN+" not in fields:
                raise OutputError("Advertised FAN capability is missing from AC profile", "unsupported")
            rows = list(reader)
        if len(rows) != count or {x["amplicon_number"] for x in rows} != set(graph_by_amp):
            raise OutputError("AC rows do not match reconstructed amplicons exactly", "missing", "ac")
        beds = sorted((root / "ac").glob("*_classification_bed_files/*_intervals.bed"))
        seen = set()
        for row in rows:
            amp = row["amplicon_number"]
            prefix, recon = graph_by_amp[amp]
            if Path(row["sample_name"]).name != prefix:
                raise OutputError("AC sample prefix does not match AA graph")
            labels = _profile_labels(row, fields, supported)
            amp_intervals = []
            found = set()
            for p in beds:
                feature_prefix = prefix + "_" + amp + "_"
                if not p.name.startswith(feature_prefix):
                    continue
                feature_id = p.name[len(feature_prefix):-len("_intervals.bed")]
                label = feature_id.rsplit("_", 1)[0]
                if label.lower().startswith("unknown"):
                    continue
                if row["amplicon_decomposition_class"] == "No-FSCNA":
                    raise OutputError("No-FSCNA profile contradicts classified feature BED")
                if label not in labels:
                    # Preserve additional simultaneously classified feature types.
                    if label not in supported:
                        raise OutputError("Unrecognized feature BED class: " + label, "unsupported")
                    labels.add(label)
                ints = _bed_intervals(p)
                if reference_lengths is not None:
                    for chrom, start, end in ints:
                        if chrom not in reference_lengths or not 0 <= start < end <= reference_lengths[chrom]:
                            raise OutputError("AC feature reaches a missing or out-of-bounds reference contig",
                                              "unsupported" if chrom in allowed_missing_contigs else "invalid", "ac")
                if not ints:
                    continue
                found.add(label)
                amp_intervals.extend(ints)
                evidence(p)
                result["features"].append(dict(amplicon_id=amp, feature_id=feature_id,
                                               label=label, intervals=ints, source=str(p.resolve())))
                seen.add(p)
            if labels - found:
                raise OutputError("Positive class lacks nonempty feature BED: " + ",".join(sorted(labels - found)), "missing", "ac")
            result["labels"].extend(labels)
            result["intervals"].extend(amp_intervals)
            result["amplicons"].append(dict(amplicon_id=amp, labels=sorted(labels),
                                           decomposition_class=row["amplicon_decomposition_class"],
                                           intervals=merge_intervals(amp_intervals), reconstruction_intervals=recon))
        # Unmatched positive BEDs indicate stale/mixed outputs, not extra calls.
        for p in beds:
            if p not in seen and "_unknown_" not in p.name and p.stat().st_size:
                raise OutputError("Unmatched classification BED: " + p.name)
        result["labels"] = sorted(set(result["labels"]))
        result["intervals"] = merge_intervals(result["intervals"])
        result.update(status="ok" if result["labels"] else "not_detected",
                      detected=bool(result["labels"]), stage="result_reading",
                      reason="Complete AA/AC output; full connected reconstruction retained")
        result["category_status"] = {x: "unsupported" if x not in supported else
                                     "detected" if x in result["labels"] else "not_detected" for x in LABELS}
    except OutputError as exc:
        result.update(status=exc.status, stage=exc.stage, reason=str(exc), detected=None)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        result.update(status="invalid", stage="result_reading", reason=f"{type(exc).__name__}: {exc}", detected=None)
    return result
