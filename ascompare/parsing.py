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
REQUIRED_PROFILE_FIELDS = {"sample_name", "amplicon_number", "amplicon_decomposition_class", "ecDNA+", "BFB+"}


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
    return _summary_count(files[0])


def _summary_count(path):
    hits = re.findall(r"^\s*#Amplicons\s*=\s*(\d+)\s*$", path.read_text(errors="replace"), re.M)
    return int(hits[0]) if len(hits) == 1 else None


def _amplicon_files(paths, artifact):
    files = {}
    for path in paths:
        match = re.fullmatch(rf"(.+)_(amplicon\d+)_{artifact}\.txt", path.name)
        if not match or match[2] in files:
            raise OutputError(f"Ambiguous AA {artifact} naming: {path.name}", stage="aa")
        files[match[2]] = path, match[1]
    return files


def _profile_header(reader, supported):
    columns = reader.fieldnames or []
    if len(set(columns)) != len(columns):
        raise OutputError("Duplicated AC profile columns", stage="ac")
    missing = REQUIRED_PROFILE_FIELDS.difference(columns)
    if missing:
        raise OutputError("AC profile header not supported: " + ",".join(sorted(missing)), "unsupported", "ac")
    if "FAN" in supported and "FAN+" not in columns:
        raise OutputError("Advertised FAN capability is missing from AC profile", "unsupported", "ac")
    return set(columns)


def _profile_amplicon(row):
    amp = row.get("amplicon_number")
    if not isinstance(amp, str) or not re.fullmatch(r"amplicon\d+", amp):
        raise OutputError("AC profile has an unrecognized amplicon identifier", "unsupported", "ac")
    return amp


def _profile_rows(path, supported):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = _profile_header(reader, supported)
        rows = {}
        for row in reader:
            amp = _profile_amplicon(row)
            if amp in rows:
                raise OutputError("Duplicate AC profile row: " + amp, stage="ac")
            rows[amp] = row
    return fields, rows


def _check_reference(intervals, lengths, allowed_missing, stage, description):
    if lengths is None:
        return
    for chrom, start, end in intervals:
        if chrom not in lengths or not 0 <= start < end <= lengths[chrom]:
            raise OutputError(description + " reaches a missing or out-of-bounds reference contig",
                              "unsupported" if chrom in allowed_missing else "invalid", stage)


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


def _read_graph_pair(graph, cycle, lengths, allowed_missing):
    if graph is None:
        raise OutputError("Missing AA graph", "missing", "aa")
    path, _ = graph
    if cycle is None or cycle[0].name != path.name.replace("_graph.txt", "_cycles.txt"):
        raise OutputError("Missing paired AA cycles file", "missing", "aa")
    cycle_path = cycle[0]
    if not cycle_path.stat().st_size:
        raise OutputError("Empty paired AA cycles file", "missing", "aa")
    text = cycle_path.read_text()
    if not re.search(r"^Segment\s", text, re.M) or not re.search(r"^Cycle=", text, re.M):
        raise OutputError("Cycles file lacks segment or cycle records", "invalid", "aa")
    try:
        intervals = _graph_intervals(path)
        _check_reference(intervals, lengths, allowed_missing, "aa", "AA extension")
    except OutputError as exc:
        raise OutputError(str(exc), exc.status, "aa") from exc
    return intervals


def _native_feature_map(path, beds, graph_by_basename):
    """Use only paths already selected by the importer, never open map entries."""
    bed_names = {}
    for bed in beds:
        bed_names.setdefault(bed.name, []).append(bed)
    mapped, seen, issues = {}, set(), []
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 2:
                issues.append(f"Malformed native feature-to-graph mapping at row {number}")
                continue
            bed_name, graph_name = (Path(value).name for value in fields)
            possible_beds = bed_names.get(bed_name, [])
            amp = graph_by_basename.get(graph_name)
            match = re.fullmatch(r"(.+)_(amplicon\d+)_(.+)_intervals\.bed", bed_name)
            if not match or len(possible_beds) != 1 or amp is None or amp != match[2]:
                issues.append(f"Native feature-to-graph mapping cannot verify row {number}")
                for bed in possible_beds:
                    seen.add(bed)
                    mapped.pop(bed, None)
                continue
            bed = possible_beds[0]
            if bed in seen:
                issues.append(f"Conflicting native feature-to-graph mapping at row {number}")
                mapped.pop(bed, None)
            else:
                mapped[bed] = amp, match[1], match[3]
                seen.add(bed)
    return mapped, seen.difference(mapped), issues


def _associate_beds(ac, profile, graphs):
    beds = sorted(ac.glob("*_classification_bed_files/*_intervals.bed"))
    map_path = profile.with_name(profile.name.replace("_amplicon_classification_profiles.tsv", "_features_to_graph.txt"))
    graph_by_basename = {path.name: amp for amp, (path, _) in graphs.items()}
    if map_path.is_file():
        mapped, ambiguous, issues = _native_feature_map(map_path, beds, graph_by_basename)
    else:
        mapped, ambiguous, issues = {}, set(), []
    features, aliases, unmatched = {}, set(), []
    for bed in beds:
        if bed in ambiguous:
            continue
        if bed in mapped:
            amp, prefix, feature = mapped[bed]
            aliases.add((amp, prefix))
        else:
            match = re.fullmatch(r"(.+)_(amplicon\d+)_(.+)_intervals\.bed", bed.name)
            if not match or match[2] not in graphs or graphs[match[2]][1] != match[1]:
                if "_unknown_" not in bed.name and bed.stat().st_size:
                    unmatched.append(bed.name)
                continue
            amp, feature = match[2], match[3]
        features.setdefault(amp, []).append((feature, bed))
    issues.extend("Unmatched classification BED: " + name for name in unmatched)
    return features, aliases, issues, map_path if map_path.is_file() else None


def _read_classification(amp, row, graph_prefix, feature_beds, alias, fields, supported, lengths, allowed_missing):
    if None in row or any(row.get(field) is None for field in fields):
        raise OutputError("Malformed AC profile row", stage="ac")
    profile_prefix = Path(row["sample_name"]).name
    if profile_prefix != graph_prefix and (amp, profile_prefix) not in alias:
        raise OutputError("AC sample prefix does not match AA graph or native mapping", stage="ac")
    labels = _profile_labels(row, fields, supported)
    found = set()
    features = []
    feature_ids = set()
    for feature_id, path in feature_beds:
        label = feature_id.rsplit("_", 1)[0]
        if label.lower().startswith("unknown"):
            continue
        if feature_id in feature_ids:
            raise OutputError("Duplicate feature IDs for AC amplicon: " + feature_id, stage="ac")
        feature_ids.add(feature_id)
        if row["amplicon_decomposition_class"] == "No-FSCNA":
            raise OutputError("No-FSCNA profile contradicts classified feature BED", stage="ac")
        if label not in supported:
            raise OutputError("Unrecognized feature BED class: " + label, "unsupported", "ac")
        labels.add(label)
        intervals = _bed_intervals(path)
        _check_reference(intervals, lengths, allowed_missing, "ac", "AC feature")
        if intervals:
            found.add(label)
            features.append(dict(amplicon_id=amp, feature_id=feature_id, label=label,
                                 intervals=intervals, source=str(path.resolve())))
    if labels - found:
        raise OutputError("Positive class lacks nonempty feature BED: " + ",".join(sorted(labels - found)), "missing", "ac")
    return labels, features


def _parse_amplicon(amp, graph, cycles, row, features, aliases, fields, supported, lengths, allowed_missing):
    record = dict(amplicon_id=amp, labels=[], intervals=[], reconstruction_intervals=[], detected=None,
                  decomposition_class=row.get("amplicon_decomposition_class") if row else None)
    stage = "aa"
    try:
        recon = _read_graph_pair(graph, cycles, lengths, allowed_missing)
        record["reconstruction_intervals"] = recon
        stage = "ac"
        if row is None:
            raise OutputError("Missing AC profile row", "missing", "ac")
        labels, parsed_features = _read_classification(amp, row, graph[1], features, aliases, fields,
                                                       supported, lengths, allowed_missing)
        record.update(status="ok" if labels else "not_detected", detected=bool(labels), labels=sorted(labels),
                      intervals=merge_intervals([i for f in parsed_features for i in f["intervals"]]),
                      stage="result_reading", reason="Complete AA/AC output for this amplicon")
        return record, parsed_features
    except OutputError as exc:
        record.update(status=exc.status, stage=stage, reason=str(exc))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        record.update(status="invalid", stage=stage, reason=f"{type(exc).__name__}: {exc}")
    return record, []


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
                  detected=None, coverage_status="unavailable", evaluable_amplicon_count=0,
                  unevaluable_amplicon_count=None, issues=[], supported_labels=sorted(supported),
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
        stage_names = ("input_check", "seed_processing", "aa", "ac")
        if not isinstance(stages, dict) or any(not isinstance(stages.get(name, {}), dict) for name in stage_names):
            raise OutputError("Stage completion evidence must contain named stage records")
        for stage in stage_names:
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
        graphs = _amplicon_files(sorted((root / "aa").glob("*_amplicon*_graph.txt")), "graph")
        cycles = _amplicon_files(sorted((root / "aa").glob("*_amplicon*_cycles.txt")), "cycles")
        if count != len(graphs) or count != len(cycles):
            result["issues"].append(dict(status="missing", stage="aa",
                reason="AA summary count does not match complete graph/cycles files"))
        if count == 0:
            if result["issues"]:
                raise OutputError(result["issues"][0]["reason"], "missing", "aa")
            result.update(status="not_detected", detected=False, stage="aa", coverage_status="complete",
                          unevaluable_amplicon_count=0,
                          reason="AA exited successfully and explicitly reported zero amplicons")
            result["category_status"] = {x: "not_detected" if x in supported else "unsupported" for x in LABELS}
            return result
        ac = stages.get("ac", {})
        profile_files = sorted((root / "ac").glob("*_amplicon_classification_profiles.tsv"))
        fields, rows, beds, aliases = set(), {}, {}, set()
        try:
            if ac.get("status") != "complete" or ac.get("returncode", 0) != 0:
                raise OutputError("Successful AC completion is not recorded", "missing", "ac")
            if len(profile_files) != 1:
                raise OutputError("Expected one AC classification profile", "missing", "ac")
            profile = profile_files[0]
            evidence(profile)
            fields, rows = _profile_rows(profile, supported)
            beds, aliases, map_issues, mapping = _associate_beds(root / "ac", profile, graphs)
            if mapping:
                evidence(mapping)
            for reason in map_issues:
                result["issues"].append(dict(status="invalid", stage="ac", reason=reason))
        except OutputError as exc:
            result["issues"].append(dict(status=exc.status, stage=exc.stage, reason=str(exc)))
            # Without a sound profile and association, individual AA outputs may
            # be recovered for review, while AC results remain indeterminate.
            rows = {}
        except (OSError, ValueError, TypeError, csv.Error) as exc:
            result["issues"].append(dict(status="invalid", stage="ac", reason=f"{type(exc).__name__}: {exc}"))
            rows = {}
        if rows and (len(rows) != count or set(rows) != set(graphs)):
            result["issues"].append(dict(status="missing", stage="ac",
                reason="AC rows do not match reconstructed amplicons exactly"))
        for amp in sorted(set(graphs) | set(cycles) | set(rows), key=lambda x: (int(x[len("amplicon"):]), x)):
            record, features = _parse_amplicon(amp, graphs.get(amp), cycles.get(amp), rows.get(amp),
                                               beds.get(amp, []), aliases, fields, supported,
                                               reference_lengths, allowed_missing_contigs)
            result["amplicons"].append(record)
            if record["reconstruction_intervals"]:
                result["reconstruction_intervals"].extend(record["reconstruction_intervals"])
                evidence(graphs[amp][0])
                evidence(cycles[amp][0])
            if record["status"] not in ("ok", "not_detected"):
                result["issues"].append(dict(amplicon_id=amp, status=record["status"], stage=record["stage"],
                                             reason=record["reason"]))
                continue
            result["evaluable_amplicon_count"] += 1
            result["labels"].extend(record["labels"])
            result["intervals"].extend(record["intervals"])
            for feature in features:
                result["evidence"].append(feature["source"])
            result["features"].extend(features)
        result["labels"] = sorted(set(result["labels"]))
        result["intervals"] = merge_intervals(result["intervals"])
        result["reconstruction_intervals"] = merge_intervals(result["reconstruction_intervals"])
        result["unevaluable_amplicon_count"] = max(count, len(result["amplicons"])) - result["evaluable_amplicon_count"]
        if result["issues"]:
            issues = result["issues"]
            # A structural AA failure must still block AC-only reclassification
            # when AC itself is pending.
            primary = next((i for i in issues if i["stage"] == "aa"), issues[0])
            result.update(status=primary["status"], stage=primary["stage"], reason=primary["reason"],
                          coverage_status="partial" if result["evaluable_amplicon_count"] else "unavailable")
            result["category_status"].update({x: "detected_in_evaluable_amplicons" for x in result["labels"]})
            return result
        result.update(status="ok" if result["labels"] else "not_detected",
                      detected=bool(result["labels"]), stage="result_reading", coverage_status="complete",
                      reason="Complete AA/AC output; full connected reconstruction retained")
        result["category_status"] = {x: "unsupported" if x not in supported else
                                     "detected" if x in result["labels"] else "not_detected" for x in LABELS}
    except OutputError as exc:
        result.update(status=exc.status, stage=exc.stage, reason=str(exc), detected=None)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        result.update(status="invalid", stage="result_reading", reason=f"{type(exc).__name__}: {exc}", detected=None)
    return result
