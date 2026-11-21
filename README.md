# AmpliconSuite Cross-Version Concordance

**Assess how AmpliconSuite results change across versions at selected regions.**

[中文说明](README.zh-CN.md) · [Examples](examples) · [Environment setup](docs/environment_setup.md) · [Statistical methods](docs/statistical_methods.md) · [Changelog](#changelog)

AmpliconSuite Cross-Version Concordance is an independent, local research tool. It compares existing AmpliconSuite outputs, orchestrates selected-seed AA/AC runs, and compares native seed discovery from whole-genome CNVkit. It produces per-sample similarity ladders, paired difference tables and compatibility reports.

**Version:** 0.5.0.dev1 · **Fixed reference:** Suite **1.6.2**, AA **1.6.r2**, AC **2.0.1**. The reference is frozen; installing a newer upstream release does not change it. This project is not an official AmpliconSuite component.

> Similarity describes observed software outputs. It is not a probability of biological truth or proof that complete workflows are interchangeable. Compatibility remains **insufficient evidence** when margins are unspecified or evidence is inadequate.

## Global CNVkit discovery with BAM

New manifests use `"cnv_scope": "global"`. Each Suite obtains its recipe's whole-genome CNVkit segments, applies its own seed filters and runs AA/AC on its native seeds at the sampled loci. Matching whole-genome measurements are reused across Suite versions. Completed discovery can find seeds even when the original seed BED was empty. Manifests that explicitly configured legacy local windows keep that scope.

[Global execution, reuse and reporting](docs/global_cnv.md) · [Manifest example](examples/concordance.with_bam.json)

To retry AA jobs that exceeded the unique-region limit, use the completed global run as `recover_aa_from` and set `aa_rescue_unique_bp` in a separate run. Recovery uses the original discovery, seed selection and verified completed AA/AC results from that run. [Recovery settings](docs/global_cnv.md#aa-读取预算恢复)

## Official Suite defaults for selected regions

New manifests use `"parameter_policy": "suite_defaults"` at the top level. AA runs in FULL/EXPLORE with the Suite default 10× target depth. The official version's licensing check selects the AA solver. Concordance leaves the inherited Python hash environment in place. AC uses its native classification thresholds. Native seed-filter defaults are gain 4.5 and minimum size 50,000 bp. Concordance checks optional coverage statistics against the locked AA reader before staging them. If AA recomputes coverage, valid graph and cycle outputs can still proceed to AC. The cache decision and observed outcome are saved in `aa_background.json`.

In the legacy local workflow the selected input regions, bounded local CNVkit design and 100 Mb access cap define the study scope. Genome-wide discovery and full-BAM QC are outside that scope. Remove any `aa_parameters` and `ac_parameters` overrides, and runtime solver or hash-seed overrides, when enabling this policy. Use a fresh output directory.

## Faster fixed-seed reconstruction (0.4.3)

The custom-parameter policy retains the 0.4.3 default AA `CYCLES` mode with `EXPLORE`. It writes the graph, cycles and summary for classification and report generation. To also generate native AA SV figures, set `aa_parameters.runmode` to `FULL`.

Matching AA reconstructions can serve multiple Suite recipes. Reuse requires the same verified AA commit, inputs, effective parameters, runtime, background and access policy. Each Suite uses its configured classifier. The per-Suite attempt records identify the producer of reused results. Use `stage_timing.tsv` for timing and cache activity. The fixed reference, seed selection, per-sample figures and cohort endpoints retain their existing definitions.

Each `invocation_id` in `stage_timing.tsv` identifies one invocation. Read its elapsed time from the `invocation_total` row. Filter `measurement_kind=subprocess` for actual AA, AC and local CNVkit subprocess times. Native AA diagnostics can overlap. Rows marked `producer_reference` give historical timing for a reused artifact. Count the invocation total once.

Start in a fresh output directory when upgrading. The first run measures cold performance. Repeat the same command to resume. To change AC parameters, choose a new output directory and optionally set the manifest field `reuse_aa_from` to a prior 0.4.3 run directory. Keep the previous output available until the new run is complete. A budget failure remains visible for any AA job that exceeds the bounded access policy.

## Sample results and source-matched CNVkit controls (0.4.2)

In legacy local mode each sample has its own local CNVkit standard: the version used to generate its input. For example, a sample originally processed with CNVkit 0.9.14 uses a local 0.9.14 run as its standard. Each target CNVkit uses the same tumor, normal-reference strategy, candidate windows, background, binning, segmentation settings and **source Suite AA seed filter**. Each distinct CNVkit version is evaluated once per sample in this comparison.

The original Suite result directory can supply the version via a single `*_run_metadata.json` recording the CNVkit command and version. Alternatively, declare the exact `source_cnvkit_version` in the sample manifest or TSV. Conflicting evidence yields an unavailable standard. Configure its exact executable in the shared `cnvkit_runtimes` registry if it is absent from the configured Suite runtimes. Optional `local_cnv.comparison_versions` adds standalone CNVkit targets. [Configuration and controls](docs/local_cnv_calibration.md)

Open `result.md`. Read the retained/evaluable counts for the standard and target alongside their **paired interval Jaccard** within original candidate regions. A high Jaccard can coexist with low retention of original candidates. `local_cnv_comparator_samples.tsv` gives paired counts, gains, losses and target-minus-standard retention change. Its cohort counterpart includes pointwise descriptive intervals with patient grouping when declared. A missing runtime, unknown source version or pair with no interval signal is explicitly identified. The Suite panel keeps its fixed reference and evaluates its own native CNVkit and AA combinations. Conditional Suite compatibility and the original-run retention check have separate evidence requirements.

Downstream AA and AC interval Jaccards and full/shared class label agreement are in `result.tsv`. Generate a report in a fresh directory from verified 0.4.3 or 0.4.2 records or earlier 0.3.0–0.4.1 supported releases. Report replay reads saved results only. Older runs display `not_recorded_in_run` for the new local control. To measure it, start a new bounded `with_bam` run.

```sh
asconcordance report --from-run runs/comparison --out reports/comparison-0.4.3
```

[Source replay and local CNV calibration protocol](docs/source_reproduction.md)

## What it does

- Checks input files, coordinates, declared versions, reference compatibility and runtime identity.
- Runs pinned components on paired loci with version-specific native seeds, records failures and resumes verified cached work.
- Preserves multiple amplification classes and distinguishes negative, failed, missing and unsupported results.
- Reports detection, additional/missed calls, class changes, interval overlap and class-specific interval overlap.
- Compares paired samples, accounts for multiple seeds within a sample, and reports effects, confidence bounds and applicable tests.
- Locates failures at input checking, seed processing, AA reconstruction, AC classification or result parsing.

## Two input modes

| Mode | Required input | Supported analysis |
|---|---|---|
| `with_bam` | Existing Suite results, source version, indexed BAM and reference resources | Global CNVkit, native Suite seeds and paired AA/AC reconstruction |
| `results_only` | Existing AA/AC outputs, source Suite version and reference-build declaration | Compare supplied multi-version results; optionally rerun different AC versions on the same source AA reconstruction |

**Without BAM, a single AA result cannot be used to rerun other AA versions.** AC-only comparisons require the same source AA version and graph identity and are reported separately from AA+AC comparisons. Missing versions remain unverified. Imported version declarations are checked against the lock; they are not proof of the historical execution environment.

Local CNVkit requires an explicit background BED. The default `matched_normal` strategy also requires a matched normal BAM and index. Tumor-only input can explicitly select `local_cnv.reference_mode="flat"`; there is no automatic fallback. Neither flanking sequence nor a flat reference is automatically a valid normal baseline.

## Quick start

The Python distribution is `ampliconsuite-cross-version-concordance`. Install it in a fresh virtual environment and run `asconcordance`.

Run these commands from the extracted repository root on Linux or macOS:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install .
asconcordance --help
asconcordance versions
asconcordance demo --out demo_output
```

Open `demo_output/report.md`. The demo creates **synthetic result files** and uses the public comparison/reporting workflow. It needs no patient data, BAM, AA or AC installation. **It does not execute AA or AC and is not a biological validation dataset.** Use a new output path for another demo.

Requirements: Python **3.10+**, NumPy, SciPy, Matplotlib and pysam. Installing this package installs the comparison dependencies; it does **not** install AmpliconSuite components, reference resources or commercial solver licenses. The Python module entry point is `python -m ascompare.concordance`.

### Configure a real comparison

Start with one of these templates:

- [`concordance.results_only.json`](examples/concordance.results_only.json): compare multiple existing outputs without BAM or algorithm environments.
- [`concordance.ac_only.json`](examples/concordance.ac_only.json): rerun AC on one existing AA reconstruction; requires the specified AC environments and reference resources.
- [`concordance.with_bam.json`](examples/concordance.with_bam.json): local CNVkit and selected-seed AA/AC; requires explicit component runtimes.

```sh
cp examples/concordance.results_only.json concordance.local.json
# Edit concordance.local.json with your local paths and version declarations.
asconcordance check --manifest concordance.local.json --out runs/comparison
asconcordance run --manifest concordance.local.json --out runs/comparison
asconcordance report --out runs/comparison
```

Paths are resolved relative to the manifest. Output must be separate from inputs and resources. Results can be specified with `aa_dir`/`ac_dir`, a `run_dir` containing `aa` and `ac`, or a recognizable Suite `result_dir`; ambiguous directories are rejected. Keep the original graph, cycles, summary, classification profiles and feature BED files where available. An image or summary label alone is insufficient.

`check` inspects input declarations, coordinates, BAM headers, existing indexes and small result-file formats. `import_format_check.tsv` lists AA/AC file counts, headings and matching issues. `run` also checks real component versions and source identities. `report` redraws saved, integrity-checked records from the same program version without accessing BAM. Use a new output directory for runs started with an earlier release. A completed command/report does not imply that every scientific job succeeded: inspect coverage and failure tables.

See the [server guide](docs/server_quickstart.md), [environment instructions](docs/environment_setup.md) and [local CNV calibration requirements](docs/local_cnv_calibration.md) for runtime configuration. These detailed technical notes are currently in Chinese; commands and configuration keys are shared across both languages.

### Batch input and cohort reports (0.4.0)

Import a local sample TSV and a shared JSON config, then analyze each sample in first-appearance order:

```sh
asconcordance import-samples --samples samples.tsv --config batch.json --out batch.manifest.json
asconcordance check --manifest batch.manifest.json --out runs/batch01
asconcordance run --manifest batch.manifest.json --out runs/batch01
```

Templates: [with BAM sample table](examples/samples.with_bam.tsv) + [config](examples/batch.with_bam.json), or [existing results table](examples/samples.results_only.tsv) + [config](examples/batch.results_only.json). These contain synthetic placeholder paths; configure your actual inputs and runtimes. Sample and output-folder names are your choice. See the [full input schema, path rules and batch guide](docs/batch_input.md).

Open `sample_reports.md` for all sample plots and `report.md` for the cohort. The new `cohort_compatibility.*` figure shows the actual statistical decision endpoints, bounds, prespecified margins, independent-unit counts and conditional final decisions. `cohort_compatibility_metrics.tsv` preserves its exact data. Patient IDs group linked samples; source versions and AC-only analyses stay separate. Per-sample plots show observed differences, not formal single-sample compatibility. `sample_progress.tsv` tracks serial processing and retains failed samples.

## Legacy local computation

This section describes explicit local-window manifests. It does not run CNVkit `batch`/`autobin`, discover seeds genome-wide, sort or index BAM, convert CRAM, run full-file BAM statistics, or read an entire BAM to calculate its hash. Supply a sorted BAM and index prepared outside this program.

Upstream analysis measures all supplied candidate regions plus finite flanks and explicit background windows, then runs the corresponding AA seed filter. Downstream analysis uses the same frozen random subset across versions. Seeds are candidates, not confirmed amplifications.

| Setting | Default | Meaning |
|---|---:|---|
| `sampling.max_seeds` | 10 per sample | Maximum downstream seed count; use all if fewer are available |
| `sampling.random_seed` | 20261117 | Reproducible selection before examining rerun results |
| `local_cnv.flank_bp` | 100,000 | Extension on each side of a candidate |
| `local_cnv.bin_size` | 1,000 | Local CNV coverage bin size |
| `local_cnv.max_window_bp` | 25,000,000 | Maximum merged local window length |
| `local_cnv.min_seed_overlap` | 1.0 | Required fraction of an original candidate covered by retained seed intervals |
| `statistics.min_samples` | 20 | Conservative operational minimum, not a guarantee of statistical power |
| Compatibility margins | Unset | No compatibility acceptance without prespecified allowed differences |

These are configurable engineering starting points, **not official compatibility thresholds or a validated optimal seed count**. Candidate coverage is a matching rule, not a statistical acceptance margin.

AA retains finite indexed background sampling and extension into connected regions. The BAM is not cropped. If reading limits are exceeded, the run fails explicitly instead of returning a truncated structure as successful. Native CNVkit `bedcov` is bounded by finite BED/query lengths and process timeouts; its internal read count is not directly observable. The guard is intended for audited scientific callers, not hostile-code isolation.

Source-Suite reruns are compared with the supplied original results. Upstream calibration uses a separate same-version source control and traceable original evidence. Target-version losses are observations. Combined decisions require the appropriate source checks and shared local analysis conditions. The recomputed local seeds are not all fed into AA again, so the combined report is **stepwise conditional evidence**, not an end-to-end discovery benchmark.

## Locked component panel

| Suite | AA | AC | BFBArchitect | Container-recipe CNVkit |
|---|---|---|---|---|
| 1.2.2 | 1.3.r8 | 1.1.2 | — | 0.9.10 |
| 1.3.3 | 1.4.r2 | 1.2.0 | — | 0.9.10 |
| 1.3.4 | 1.4.r2 | 1.2.1 | — | 0.9.10 |
| 1.3.7 | 1.5.r2 | 1.3.2 | — | 0.9.10 |
| 1.4.0 | 1.5.r3 | 1.4.2 | — | 0.9.10 |
| 1.5.0 | 1.5.r4 | 1.5.1 | — | 0.9.10 |
| 1.5.2 | 1.5.r6 | 1.5.2 | — | 0.9.10 |
| 1.5.3 | 1.5.r7 | 1.5.2 | — | 0.9.13 |
| 1.6.0 | 1.6.r0 | 2.0.0 | 1.0.1 | 0.9.13 |
| 1.6.1 | 1.6.r1 | 2.0.0 | 1.0.1 | 0.9.13 |
| **1.6.2 — reference** | **1.6.r2** | **2.0.1** | **1.0.2** | **0.9.13** |

These are explicit comparison component recipes. Source commits, archive hashes and source-manifest hashes are recorded in [`versions.lock.json`](versions.lock.json); CNVkit build-recipe evidence is in [`builtin_cnvkit_map.json`](ascompare/builtin_cnvkit_map.json). A Suite label or image name alone does not prove its installed components. Historical versions are never silently replaced. See [version evidence](docs/version_evidence.md) for qualifications, including Suite 1.3.4 and the separate legacy diagnostic recipe.

CNVkit 0.9.14 can serve as the local source standard when it matches the original input, or as an additional controlled target. Historical AA requires an appropriate Mosek installation/license. See the validation record for the separate evidence levels: source lock, runtime launch, synthetic component test and real-sample workflow. For a downsampling reproducibility comparison, see the [cohort study guide](docs/cohort_study.md).

## Similarity and compatibility

The main ladders compare each input **source version → target version**, with source strata kept separate. Reports also contain comparisons against the **fixed reference** and other version pairs. Source self-comparisons are identity anchors, not independent validation.

Four downstream endpoints are kept separate: detection agreement, full label-set agreement, interval Jaccard and class-specific interval Jaccard. Unsupported classes are explicit; partial shared-class comparisons cannot appear as complete classification agreement. Upstream reports candidate retention separately.

The report opens with paired direction (`target − source`), detection change and an interval, added/missed/changed counts, and eligible denominators. Upstream also shows paired candidate retention on the candidates measurable by both versions. Add `patient_id` to every sample for patient-level grouping. Seeds are averaged within samples, then samples within each patient, and patients are weighted equally. Agreement decisions use prespecified loss margins and conservative one-sided bounds with multiple-comparison control. For a learning cohort, set `statistics.primary_pairs` to `[]` to produce descriptive results. For validation, prespecify the selected version pairs and margins on a separate set of patients. The [cohort study guide](docs/cohort_study.md) explains how to plan margins and avoid choosing thresholds on the validation data.

| Decision | Interpretation |
|---|---|
| `within_margin` | Evidence supports agreement within prespecified differences for the evaluated scope |
| `exceeds_margin` | Evidence supports differences exceeding the prespecified range for the stated endpoint |
| `insufficient_evidence` | Missing margins, too few independent samples, incomplete outputs or other unmet evidence conditions |

Execution failures, missing files and unsupported classes remain distinct from these statistical decisions. High similarity cannot establish identical path order, breakpoint topology, biological truth or general software interchangeability. See [statistical methods](docs/statistical_methods.md).

## Reports

![Synthetic per-sample similarity ladder](docs/assets/synthetic_sample_ladder.png)

*Illustration generated from invented result files. No patient data or genuine AA/AC execution is represented.*

![Synthetic cohort compatibility figure](docs/assets/synthetic_cohort_compatibility.png)

| Output | Contents |
|---|---|
| `report.md`, `summary.json` | Paired effects, denominators, source/reference identities and evidence |
| `import_format_check.tsv` | Small-file format audit from `asconcordance check` |
| `similarity_ladder.*`, `sample_ladders/` | Aggregate and per-sample ladders in PNG/PDF/SVG |
| `sample_plot_index.tsv`, `sample_similarity.tsv` | Plot paths, observed similarities and evaluation status |
| `compatibility_table.*`, `compatibility_combined.tsv` | Three-state compatibility summaries |
| `upstream_retention_ladder.*`, `upstream_*.tsv` | Local retention plots, seed/sample statistics and calibration state |
| `cohort_compatibility.*`, `cohort_compatibility_metrics.tsv` | Cohort endpoint differences, bounds, margins, independent n and final decisions |
| `sample_progress.json`, `sample_progress.tsv` | Serial sample order, state and result status counts |
| `downstream/` | Stratified comparisons, differences, statistics and compatibility decisions |
| `local_cnv_comparator_{seeds,samples,cohort}.tsv`, `local_cnv_comparator_cohort.*` | Source-matched CNVkit contrasts, evaluable counts, descriptive intervals and cohort figure |
| `sampling/`, `source_bridge/` | Frozen seed selection and source-result reproducibility checks |
| `imports/`, `reclassifications/`, `*_jobs/`, `cnv_cache/` | Staged small results, intermediate files, logs and cache evidence |
| `provenance.json`, `input_identity.json`, `runtime_check.json`, `versions.lock.json`, `run_metadata.tsv` | Actual versions, parameters, input/resource identities and execution metadata |

Repeat the same `run` command to resume. Successful caches are checked against raw output hashes and reparsed; changed inputs, code identities, dependencies or parameters prevent inappropriate reuse. Changed manifests require a new output directory. BAM identities use three 64 KiB samples and file metadata; they do not certify whole-file integrity.

## Changelog

| Version | Main changes |
|---|---|
| 0.5.0.dev1 | Retry AA jobs that exceeded a unique-region budget using a separate bounded policy and verified results from the previous global run. |
| 0.5.0.dev0 | Discover seeds with global CNVkit, reuse verified stages, pair loci and reconstruct each Suite’s own native seeds. |
| 0.4.4 (development) | Add an official Suite parameter policy for selected-region reconstruction and classification. |
| 0.4.3 | Time bounded stages, use CYCLES for reconstruction, cache AA and AC separately and share verified identical AA jobs across Suite recipes. |
| 0.4.2 | Use the sample’s original CNVkit version as the local baseline. Report seed retention, gains, losses and interval overlap. Bundle CNVkit 0.9.14 in the Linux image. |
| 0.4.1 | Add per-sample result summaries for seed retention, AA and AC interval overlap, and label agreement. Reuse verified background statistics for eligible older AA versions. |
| 0.4.0 | Import batches from a TSV sample sheet. Add an index of sample reports, cohort compatibility plots and progress tracking. |
| 0.3.2 | Adopt the MIT license and enable regenerating reports from verified 0.3.1 results. |
| 0.3.1 | Rename the project to AmpliconSuite Cross-Version Concordance and introduce the `asconcordance` command. |
| 0.3.0 | Add Suite 1.5.3 to the version panel. Add patient grouping, paired effect estimates and independent source calibration. |
| 0.2.0 | Introduce local CNVkit and fixed-seed AA/AC comparisons with BAM, plus comparisons of existing results without BAM. |
| 0.1.0 | Initial prototype for running several Suite versions and comparing their results with a fixed reference. |

0.1.0 and 0.4.1 were local development milestones.

## Validation status

The 0.4.3 public `with_bam` test used native AA and AC for all eleven Suite recipes on an artificial positive BAM. The FULL and CYCLES runs produced identical graph, cycles, summary, classification and feature BED files for every recipe. Reusing matching AA reduced the number of reconstructions from eleven to ten. [Validation record](docs/validation_status.md)

Local 0.4.2 regression passed all 311 tests. A separate integration run used native CNVkit 0.9.14 as the local standard for 0.9.10 and 0.9.13, with the same native AA 1.6.r1 seed filter on an artificial tumor and matched normal. See [validation status](docs/validation_status.md) for the four historical native AA controls.

The 0.4.0 batch update passed **257 local tests with one optional integration skipped**. Public batch import, serial toy-BAM execution, failure continuation, interruption/resume, patient grouping and exact cohort-plot bounds were checked with synthetic inputs. This is software validation, not a patient compatibility result.

The 0.3.0 regression exercises independent upstream calibration, patient-level aggregation, exploratory comparisons, the paired change summary and partial native results. Actual CNVkit 0.9.10 and 0.9.13 completed synthetic BAM tests. Suite 1.5.3 has an audited source lock and local component launch evidence.

GitHub Actions passed the Python 3.10 and 3.12 synthetic regressions and the Linux amd64 0.3.0 environment checks. Real-patient bounded AA reconstruction and patient-level compatibility remain unverified. [Validation details](docs/validation_status.md) · [Setup and Docker image](https://github.com/YenningC/AmpliconSuite-Cross-Version-Concordance-Environments)

```sh
python -m unittest discover -s tests -v
```

Default tests create temporary synthetic inputs and do not download upstream tools or access patient data. Real CNVkit integration tests require explicitly configured external environments. Test configuration is in [`.github/workflows/tests.yml`](.github/workflows/tests.yml).

## Official comparison tools

Independent adapters for [EC-comparator](https://github.com/AmpliconSuite/EC-comparator) and AC's `results_comparator.py` are retained for optional structural diagnostics. `asconcordance` does not automatically invoke them or turn their scores into compatibility thresholds. See [tool evaluation and commands](docs/official_comparison_tools.md).

## Repository and licensing

`ascompare/` contains the implementation; `tests/` contains regression tests; `examples/` contains templates and synthetic fixtures; `scripts/` contains runtime preparation and packaging utilities. Some internal modules are retained for regression coverage, but the old `ascompare run/check` genome-wide workflow is disabled. Use `asconcordance` for this release.

This repository contains no patient data, private configurations, reference genomes, upstream algorithm environments or solver licenses. Keep actual analysis manifests, inputs and generated results outside version control.

The project's original code and documentation are available under the **MIT License**. Third-party components follow their own licenses. [MIT License](LICENSE) · [Repository upload instructions](GITHUB_UPLOAD.md).
