# Validation status

## 0.5.0.dev0 全局发现

人工肿瘤与正常 BAM 未提供原种子 BED。11 套锁定 Suite 全部完成全局筛种、原生 AA 和 AC。两个不同的 CNVkit 配方各运行一次，其余九套使用校验通过的 CNS。第二版的肿瘤及正常原始覆盖度查询均命中第一版的计数。在独立人工对照上，每一版的肿瘤覆盖度、正常覆盖度、参考、CNR 及 CNS 都与同版本直接执行的原生 batch 和 CBS 逐字节相同。另用旧、新两版完整官方 Suite 入口复核筛种，最终坐标、顺序和 CN 一致，BED 的来源路径随运行目录变化。单个合成样本的兼容结论为证据不足。[逐版验证结果](validation/global_wgs.synthetic.json)


## 0.4.4.dev2 empty candidate controls

A valid empty BED produces no-seed records for native Suite comparisons, fixed-source CNVkit comparisons and AA/AC comparisons, even when the component runtime is unavailable. The installed CLI produced a complete report for an explicitly artificial empty BED with unavailable AA. No runtime probes or scientific subprocesses ran. Descriptive expected and evaluable seed counts were both zero, and the compatibility table reported insufficient evidence. A separate regression kept a missing BED file as an input error.

## 0.4.4.dev1 coverage-statistics checks

The wrapper checks whether the locked AA reader can accept a verified seed-stage measurement. The paired-read support stored by the upstream stage is preserved. When native AA recomputes coverage, audited and structurally valid reconstruction can proceed to AC. Lack of a reuse marker is recorded as an unverified reuse outcome.

In an artificial BAM control, Suite 1.3.7, 1.4.0, 1.5.0 and 1.5.2 skipped staged measurements with an artificial stored support of 11, recomputed coverage, and finished native AA and AC within the access budget. Suite 1.3.3 reused a genuine compatible measurement from the artificial fixture and finished both stages. All five kept FULL/EXPLORE, 10× target depth and their native thresholds. The local suite passed 369 regressions with the optional native CNVkit case enabled. [Synthetic result](validation/aa_coverage_fallback.synthetic.json)

## 0.4.4 development Suite-default controls

The `suite_defaults` policy was checked against all eleven locked Suite sources. It retains the default FULL/EXPLORE/10× invocation, uses the locked Suite's license check for recent solver selection and lets the native classifier supply its numerical defaults. It treats coverage-statistic recomputation as a normal native outcome when AA succeeds. A run with an inherited non-fixed Python hash seed does its own AA reconstruction per Suite, with resume restricted to the same job.

The local regression suite passed 366 tests with the native CNVkit 0.9.13 integration enabled. It covers unchanged native BAM paths containing spaces, official classifier defaults in both workflows, license selection, read scope and AA background recomputation. Replaying an existing 0.4.3 synthetic report produced a byte-identical `result.tsv`.

The native CLI completed AA and AC for all eleven recipes on one artificial linear fixture, with local CNVkit. Captured commands matched FULL/EXPLORE/10×, left the Python hash seed unset in the test environment and used native AC numerical defaults. [Synthetic record](validation/official_suite_defaults.synthetic.json)

The selected-region workflow uses supplied seed regions and a finite case-control CNVkit design. Full-pipeline genome-wide discovery and BAM-wide QC rules are outside this validation scope.

## 0.4.3 AA execution and output controls

We used the public `with_bam` CLI on an artificial tumor/normal pair with a 110 kb amplified candidate. The control used the complete eleven-recipe panel and fixed reference Suite 1.6.2. For each recipe the released 0.4.2 workflow ran its own AA in FULL/EXPLORE and its own AC. The 0.4.3 workflow used CYCLES/EXPLORE and the same seeds, filters, thresholds and resource limits. Native CNVkit, AA and AC were executed locally. This control has one Linear feature.

- The AA breakpoint graphs, cycles, summaries, AC classification profiles and feature BED files matched byte for byte for all eleven recipes. The comparison retained all eleven Suite records. Suite 1.3.3 and 1.3.4 shared one verified AA 1.4.r2 reconstruction and ran their distinct AC versions. This reduced the eleven AA executions to ten on this control.
- The six Concordance report PNGs matched byte for byte. Five named result tables also matched byte for byte. The local CNVkit seed table preserved all cells after making the two provenance paths relative to their respective run directories. [Machine-readable comparison](validation/aa_optimization.synthetic.json)
- A second artificial BAM contained two amplicons. With native AA 1.4.r2 and 1.6.r2, graph/cycles files from FULL and CYCLES matched at both amplicons. Holding AC 2.0.1 fixed for this mode comparison, the classification profiles and BEDs matched for a cyclic ecDNA-positive amplicon and a Linear amplicon. Both AA revisions matched with downsampling disabled, with requested depth 10 (all reads retained on the low-depth fixture), and with requested depth 0.25 (roughly half the reads retained). [Two-amplicon measurements](validation/aa_two_amplicons.synthetic.json)
- A separate synthetic run lowered the unique-access cap to 1,000,000 bp. Both tested recipes (Suite 1.6.1 and 1.6.2) reported AA failure, recorded the rejected request, skipped AC and produced no successful AA cache entry. The configured 100 Mb patient-run limit and performance on patient structural variation require a separate run.

All 355 local regression tests passed with the optional native CNVkit 0.9.13 integration enabled. Tests include AC-only retry, changed parameters and inputs, corrupted artifacts, partial multi-amplicon output, missing dependency inventories, resumed timing and failed access budgets.

To reproduce the comparison, use two independent output directories for the same positive synthetic manifest. Run the baseline with the 0.4.2 CLI, the optimized run with 0.4.3, then run the following against the completed directories. The verifier checks the named scientific artifacts and report figures and returns an error on a difference. [Artificial BAM fixture instructions](local_cnv_calibration.md#synthetic-control-fixture)

```sh
python scripts/check_aa_optimization_control.py runs/full-0.4.2 runs/cycles-0.4.3
```

## 0.4.2 local CNVkit version controls

The source-matched comparison was exercised with an artificial tumor and matched normal, one 110 kb candidate and a measured background of 800 kb. CNVkit 0.9.14 supplied the local standard. Native CNVkit 0.9.10, 0.9.13 and 0.9.14 ran with the same source Suite 1.6.1 / AA 1.6.r1 seed filter. Each version retained the full candidate. Both cross-version comparisons yielded interval Jaccard 1 and retention difference 0. The original-candidate calibration remained unverified because the fixture had no original-run identity evidence. This validates the positive software path on one artificial case. [Synthetic control results](validation/local_cnv_versions.synthetic.json)

The public `with_bam` CLI was also exercised with both Suite 1.6.1 and 1.6.2 using the locked 0.9.13 recipe. The original 0.9.14 environment was registered independently, and 0.9.10 was an extra comparison. The program reused the measured native source job for 0.9.13, ran the independent standard and extra version, and produced sample and cohort figures. AA/AC completed for both selected Suite versions.

Regression cases also cover missing version environments, a declared version conflicting with the original execution record, partially evaluable panels, empty intervals in both versions, patient grouping and integrity checks on replay.


## 0.4.1 similarity and bounded AA update

The new tests cover local final seed overlap, retained candidate counts, complete versus missing classes/coordinates, patient grouping, result-only replay, and AA interval evidence independent of AC features. The CNV control script checks existing CNR/CNS tables on mutually observed bases.

Four locked historical Suite entries (1.2.2, 1.3.3, 1.3.4 and 1.3.7, with three distinct AA source revisions) were run on **synthetic BAM**. Each entry ran its native seed filter and two pairs of AA attempts (background only, and a nonempty fixed seed with FULL/EXPLORE). On the nonempty panel with active downsampling and a fixed Python hash seed, reusing native coverage statistics saved approximately 20–21% of unique queried bases; graph and cycles files were byte-identical within each cold/reuse pair. These figures describe the artificial fixture. Historical AA patient completion under the 100 Mb cap remains to be measured. [Recorded native measurements](validation/legacy_aa_background.synthetic.json)

## 0.4.0 batch and cohort update

Local regression discovered **258 tests: 257 passed, one optional external-CNVkit integration skipped**. The new public TSV importer was exercised through `import-samples`, `check`, and `run` on three explicitly synthetic samples. A deliberately invalid middle sample remains recorded while the last sample is analyzed. Interruptions, completed-result preservation and manifest identity rejection were tested through the public CLI.

A separate two-sample with-BAM test executes guarded toy AA/AC subprocesses on temporary synthetic BAMs. It verifies serial sample/version order, paired fixed seeds, unchanged original toy inputs, and one independent patient for two linked samples. Component identity verification and CNVkit in that test are explicit synthetic stubs; it does not validate official AA/AC on patient data.

Per-sample figures and the new cohort compatibility figure were generated in PNG/PDF/SVG and visually checked. Plot-table tests check exact existing decision limits, margins, patient counts, source/analysis strata, shared-class annotations, and separation of missing/failed/all-negative observations. The scientific acceptance formulas are unchanged; the existing lower decision bound is now also exported.

The update also includes the native AC summary-map writer with synthetic subprocess and parser tests. Local tests do not establish real-patient cross-version compatibility or successful full reconstruction across all eleven official runtimes. Current-release hosted CI and environment publishing are tracked separately in GitHub Actions; prior release evidence below remains historical.

## 0.3.0 baseline evidence

The fixed reference is Suite 1.6.2. The 0.3.0 update was developed and tested with synthetic inputs. The local regression suite checks patient clustering, exploratory cohorts, signed paired changes and their direction, an independent source calibration, the reachable upstream exceedance branch, partial amplicons and format-only preflight. Real CNVkit 0.9.10 and 0.9.13 were separately exercised on temporary synthetic BAMs. Suite 1.5.3 source identities were checked against the official image and source locks. Its AA/AC launch and a small local Mosek solver probe passed.

GitHub Actions passed the 0.3.0 synthetic regression on Python 3.10 and 3.12, with 231 tests passed and one optional integration case skipped on each version. The Linux amd64 build separately tested both pinned CNVkit versions on synthetic inputs, reference AA launch, and AC launches across the eleven-version panel. The public 0.3.0 image uses the pinned source commit `61a720bcf4de25ac5103961cbd80dd4bef315b34`. Real-patient bounded AA and cohort-level compatibility remain unverified. [Python regression](https://github.com/YenningC/AmpliconSuite-Cross-Version-Concordance/actions/runs/44041256493) · [Linux environment test](https://github.com/YenningC/AmpliconSuite-Cross-Version-Concordance-Environments/actions/runs/44041555708) · [Image publication](https://github.com/YenningC/AmpliconSuite-Cross-Version-Concordance-Environments/actions/runs/44042036496)

Machine-readable results are in [validation_summary.json](../validation_summary.json).

## 0.2.0 baseline evidence

The fixed comparison reference was Suite 1.6.2 / AA 1.6.r2 / AC 2.0.1. The initial validation contains **219 discovered regression tests: 218 passed, one optional actual-CNVkit test skipped**. The optional workflow was also run separately with explicitly configured real components. The current machine-readable file retains this baseline evidence.

| Area | Evidence | Limit |
|---|---|---|
| Input/import handling | Synthetic native-format files; coordinate/reference checks, source declarations, read-only input and directory-overlap tests | Version declarations do not prove the original execution environment; sampled BAM identities are not full-file checksums |
| Results-only workflow | Three synthetic samples and two supplied version outputs through the public CLI | Invented AA/AC outputs; no biological conclusions |
| AC-only orchestration/cache | Synthetic subprocesses, fixed-source-graph checks, raw-output reparsing and cache invalidation | Does not establish classification accuracy on real samples |
| Real AC/BFB dependencies | AC2.0.0/BFB1.0.1 and AC2.0.1/BFB1.0.2 dependency and small solver probes | Small solver problems do not validate complete classification workloads; backend selection is environment-dependent |
| Real local CNVkit | 0.9.10, 0.9.13 and 0.9.14: finite synthetic coverage/reference/fix/CBS; 0.9.13: additional flat-reference test | No real-patient seed-retention validation; 0.9.14 is a controlled condition |
| Bounded BAM access | Temporary synthetic indexed BAM, finite fetch/bedcov, connected-region queries, budgets, rejection of unbounded/whole-contig access | Native bedcov read counts are not observable; finite BED and time/query budgets apply |
| AA/AC scheduling | Synthetic subprocess success/failure/timeout, multi-seed projection, linked intervals, source-result bridging and resume | **Real-patient AA reconstruction under the new bounded executor has not been completed** |
| Statistics/reporting | Paired/sample-cluster behavior, negative exclusion, incomplete categories, missing seeds, insufficient evidence and calibration/bridge gates | Software tests do not establish biological truth or adequate statistical power |
| Packaging/install | Wheel imported from a new virtual environment outside the source directory; eight installation/CLI checks | Scientific dependencies were explicitly reused; not a clean Linux-server deployment test |

The real CNVkit tests used a finite synthetic BED of approximately32kb. Matched-normal cases produced seven guarded processes including two identity probes; flat produced six. Passing attempts recorded no blocked operations. CNVkit0.9.10 used explicit compatibility dependencies for CBS; HMM was not validated. See [environment setup](environment_setup.md).

## 0.2.0 GitHub distribution checks

The GitHub distribution preserves the scientific Python implementation. Public-facing README files, packaging utilities, repository exclusions and CI configuration were added; internal operational history and real-study records were excluded. Source-lock descriptive notes were generalized without changing component versions, commits, archive hashes or source-manifest hashes. A timeout regression test has a more generous child-process startup allowance for shared CI hosts; product timeouts are unchanged.

The default test suite uses temporary synthetic inputs and does not fetch upstream code or access patient files. The current regression workflow is configured for Linux/Python 3.10 and 3.12. Native Windows support has not been verified.

## Still required for research use

Provide traceable source outputs/version declarations, matching resources and actual locked runtimes. With BAM, also provide the final seed BED, existing index, explicit local background and the inputs for the selected normal/flat strategy. Historical AA can require a Mosek license. Component environments, reference genomes and licenses are not distributed.

Validate real bounded AA reconstruction, native historical seed filtering and relevant complex amplification structures on the intended platform. Local upstream calibration needs original provenance, and source-Suite seed reruns must be checked against original outputs. Missing historical evidence remains unverified; it cannot be manufactured retrospectively.

Predefine acceptable differences and independent sample units before formal acceptance testing. If parameters are trained, use separate validation data. The current release does **not** establish patient-level compatibility, genome-wide sensitivity, candidate discovery outside the supplied set, biological truth or complete Suite interchangeability.
