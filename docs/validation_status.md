# Validation status — 0.3.0

The fixed reference remains Suite 1.6.2. The 0.3.0 update was developed and tested with synthetic inputs. The local regression suite checks patient clustering, exploratory cohorts, signed paired changes and their direction, an independent source calibration, the reachable upstream exceedance branch, partial amplicons and format-only preflight. Real CNVkit 0.9.10 and 0.9.13 were separately exercised on temporary synthetic BAMs. Suite 1.5.3 source identities were checked against the official image and source locks. Its AA/AC launch and a small local Mosek solver probe passed.

GitHub Actions passed the 0.3.0 synthetic regression on Python 3.10 and 3.12, with 231 tests passed and one optional integration case skipped on each version. The Linux amd64 build separately tested both pinned CNVkit versions on synthetic inputs, reference AA launch, and AC launches across the eleven-version panel. The public 0.3.0 image uses the pinned source commit `61a720bcf4de25ac5103961cbd80dd4bef315b34`. Real-patient bounded AA and cohort-level compatibility remain unverified. [Python regression](https://github.com/YenningC/AmpliconSuite-Benchmarks/actions/runs/44041256493) · [Linux environment test](https://github.com/YenningC/AmpliconSuite-Benchmarks-Environments/actions/runs/44041555708) · [Image publication](https://github.com/YenningC/AmpliconSuite-Benchmarks-Environments/actions/runs/44042036496) · [Public environment package](https://github.com/YenningC/AmpliconSuite-Benchmarks-Environments/pkgs/container/ampliconsuite-benchmarks-environments)

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

The default test suite uses temporary synthetic inputs and does not fetch upstream code or access patient files. The current Benchmarks workflow is configured for Linux/Python 3.10 and 3.12. Native Windows support has not been verified.

## Still required for research use

Provide traceable source outputs/version declarations, matching resources and actual locked runtimes. With BAM, also provide the final seed BED, existing index, explicit local background and the inputs for the selected normal/flat strategy. Historical AA can require a Mosek license. Component environments, reference genomes and licenses are not distributed.

Validate real bounded AA reconstruction, native historical seed filtering and relevant complex amplification structures on the intended platform. Local upstream calibration needs original provenance, and source-Suite seed reruns must be checked against original outputs. Missing historical evidence remains unverified; it cannot be manufactured retrospectively.

Predefine acceptable differences and independent sample units before formal acceptance testing. If parameters are trained, use separate validation data. The current release does **not** establish patient-level compatibility, genome-wide sensitivity, candidate discovery outside the supplied set, biological truth or complete Suite interchangeability.
