# Version identity and panel scope

The reference is frozen at Suite **1.6.2 / AA 1.6.r2 / AC 2.0.1 / BFBArchitect 1.0.2**, selected on 2026-11-17. This is a reproducible benchmark reference, not a promise to follow future releases automatically.

The eleven Suite recipes and component versions are listed in the root README. The exact source repositories, tags, 40-character commits, archive SHA-256 values, source-manifest SHA-256 values and evidence URLs are recorded in `versions.lock.json`. Its packaged copy is `ascompare/versions.lock.json`.

## What the lock does and does not establish

Source preparation verifies the selected public archives and source inventories. Execution additionally verifies real versions, import locations, interpreters and dependencies. A directory name, image tag or imported result declaration alone is not execution evidence. Missing historical runtimes remain missing; the program never substitutes the current AA/AC silently.

Suite **1.3.4 / AA 1.4.r2 / AC 1.2.1** is an explicitly selected benchmark recipe supported by recorded component pairing. The Suite release alone does not pin that exact combination, so it must not be inferred for every independently installed Suite 1.3.4 environment. Input declarations that disagree with the selected component recipe cannot be treated as matching by Suite label alone.

The lock also retains a separate legacy diagnostic combination, **AA 1.3.r2 + AC 1.2.1**, for existing preparation/regression interfaces. This recipe has no Suite identity and sits outside the eleven-Suite `asbench` panel. The original classifier of a publication must be established separately.

## CNVkit mapping

`ascompare/builtin_cnvkit_map.json` records the Docker/Singularity requirements at the exact locked Suite commits. The first seven Suite recipes pin CNVkit **0.9.10**. Suite 1.5.3 and the three 1.6.x releases pin **0.9.13**. These pins describe the container build recipes. Conda or independent installation constraints can lead to other versions, which must be reported as actually installed.

CNVkit 0.9.14 is an explicitly identified controlled diagnostic condition. Shared dependencies and runtime adaptations are recorded to distinguish a working environment from an exact historical container recreation.

### Evidence for Suite 1.5.3

The published `jluebeck/ampliconsuite-pipeline:v1.5.3` manifest resolves to `sha256:5442650a0e2c8ed036a2c4f29d3ceba564c20383d4467b5dcffeb003922c71eb`. The hashed image source layers contain GitHub ZIP archives with the Suite commit `3820e75013109c59147383cae49a1ee26ad51428`, AA commit `cc84a765c688c9d2d5a77f58afd89a6161dae038` and AC commit `b8a0abf1d4ef144ac8ba333052b376db503374e6`. Those commits match the Suite 1.5.3, AA 1.5.r7 and AC 1.5.2 release tags. This method establishes the actual source pairing even though the Dockerfile fetched AA and AC branch heads. Both published build requirement files pin CNVkit 0.9.13 on line 6. [Suite 1.5.3 release](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/releases/tag/v1.5.3) · [Published image tag](https://hub.docker.com/r/jluebeck/ampliconsuite-pipeline/tags?name=v1.5.3)

AA 1.5.r7 uses the historical Mosek solver. AC 1.5.2 has the ecDNA and BFB fields and uses the same supported classes as Suite 1.5.2 for this panel. The recipe predates BFBArchitect integration and FAN output. Strict source acquisition, AA/AC version and help probes, and a minimal Mosek optimization completed in the local runtime. The published container's installed Python dependencies and real sample execution remain unverified.

## Important algorithm generations

- [Suite 1.4.0](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/releases/tag/v1.4.0) covers the AA1.5/AC1.4 generation.
- [AC 1.5.1](https://github.com/AmpliconSuite/AmpliconClassifier/releases/tag/v1.5.1) includes coordinate/context fixes relevant to result comparison.
- [Suite 1.6.0](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/releases/tag/v1.6.0) introduces the coordinated Clarabel/FAN/BFBArchitect generation.
- [Suite 1.6.1](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/releases/tag/v1.6.1) includes stable downsampling changes. The benchmark defaults to disabled downsampling and therefore does not assess those changes under default settings.
- [Suite 1.6.2](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/releases/tag/v1.6.2) is the fixed current benchmark recipe.

### Controlled check for Suite 1.6.1 downsampling

Suite 1.6.0's AA uses Python `hash()` for read-name selection. Suite 1.6.1 switches to `crc32()`, which addresses variation across processes. Compare Suite 1.6.0 with 1.6.1 under the same input, reference, selected seeds, solver backend and AC version. Configure both 1.6.0 and 1.6.1 as targets. The fixed global 1.6.2 reference remains in the report, with the historical pair available in the pairwise table. [AA 1.6.r1 release](https://github.com/AmpliconSuite/AmpliconArchitect/releases/tag/v1.6.r1)

Start with the default `aa_parameters.downsample: -1` as a control. Use a separate set of runs with `aa_parameters.downsample` set to an explicit positive coverage target below measured coverage, for example `10` for eligible inputs above 10×. Confirm in the AA log that the resulting ratio is below one. Keep `sampling.random_seed` fixed to select the same candidate regions. Use a fresh output directory and independent AA process for each repeat so cache hits do not count as repeats. Keep other runtime settings fixed and record `PYTHONHASHSEED`. For a cross-process stress test, deliberately vary it in the declared runtime environment across repeated pairs, with the same values for both versions. Compare within-version reproducibility using scientific AA/AC fields after excluding timestamps and other volatile metadata. Repeats measure reproducibility on the same input and do not add independent patients.

This is a selected historical panel, not every published research combination. Source availability, version/help probes and solver smoke tests are distinct from real-sample completion. Older AA needs a valid compatible Mosek setup; modern BFB solver selection can also depend on visible licenses and runtime dependencies. Capture actual backend evidence rather than assuming an installed solver was used.

The bounded workflow does not validate genome-wide discovery, every reference build, viral sequences, all structural classes or all platforms. See [validation status](validation_status.md) and [environment setup](environment_setup.md).
