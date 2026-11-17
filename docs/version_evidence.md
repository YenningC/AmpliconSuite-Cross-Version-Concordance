# Version identity and panel scope

The reference is frozen at Suite **1.6.2 / AA 1.6.r2 / AC 2.0.1 / BFBArchitect 1.0.2**, selected on 2026-11-17. This is a reproducible benchmark reference, not a promise to follow future releases automatically.

The ten Suite recipes and component versions are listed in the root README. The exact source repositories, tags, 40-character commits, archive SHA-256 values, source-manifest SHA-256 values and evidence URLs are recorded in `versions.lock.json`. Its packaged copy is `ascompare/versions.lock.json`.

## What the lock does and does not establish

Source preparation verifies the selected public archives and source inventories. Execution additionally verifies real versions, import locations, interpreters and dependencies. A directory name, image tag or imported result declaration alone is not execution evidence. Missing historical runtimes remain missing; the program never substitutes the current AA/AC silently.

Suite **1.3.4 / AA 1.4.r2 / AC 1.2.1** is an explicitly selected benchmark recipe supported by recorded component pairing. The Suite release alone does not pin that exact combination, so it must not be inferred for every independently installed Suite 1.3.4 environment. Input declarations that disagree with the selected component recipe cannot be treated as matching by Suite label alone.

The lock also retains a separate legacy diagnostic combination, **AA 1.3.r2 + AC 1.2.1**, for existing preparation/regression interfaces. It has no Suite identity and is not part of the ten-Suite `asbench` panel. It is not claimed to reproduce the original classifier of a publication.

## CNVkit mapping

`ascompare/builtin_cnvkit_map.json` records the Docker/Singularity requirements at the exact locked Suite commits. The first seven Suite recipes pin CNVkit **0.9.10**; the three 1.6.x recipes pin **0.9.13**. These are build-recipe observations, not inspected published-image contents. Conda or independent installation constraints can lead to other versions, which must be reported as actually installed.

[Suite 1.5.3](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/releases/tag/v1.5.3) documents a CNVkit update but is not included in this panel. CNVkit0.9.14 is an explicitly identified controlled diagnostic condition, not a historical replacement. Shared dependencies and runtime adaptations are recorded; a working environment is not a recreation of every dependency in an old container.

## Important algorithm generations

- [Suite 1.4.0](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/releases/tag/v1.4.0) covers the AA1.5/AC1.4 generation.
- [AC 1.5.1](https://github.com/AmpliconSuite/AmpliconClassifier/releases/tag/v1.5.1) includes coordinate/context fixes relevant to result comparison.
- [Suite 1.6.0](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/releases/tag/v1.6.0) introduces the coordinated Clarabel/FAN/BFBArchitect generation.
- [Suite 1.6.1](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/releases/tag/v1.6.1) includes stable downsampling changes. The benchmark defaults to disabled downsampling and therefore does not assess those changes under default settings.
- [Suite 1.6.2](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/releases/tag/v1.6.2) is the fixed current benchmark recipe.

This is a selected historical panel, not every published research combination. Source availability, version/help probes and solver smoke tests are distinct from real-sample completion. Older AA needs a valid compatible Mosek setup; modern BFB solver selection can also depend on visible licenses and runtime dependencies. Capture actual backend evidence rather than assuming an installed solver was used.

The bounded workflow does not validate genome-wide discovery, every reference build, viral sequences, all structural classes or all platforms. See [validation status](validation_status.md) and [environment setup](environment_setup.md).
