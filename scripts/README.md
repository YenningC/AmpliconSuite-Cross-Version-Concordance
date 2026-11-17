# Auxiliary scripts

Run these from the repository root after installing the comparison package. Inputs remain local.

| Script | Purpose |
|---|---|
| `prepare_panel.py` | Fetch and verify the locked public source snapshots; does not process samples |
| `build_concordance_manifest.py` | Combine existing-result input with prepared sources and explicit runtime templates for `asconcordance` |
| `create_runtime_envs.py` | Create explicit AA/AC environments from a scientific base environment; its historical CNVkit0.9.14 template is not the native CNVkit panel |
| `prepare_cnvkit_smoke_envs.py` | Prepare isolated exact-version CNVkit test environments |
| `run_cnvkit_smoke.py` | Run the explicitly configured finite synthetic CNVkit integration test |
| `probe_cnvkit_record.py`, `probe_clarabel.py`, `probe_mosek.py` | Component identity/solver diagnostics |
| `adapt_native_cbc.py` | Explicit, recorded native CBC runtime adaptation when required; not an algorithm replacement |
| `package_program.py` | Build this public source ZIP/TAR.GZ and optional wheel, including bilingual docs and hidden GitHub configuration |
| `build_panel_config.py`, `provision_bash.py` | Retained legacy preparation utilities used by regression tests or historical environments; not the current sample workflow |

For sample analysis, use `asconcordance check/run/report`. The old `ascompare run/check` workflow is disabled. Do not use a legacy generated configuration to bypass the bounded workflow.

Packaging example; use a new directory outside the repository:

```sh
python -m pip install 'setuptools>=77.0.3' wheel
python scripts/package_program.py --out ../concordance-release --wheel
```

The packaging allowlist does not constitute an automatic sensitive-data detector. Keep actual inputs, manifests and outputs outside the source tree. Both the source archive and the wheel include the repository’s MIT license.
