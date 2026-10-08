# Software evidence — Linux clean-clone re-verification (2026-10-08)

**Evidence type:** software only. This record says nothing about scientific validity and does not reopen T032 acceptance.

**Tested commit:** `46e0f6a` (branch `chore/install-agentic-bioinformatics-stack`; base `44bfdb8` plus two test-infrastructure commits listed below). Working tree clean at test time.

## Environment

| Item | Value |
|---|---|
| OS / kernel | Ubuntu 24.04.5 LTS, Linux 6.18 (x86_64, ephemeral cloud container) |
| Environment manager | Miniforge, mamba 2.9.0, channel `conda-forge` only |
| Python | 3.10.22 (conda-forge) |
| Open Babel | 3.2.1 (`obabel -V`, `OBReleaseVersion()`); CLI from conda-forge, Python module from the PyPI 3.2.1 wheel installed with PLIP |
| numpy / pandas / scipy | 2.2.6 / 2.3.3 / 1.15.2 |
| biopython / plip / networkx | 1.88 / 3.0.1 / 3.4.2 |
| scikit-learn / matplotlib / seaborn | 1.7.2 / 3.10.9 / 0.13.2 |
| Flask / PyYAML / openpyxl | 3.1.3 / 6.0.3 / 3.1.5 |
| pdb-tools / questionary / pytest | 2.7.0 / 2.1.1 / 9.1.1 |
| RDKit | not installed (not a declared dependency) |

**Deviation from `environment.yml`:** Python 3.10 instead of 3.9. On linux-64, conda-forge publishes `openbabel=3.2.1` builds only for Python ≥ 3.10, so the 3.9 + 3.2.1 combination does not resolve. The Open Babel version pin was kept; the Python pin was relaxed for this run only. `environment.yml` is unchanged.

### Reproduction commands

```bash
# Miniforge installed to a scratch prefix, then:
mamba create -y -p <env> -c conda-forge --override-channels \
  python=3.10 openbabel=3.2.1 numpy pandas biopython matplotlib seaborn scipy \
  scikit-learn pyyaml openpyxl "networkx>=3.2,<4" pytest flask
<env>/bin/python -m pip install "plip>=2.2.0" "questionary>=2.0.1" "pdb-tools>=2.5.0"
<env>/bin/python -m pip install -r requirements.txt
<env>/bin/python -m pip install -e .
```

## Results

| Command | Result |
|---|---|
| `python -m pytest test/ -q -rs` at `44bfdb8` (fresh clone, before fixes) | 168 passed, 10 failed, 0 skipped |
| `python -m pytest test/ -q -rs` at `46e0f6a` | **177 passed, 1 skipped, 0 failed** |
| `python -m pytest test/test_spec032_md_inputs.py test/test_spec031_*.py -q` at `46e0f6a` | 69 passed, 1 skipped |
| `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix` at `46e0f6a` | exit 0 (2026-10-08 08:20–08:23 UTC) |
| `python -m compileall -q .` | exit 0 (syntax/bytecode only, not lint or type checking) |
| `git diff --check 44bfdb8 46e0f6a` | clean |

No lint or type-check command is configured in the repository, so none was run.

### Failures at `44bfdb8` and their causes (none scientific)

1. **9 × `test/test_webui.py`.** The synthetic web UI fixture `test/fixtures/webui_project/` had never been committed in full. The repository `.gitignore` rules `*.csv` and `.workflow/` excluded `consensus_ranked.csv`, `affinity_summary.csv`, `barcode_qc_rep1.csv` and `.workflow/state.json`, so these files existed only on the original development machine. **Fix `62e1f27`:** added `.gitignore` negations scoped to that fixture and committed minimal synthetic placeholder files (ligand names `ligA`–`ligC`, zero scores), built from the test assertions and the `workflow.state` default schema. No source code changed.
2. **1 × `test/test_spec031_remediation.py::test_vina_1iep_fixture_hashes_and_direct_parser_rows`.** The test reads the 1IEP artifacts from `docs/appraisal_20260819/`, which is intentionally gitignored and private, so a fresh clone does not have them. **Fix `46e0f6a`:** when any of the four recorded artifacts is absent, the test now skips with an explicit reason. When they are present, it runs unchanged, including the hash and parser assertions. The files were not copied or committed.

### Skips at `46e0f6a`

- `test_vina_1iep_fixture_hashes_and_direct_parser_rows`: skipped because the private 1IEP appraisal package was not present. **The 1IEP hash and parser regression was not executed in this environment.** The fixture hashes and the parsed Vina rows (−11.702, −8.898) are not re-established by this record.

The Spec 032 Open Babel-backed test was **not** skipped. It executed against Open Babel 3.2.1 on Linux.

## Not established by this record

- The scientific correctness of any method, score, pose, or protonation state.
- The 1IEP fixture regression in this environment (see the skip above).
- NMRbox evidence (T026 remains blocked).
- Behaviour under Python 3.9 on Linux.
