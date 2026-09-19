# Independent scientific verification

Recorded 2026-09-14 by an independent test worker. The pinned baseline is
`8d4434d3a33eb83b1e12cad82944b02c83270e47` in
`work/OmniDock`. The shared checkout contains concurrent scientific changes;
this report owns no implementation files and records only verification results.

## Verification environment

Tests and packaging use the sibling interpreter:

```text
C:\Users\salma\Documents\Codex\2026-09-14\omnidock-scientific-integrity\work\verification-venv\Scripts\python.exe
```

The venv was created with `--system-site-packages` and contains only the
ordinary dependencies needed for this pass (`build==1.6.1`,
`setuptools==84.0.0`, `wheel==0.48.0`, `biopython==1.88`, `openpyxl==3.1.5`,
`questionary==2.1.1`, and `pdb-tools==2.7.0`). The system RDKit was not
replaced. `rdkit.Chem` remains unavailable because Windows Application Control
blocks its `rdinchi` DLL; the chemistry-only pose and ligand tests are therefore
excluded from the filtered run below.

## Baseline replay

To keep the shared working tree unchanged, a source archive of the pinned commit
was extracted to `work/verification-archive`, and the seven new
`test_scientific_*.py` files were copied into its `test/` directory. The archive
contains the baseline implementation; no current implementation files were
copied into it. One historical absolute symlink under `4-K_pneumoniae_` was
excluded during Windows extraction because it cannot be materialized safely; it
is unrelated to these tests.

The five new nonchemistry contract files were run individually from the archive
outside the sandbox with pytest temporary-directory access:

```powershell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD='1'; ..\verification-venv\Scripts\python.exe -m pytest test\<file> -q
```

The exact per-file results were:

| Test file | Result on pinned baseline |
| --- | --- |
| `test_scientific_consensus_contract.py` | `1 passed, 8 failed` |
| `test_scientific_pose_records.py` | `3 passed, 10 failed` |
| `test_scientific_receptor_contract.py` | `2 passed, 13 failed` |
| `test_scientific_reference_contract.py` | `0 passed, 1 failed` |
| `test_scientific_preparation_config.py` | `9 passed, 5 failed` |

The aggregate baseline replay was `15 passed, 37 failed`; all 52 tests
collected. There were no missing-symbol or import-collection failures in these
five files. The failures are behavioral evidence against the pinned code:

- Consensus accepted duplicate engine/tag rows, inconsistent tag identity,
  mixed scoring functions, duplicate or out-of-scope expected engines, and
  malformed expected-engine values.
- Pose-record parsing accepted empty records or non-standalone separators and
  did not fail closed for malformed model boundaries.
- Receptor selection/conservation did not enforce crossing connectivity records,
  source-scoped frame IDs, atom maps/identity/coordinates, or early invalid
  input checks and explicit PDB2PQR titration options.
- Reference validation accepted duplicate reference input rows.
- Preparation configuration accepted nonfinite/boolean pH and silently
  defaulted an unknown profile.

These are original baseline failures, not environment exclusions. The two
chemistry-aware files (`test_scientific_pose_chemistry.py` and
`test_scientific_ligand_contract.py`) were copied into the archive but were not
run because their top-level RDKit imports are blocked on this host.

## Current frozen nonchemistry contracts

The currently frozen consensus, pose-record, receptor, reference, and
preparation-config files were then run from the shared checkout with the same
venv interpreter and environment:

```powershell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD='1'; ..\verification-venv\Scripts\python.exe -m pytest `
  test/test_scientific_consensus_contract.py `
  test/test_scientific_pose_records.py `
  test/test_scientific_receptor_contract.py `
  test/test_scientific_reference_contract.py `
  test/test_scientific_preparation_config.py -q
```

Result: `52 passed in 2.32s`, exit code `0`, with no skips. This is the earlier
targeted checkpoint, recorded before the final preparation-config guard was
added; the final current filtered regression below contains 53 tests.

## Final frozen current gates

The ligand worker formally handed off the frozen implementation before these
commands. The strict full collection command was:

```powershell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD='1'; ..\verification-venv\Scripts\python.exe -m pytest test -q
```

It exited `2` after 9.26 seconds with exactly five collection errors, all caused
by the host's Windows Application Control policy blocking RDKit's
`rdCIPLabeler` DLL. The excluded modules are:

- `test/test_postdocking_correctness.py`
- `test/test_preparation_correctness.py`
- `test/test_reference_chemistry.py`
- `test/test_scientific_ligand_contract.py`
- `test/test_scientific_pose_chemistry.py`

No test assertions ran in those five modules. The collection result is an
environment limitation, not evidence that their assertions pass or fail.

The explicit local filtered regression command was:

```powershell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD='1'; ..\verification-venv\Scripts\python.exe -m pytest test -q `
  --ignore=test/test_postdocking_correctness.py `
  --ignore=test/test_preparation_correctness.py `
  --ignore=test/test_reference_chemistry.py `
  --ignore=test/test_scientific_pose_chemistry.py `
  --ignore=test/test_scientific_ligand_contract.py
```

It exited `1` after 272.67 seconds with `195 passed, 5 failed, 546 warnings`.
The two redocking smoke failures remain chemistry-limited:

- `test/test_dockforge_smoke.py::_smoke_redocking_multi_pose_sdf_parser_contract`
- `test/test_dockforge_smoke.py::_smoke_redocking_validation_multi_pose_end_to_end_contract`

Both report `rdkit_required_for_chemical_mapping`. The three unexpected
nonchemistry failures are:

- `_smoke_sqlite_parity_and_skip_rationale`
- `_smoke_post_docking_allscore_contracts`
- `_smoke_biology_unresolved_mapping_reporting`

Each reaches `build_consensus_rankings` with `expected_engines=[]`, which now
raises `ValueError("empty_expected_engine_scope")`. These are implementation
regressions requiring review; they are not Windows temp-directory failures.

With the same frozen sibling interpreter, the current build command was:

```text
..\verification-venv\Scripts\python.exe -m build
```

It exited `0` and built `pdb_prepare_wizard-3.0.1.tar.gz` and
`pdb_prepare_wizard-3.0.1-py3-none-any.whl`. The build emitted the existing
setuptools license-classifier deprecation and `post_docking_analysis.config`
package-discovery warnings. The current wheel check was:

```text
..\verification-venv\Scripts\python.exe scripts/check_wheel.py
```

It exited `0`, importing every console entrypoint from the extracted wheel
outside the checkout.

The final shell syntax command used the installed Git Bash executable:

```powershell
$gitbash='C:\Program Files\Git\bin\bash.exe'; $files=@('autodock/prep_autodock.sh','autodock/prep_ligands_custom.sh','install_vina_dock_dependencies.sh','prep_autodock_enhanced.sh'); foreach ($f in $files) { & $gitbash -n $f; Write-Output ("$f exit=" + $LASTEXITCODE) }
```

The inherited results are unchanged: `autodock/prep_autodock.sh` exits `2`
(line 77), `autodock/prep_ligands_custom.sh` exits `2` (line 16), and
`install_vina_dock_dependencies.sh` plus `prep_autodock_enhanced.sh` exit `0`.
These two parse errors belong to the existing shell cleanup scope.

Both `git diff --check` and `git diff --cached --check` exited `0` after the
verification run. No source or test file was modified by this worker.

After the final preparation-config guard and formal ligand handoff, the five
new nonchemistry contract files were rerun with the same sibling interpreter:

```powershell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD='1'; ..\verification-venv\Scripts\python.exe -m pytest `
  test/test_scientific_consensus_contract.py `
  test/test_scientific_pose_records.py `
  test/test_scientific_receptor_contract.py `
  test/test_scientific_reference_contract.py `
  test/test_scientific_preparation_config.py -q
```

The definitive focused result was `53 passed in 2.37s`, exit code `0`, with no
skips.

## Post-consensus correction rerun — 2026-09-14

The final reviewed source state was commit `48c0771` (`Preserve inferred
consensus scope for existing callers`). This rerun preserves the earlier
`195 passed, 5 failed, 546 warnings` discovery result above and measures the
same filtered regression after the `expected_engines=[]` compatibility fix.

The command was:

```powershell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD='1'; ..\verification-venv\Scripts\python.exe -m pytest test -q `
  --ignore=test/test_postdocking_correctness.py `
  --ignore=test/test_preparation_correctness.py `
  --ignore=test/test_reference_chemistry.py `
  --ignore=test/test_scientific_pose_chemistry.py `
  --ignore=test/test_scientific_ligand_contract.py
```

It exited `1` after 263.23 seconds with `199 passed, 2 failed, 570
warnings`. The only failures were the two previously known multi-pose
redocking smoke contracts:

- `test/test_dockforge_smoke.py::_smoke_redocking_multi_pose_sdf_parser_contract`
- `test/test_dockforge_smoke.py::_smoke_redocking_validation_multi_pose_end_to_end_contract`

Both remain host-chemistry limited: the parser reports
`rdkit_required_for_chemical_mapping`, and the end-to-end validation therefore
returns `not_evaluable`. The three earlier nonchemistry consensus failures no
longer reproduce after commit `48c0771`.

After that source state, the sibling interpreter build command
`..\verification-venv\Scripts\python.exe -m build` exited `0`, producing the
sdist and wheel. `..\verification-venv\Scripts\python.exe
scripts/check_wheel.py` also exited `0`. Existing setuptools license
classifier and package-discovery warnings were emitted during the build.

## Published CI verification — 2026-09-19

The reviewed branch was subsequently published at
`f5760bf45c69d42a4420a866c77aa8f6e74875cc`. The mandatory GitHub Actions
matrix then completed successfully on Ubuntu and Windows with Python 3.10 and
3.12. In every job, chemistry extras installed, full pytest reported `305
passed`, and the build and wheel checks passed. The authoritative job record is
[GitHub Actions run 34876658350](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/actions/runs/34876658350).

This supersedes no local result above: Windows Application Control still
prevented RDKit loading in the independent local environment, and the earlier
filtered failures and collection limits remain historical evidence for that
environment. The successful CI matrix verifies its declared full-suite and
packaging gates; it does not establish real-backend execution or benchmark
performance.

The published branch is represented by [draft PR #2](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/pull/2).
It has not been merged; remote `main` remains at the recorded baseline
`8d4434d3a33eb83b1e12cad82944b02c83270e47`.

## Historical final-head intermittent evidence — 2026-09-19

For final head `1cf6492`, push run
[35433642487](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/actions/runs/35433642487)
passed all four matrix jobs. The same-head PR run
[35433644917](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/actions/runs/35433644917)
had one failure and three passing jobs: Ubuntu Python 3.10 failed
`full_parallel_execution_contract` after `304 passed, 1 failed`.

An independent Luna replay of archived `1cf6492` on Windows Python 3.12.10
found four passing and one failing valid runs (02–06). In run 06,
`polypharmacology` failed with the exact detail `bad allocation`; reports
were blocked, with `9 completed, 1 failed, 1 blocked_by_failure`.
Run 01 is excluded because its wrapper hit a CP1252 print
error although its DAG succeeded. Captures are retained at
`work/repro-105872532281-captures/run06`. These observations are historical
evidence of an intermittent issue under investigation; no cause or fix is
claimed here.
