# Scientific baseline

Baseline recorded on 2026-09-14 from the requested checkout:

- Repository: `OmniDock`
- Commit: `8d4434d3a33eb83b1e12cad82944b02c83270e47`
- Branch: `codex/scientific-integrity`
- Working directory: `C:\Users\salma\Documents\Codex\2026-09-14\omnidock-scientific-integrity\work\OmniDock`
- Reference instructions read: `docs/correctness-update.md`, `pyproject.toml`, `setup.py`, and `pytest.ini`.
- No source files or tests were changed for this baseline. Another concurrent worker created the untracked `audit/scientific-integrity-plan.md`; it was preserved.

## Environment

Commands were run from PowerShell 7.6.5 on 64-bit Windows (`Microsoft Windows NT 10.0.26200.0`, Python process architecture `X64`). Python is CPython 3.12.10 at:

```text
C:\Users\salma\AppData\Local\Programs\Python\Python312\python.exe
```

The process had no active `CONDA_PREFIX` or `VIRTUAL_ENV`. `python -m pip --version` reported pip 25.0.1. Git reported 2.55.0.windows.4. Git Bash is installed at `C:\Program Files\Git\bin\bash.exe` and reported GNU bash 5.3.15(1)-release when run outside the sandbox; `bash` is not on `PATH`.

Selected installed Python package metadata:

| Package | Version/status |
| --- | --- |
| numpy | 2.5.2 |
| pandas | 3.0.5 |
| matplotlib | 3.11.1 |
| seaborn | 0.13.2 |
| scipy | 1.18.1 |
| PyYAML | 6.0.3 |
| requests | 2.34.2 |
| scikit-learn | 1.9.0 |
| pytest | 9.1.1 |
| rdkit | 2026.3.6; top-level `rdkit` import works, but chemistry imports are blocked below |
| biopython, openpyxl, questionary, pdb-tools | not installed |
| build, setuptools, wheel | not installed/available to `importlib.metadata` |
| meeko, gemmi, plip, prolif, MDAnalysis, openbabel, vina | not installed |

The following executable identities were checked with `Get-Command`; each was unavailable on `PATH`: `vina`, `vina_split`, `autodock4`, `autogrid4`, `gnina`, `smina`, `obabel`, `babel`, `pythonsh`, `prepare_ligand4.py`, `prepare_receptor4.py`, `conda`, `mamba`, `micromamba`, `apptainer`, `singularity`, `docker`, `podman`, `sbatch`, `srun`, and `qsub`. `wsl.exe` exists, but `wsl.exe --status` and `wsl.exe --list --quiet` reported that Windows Subsystem for Linux is not installed. No remote job or benchmark was invoked.

Chemistry import probes produced:

```text
rdkit                         OK (top-level package only)
rdkit.Chem                    ImportError: DLL load failed while importing rdinchi: An Application Control policy has blocked this file.
rdkit.Chem.rdinchi            same Application Control block
rdkit.Chem.rdDetermineBonds   same Application Control block
meeko, gemmi, plip, prolif,
MDAnalysis, openbabel, vina   ModuleNotFoundError (not installed)
```

The RDKit result is an operating-environment restriction, not evidence that the source assertions passed or failed. Do not replace RDKit or bypass the Windows policy to make this baseline green.

## Requested baseline commands

The requested POSIX spelling of the test command is:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test -q
```

On this Windows host the equivalent PowerShell invocation was:

```powershell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD='1'; python -m pytest test -q
```

An in-sandbox run exited `1` after 19.47 seconds. Collection stopped with three import errors:

- `test/test_postdocking_correctness.py`: `rdkit.Chem` could not load `rdinchi` because Windows Application Control blocked its DLL.
- `test/test_preparation_correctness.py`: same `rdinchi` block.
- `test/test_reference_chemistry.py`: same `rdinchi` block.

That run also emitted two `PytestCacheWarning` messages because the sandbox could not create pytest cache temporary directories (`WinError 5`, access denied).

The same invocation was then run outside the sandbox, with the normal test permission escalation permitted for this baseline. It exited `1` after 19.29 seconds with the same three collection errors and no cache warnings. Thus the RDKit block persists independently of sandbox filesystem restrictions. No test assertions ran and this is not a passing correctness baseline.

The build command was run exactly as requested:

```text
python -m build
```

It exited `1` immediately with:

```text
C:\Users\salma\AppData\Local\Programs\Python\Python312\python.exe: No module named build
```

No wheel was produced. The requested wheel check was still run:

```text
python scripts/check_wheel.py
```

It exited `1` with `Build a wheel before running this check`. This is downstream of the missing build tool and does not establish a wheel defect.

## Shell syntax

`bash -n` could not start inside the sandbox because Git Bash itself failed with `CreateFileMapping ... Win32 error 5`. The syntax checks were repeated outside the sandbox using `C:\Program Files\Git\bin\bash.exe -n` for each tracked Bash script:

| Script | `bash -n` result |
| --- | --- |
| `autodock/prep_autodock.sh` | exit `2`; syntax error at line 77 near `2` in `for mol in "$RAW_LIG"/*.{mol2,sdf} 2>/dev/null; do` |
| `autodock/prep_ligands_custom.sh` | exit `2`; syntax error at line 16 near `2` in the same loop form |
| `install_vina_dock_dependencies.sh` | exit `0` |
| `prep_autodock_enhanced.sh` | exit `0` |

The two parse errors are source-level baseline findings. Each affected file contains a later `.mol` loop with the same redirection form, but Bash stops at the first parse error; no claim is made here about later diagnostics.

## Reproduction path for a usable chemistry baseline

Use a policy-approved Python environment with the declared chemistry and development extras, for example:

```bash
python -m pip install -e '.[chemistry,dev]'
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test -q
python -m build
python scripts/check_wheel.py
```

The host must also provide the external docking/preparation binaries separately and record their executable paths and version output. If Windows Application Control still blocks `rdinchi` or `rdDetermineBonds`, rerun on an approved host/environment where the installed RDKit DLLs are allowed (such as the configured CI/Linux environment). The current machine has no WSL or Conda environment, so shell execution and real engine identity checks remain unavailable here. The initial baseline did not install system-level packages, bypass policy, invoke remote benchmarks, or modify source/tests.

## Bounded verification venv follow-up — 2026-09-14

To separate ordinary packaging/runtime gaps from the blocked chemistry backend, a sibling virtual environment was created at:

```text
C:\Users\salma\Documents\Codex\2026-09-14\omnidock-scientific-integrity\work\verification-venv
```

It was created with `python -m venv --system-site-packages` outside the sandbox and retains the system RDKit installation. Only these missing ordinary packages were installed into that venv: `build==1.6.1`, `setuptools==84.0.0`, `wheel==0.48.0`, `biopython==1.88`, `openpyxl==3.1.5`, `questionary==2.1.1`, and `pdb-tools==2.7.0`. No alternate RDKit, Meeko, Gemmi, Open Babel, or docking engine was installed.

The filtered nonchemistry command was:

```powershell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD='1'; & verification-venv\Scripts\python.exe -m pytest test -q `
  --ignore=test/test_postdocking_correctness.py `
  --ignore=test/test_preparation_correctness.py `
  --ignore=test/test_reference_chemistry.py
```

The three ignored modules are the files that fail collection on the blocked RDKit `rdinchi` import in the full baseline. A sandbox run of the filtered command was unusable (`75 failed, 26 passed, 46 errors`) because pytest and project fixtures could not create temporary directories (`WinError 5`). The same filtered command was run once outside the sandbox and completed in 240.45 seconds with:

```text
2 failed, 145 passed, 570 warnings
```

The two failures were `test/test_dockforge_smoke.py::_smoke_redocking_multi_pose_sdf_parser_contract` and `test/test_dockforge_smoke.py::_smoke_redocking_validation_multi_pose_end_to_end_contract`. Both report `rdkit_required_for_chemical_mapping`; they remain chemistry/backend-limited even though the three RDKit collection modules were excluded. The other 145 collected tests passed in this environment. This filtered run is therefore a partial gate, not a full scientific-correctness pass.

Using the same venv interpreter, the build command was rerun outside the sandbox:

```text
verification-venv\Scripts\python.exe -m build
```

It exited `0` and built `pdb_prepare_wizard-3.0.1.tar.gz` and `pdb_prepare_wizard-3.0.1-py3-none-any.whl`. The earlier sandbox attempt with this interpreter failed only while creating `build`'s temporary isolated environment under `%TEMP%` (`WinError 5`). The wheel checker was then run outside the sandbox:

```text
verification-venv\Scripts\python.exe scripts/check_wheel.py
```

It exited `0` and reported `Imported every console entrypoint from the extracted wheel outside the checkout` and `Verified runtime files and console entrypoints: pdb_prepare_wizard-3.0.1-py3-none-any.whl`.

The venv is outside the repository and remains untracked. Build outputs under `build/`, `dist/`, and `pdb_prepare_wizard.egg-info/` are ignored by the repository. These packaging results are independent of the RDKit/Application Control limitation and do not certify chemistry execution or real docking engines.
