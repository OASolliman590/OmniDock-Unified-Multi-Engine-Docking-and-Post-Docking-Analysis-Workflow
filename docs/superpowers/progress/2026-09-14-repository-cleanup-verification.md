# OmniDock cleanup verification evidence

- **Date:** 2026-09-20
- **Fixed base:** `8d4434d3a33eb83b1e12cad82944b02c83270e47`
- **Verified implementation/review-fix head:** `1ca60902254c370beb71339b47a1fc34af39a439`
- **Host:** Microsoft Windows NT 10.0.26200.0, Python 3.12.10
**Working directory:** repository root unless a command says otherwise

This record replaces pass/fail shorthand with reproducible commands, exit
statuses, skips, warnings, and host limitations. Temporary directories are
shown as repository-relative paths so the record contains no personal absolute
paths.

## Environment

The pytest commands used the same relevant environment as CI:

```powershell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD = '1'
$env:PYTHONUTF8 = '1'
$env:MPLBACKEND = 'Agg'
```

## Repository and static gates

```powershell
python scripts/check_repository.py
$env:PYTHONPYCACHEPREFIX = '..\omnidock-compileall-final'
python -m compileall -q .
git diff --check 8d4434d3a33eb83b1e12cad82944b02c83270e47...HEAD
```

All three commands exited 0. `compileall` reported that it could not list stale
sandbox-owned `pytest-cache-files-*` directories; source compilation still
completed successfully. The working tree contained no tracked generated
artifacts or dangling symlinks.

Git Bash was run outside the managed sandbox because sandboxed Git Bash could
not create its Windows signal pipe:

```bash
bash -n install_vina_dock_dependencies.sh
bash -n prep_autodock_enhanced.sh
```

Both commands exited 0.

The expanded documentation audit and public help checks were:

```powershell
$audit = (python ..\docs_audit.py | ConvertFrom-Json)
python main.py --help
python main.py workflow --help
python main.py pdb --help
python main.py prep --help
python main.py dock --help
python main.py analyze --help
python -m post_docking_analysis --help
```

The audit checked 125 Markdown files and found zero broken relative links. All
seven help commands exited 0. Git reported its configured LF-to-CRLF working-copy
notices while staging Markdown; `git diff --check` itself exited 0.

## Focused cleanup regression

The final repository/deployment contract run was:

```powershell
python -m pytest test/test_repository_contracts.py test/test_docking_correctness.py `
  -q -p no:cacheprovider --basetemp ..\pytest-task6-core-final
```

Exit status 0: 44 passed, 1 skipped, 31 warnings in 14.88 seconds. The skip is
the native-symlink capability check on this Windows host. The warnings are the
expected legacy `5-Post_Docking_Analysis` compatibility notices emitted while
tests verify migration to `5-Analysis`.

The final parser/workflow and representative CLI/deployment smoke runs were:

```powershell
python -m pytest test/test_workflow_correctness.py `
  -q -p no:cacheprovider --basetemp ..\pytest-task6-workflow-final
python -m pytest test/test_dockforge_smoke.py `
  -q -k "deployment or hpc or condor or slurm or cli or workflow" `
  -p no:cacheprovider --basetemp ..\pytest-task6-smoke-final
```

Both exited 0. The workflow run reported 15 passed. The smoke run reported 14
passed, 86 deselected, and 19 expected legacy-layout compatibility warnings.
Earlier task-local reviews also compared all 759 recursive argparse actions and
found their option strings, destinations, defaults, required flags, choices,
types, and nested command structure identical across the parser extraction.

## Full regression and host limits

The full command was:

```powershell
python -m pytest test -q --basetemp ..\pytest-task6-full-final
```

Collection stopped with two errors because Windows Application Control blocked
RDKit's `rdForceField` DLL while importing
`test/test_postdocking_correctness.py` and
`test/test_preparation_correctness.py`:

```text
ImportError: DLL load failed while importing rdForceField:
An Application Control policy has blocked this file.
```

The non-RDKit command was then run outside the managed sandbox:

```powershell
python -m pytest test -q `
  --ignore=test/test_postdocking_correctness.py `
  --ignore=test/test_preparation_correctness.py `
  -p no:cacheprovider --basetemp ..\pytest-task6-nonrdkit-final
```

Exit status 1: 162 passed, 1 skipped, 1 failed, 572 warnings in 219.89 seconds.
The failure is
`_smoke_redocking_validation_multi_pose_end_to_end_contract`, where a synthetic
multi-pose SDF is classified `not_evaluable`. The failing smoke test and the
candidate redocking implementation paths are unchanged from the fixed base.
This cleanup did not change the scientific behavior, skip the test, or claim
that it passed. Clean Linux/Windows chemistry CI remains required.

## Distribution

The system interpreter lacked the `build` frontend, so it was installed into a
task-local scratch directory. The successful build command was run from the
parent directory to keep the frontend distinct from setuptools' local `build/`
output directory:

```powershell
$env:PYTHONPATH = '..\build-deps-final'
Push-Location ..
python -m build OmniDock
Pop-Location
```

Exit status 0: built `pdb_prepare_wizard-3.0.1.tar.gz` and
`pdb_prepare_wizard-3.0.1-py3-none-any.whl`. Setuptools emitted non-blocking
warnings about the legacy license classifier and the data-only
`post_docking_analysis.config` namespace; the required YAML files were present
in the wheel.

The wheel and all declared runtime dependencies were installed into a fresh
task-local virtual environment outside the source checkout. Verification used:

```powershell
..\wheel-venv-final\Scripts\python.exe scripts\check_wheel.py
```

Exit status 0: every console entrypoint imported from the extracted wheel, and
all required runtime files were present. Installed `--help` exited 0 for
`pdb-prepare-wizard`, `pdb-wizard-workflow`, `pdb-wizard-cli`,
`pdb-wizard-batch`, `pdb-wizard-prepare-docking`, `pdb-wizard-dock`, and
`post-docking-analysis`. The two interactive entrypoints were verified as
installed callables because they enter `questionary` immediately rather than
implementing `--help`; a non-console invocation produces the expected
`NoConsoleScreenBufferError`. The installed schema and preparation shell
resources were also present.

## Before and after

Worktree sizes are used at both endpoints so Windows line endings are compared
consistently.

| Measure | Fixed base | Verified head | Change |
| --- | ---: | ---: | ---: |
| Tracked files | 300 | 305 | +5 |
| Tracked bytes | 13,327,496 | 13,152,966 | -174,530 |
| Generated verification artifacts | 5 | 0 | -5 |
| Generated verification bytes | 184,736 | 0 | -184,736 |
| Broken relative Markdown links | 5 | 0 | -5 |
| Active Markdown documents | 16 | 19 | +3 |
| Historical Markdown documents | 103 | 104 | +1 |
| Specialized research documents | 2 | 2 | 0 |

The research/manuscript and post-docking archive moves are organization-only
and preserve blob content. Actual footprint reduction comes from removing the
generated audit artifacts, dangling private symlink, and two broken obsolete
shell scripts, partly offset by the checker/tests, evidence indexes, and
extracted orchestration modules.

Largest affected top-level functions changed as follows:

- parser construction: 465-line `workflow.cli.build_parser` to a 23-line
  composer in `workflow.cli_parser`, with a 223-line largest declaration group;
- `generate_condor_deployment`: 217 to 135 lines;
- `generate_slurm_deployment`: 208 to 127 lines.

## Review and external status

Independent Standards and Spec reviewers examined the fixed-base diff. Their
privacy and evidence-record findings produced the final documentation-only
fixes; scoped re-review is recorded in the task ledger. Draft pull-request CI
and its Linux/Windows Python 3.10/3.12 matrix are external pending evidence until
the branch is pushed.
