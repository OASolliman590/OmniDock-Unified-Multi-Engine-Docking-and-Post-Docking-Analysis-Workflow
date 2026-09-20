# Testing and evidence

## Local gates

From a clean checkout with Python 3.10+ and the development/chemistry extras:

```bash
python -m pip install -e '.[dev,chemistry]'
export PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
export PYTHONUTF8=1
export MPLBACKEND=Agg
python -m pytest test -q
python scripts/check_repository.py
python -m build
python scripts/check_wheel.py
git diff --check
```

PowerShell equivalent:

```powershell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD = '1'
$env:PYTHONUTF8 = '1'
$env:MPLBACKEND = 'Agg'
python -m pytest test -q
python scripts/check_repository.py
python -m build
python scripts/check_wheel.py
git diff --check
```

These variables isolate pytest from unrelated globally installed plugins, force
UTF-8 mode, and select a headless Matplotlib backend. They match the CI test
environment and should remain set for the complete local gate.

Run all help surfaces named in changed documentation. For shell changes, run
`bash -n` on each tracked shell script from Bash/WSL/Linux; PowerShell is not a
Bash parser. The repository checker enforces repository hygiene, including
documentation/link rules. The wheel checker verifies distribution contents and
entrypoint imports; it does not run docking engines.

## CI expectations

`.github/workflows/correctness.yml` defines a clean matrix for Ubuntu and
Windows on Python 3.10 and 3.12. Each job installs `.[dev,chemistry]`, runs the
test suite and repository checker, builds the distribution, and checks wheel
contents. A workflow definition is an expectation, not evidence that a remote
run occurred; link the actual run when claiming CI is green.

## Windows limitations

The core suite is intended to run on Windows, but optional chemistry binaries
may not. The dated 2026-09-12 record reports that Windows Application Control
blocked RDKit `rdDetermineBonds`, so the real Meeko macrocycle export integration
was skipped on that host. Use an allowed RDKit installation or a clean Linux
environment to verify that integration. Native Windows support also does not
imply native support for Bash preparation, CUDA, every external binary, or a
cluster scheduler.

## What counts as current evidence

Current evidence names the commit/tree, date, operating system, Python version,
exact command, exit status, and important skips/warnings. Engine, scheduler,
GPU, preparation-backend, or interaction claims require runs of those real
components. A mocked/unit test, dry run, generated manifest, or old log does not
prove them.

[Verification history](../audit/verification-history.md) preserves a dated
2026-09-12 local record. It is not a rerun for the current tree. Likewise,
numbered specs, release checklists, and old generated reports are historical
evidence only. Missing coverage and `not_evaluable` outcomes must be reported,
not silently counted as passes.
