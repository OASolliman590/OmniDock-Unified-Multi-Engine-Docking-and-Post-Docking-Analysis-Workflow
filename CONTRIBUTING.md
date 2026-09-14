# Contributing

Python **3.10+** is the supported package contract (`python_requires`, CI).
Older README/installation text that mentions Python 3.8 or 3.9 is historical.

Optional chemistry extras, interaction extras, docking engines, Open Babel,
cluster schedulers, CUDA, and ProLIF are **not** assumed to work on every
developer machine. Do not claim a local checkout exercises all of them.

## Development install

```bash
python -m pip install -e '.[dev,chemistry]'
```

Use `.[dev]` when chemistry wheels cannot load on the host. Installing an extra
does not install external binaries.

## Tests

Always disable plugin autoload.

Focused:

```bash
# POSIX
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test/test_repository_contracts.py -q

# PowerShell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD='1'
python -m pytest test/test_repository_contracts.py -q
```

Full suite (requires an environment where RDKit can load; some hosts block it):

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test -q
```

## Repository checker

```bash
python scripts/check_repository.py
```

Exit code `0` means no findings. Do not weaken the checker to hide drift.

## Build and wheel

```bash
python -m build
python scripts/check_wheel.py
```

## Shell syntax

On Linux or Git Bash, syntax-check every tracked `*.sh` file:

```bash
bash -n prep_autodock_enhanced.sh
bash -n install_vina_dock_dependencies.sh
```

## Public interfaces

Preserve public CLI flags, defaults, console-script names, and deployment
function signatures unless a written specification names an additive,
compatibility-safe change. The installed distribution name remains
`pdb-prepare-wizard` until a later, explicit naming migration.

## Scientific invariants

Do not regress authoritative chemical graphs and atom mappings; native
receptor-frame RMSD (no ligand superposition that hides displacement); explicit
score direction; explicit `not_evaluable` status and coverage accounting;
reference provenance; fail-closed execution and QC gates; or content-verified
resume/cache identity. Do not invent benchmark, affinity, efficacy, or
production-readiness claims.

## Historical documents

Numbered trees under `specs/` (`001`–`029` and `_archive/`), curated `audit/`
reports, and dated files under `docs/` are historical snapshots. Preserve their
claims; do not rewrite them as current validation. Current cleanup work is in
`docs/superpowers/`.

## Secrets and private HPC details

Do not commit credentials, private host/account/path facts, or filled cluster
profiles. Public-safe templates live in `examples/hpc_profiles/`. Keep real
site settings in ignored local configuration such as `.workflow/hpc_profiles/`.
