# Installation

OmniDock requires Python 3.10 or newer. Package metadata still uses the legacy
distribution name `pdb-prepare-wizard`; this is expected.

## Development checkout

Create and activate an isolated environment, then install the editable package:

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e '.[chemistry,dev]'
```

On shells that interpret single quotes differently, quote the extras expression
appropriately. The `environment.yml` file is an alternative Conda environment
definition and currently selects Python 3.12.

Core-only installation is:

```bash
python -m pip install -e .
```

Install optional interaction libraries with `.[interactions]`. Extras install
Python packages only; external programs and schedulers remain separate.

## Verify the installation

```bash
python main.py --help
python main.py workflow --help
python main.py pdb --help
python main.py prep --help
python main.py dock --help
python main.py analyze --help
python -m post_docking_analysis --help
python scripts/check_repository.py
```

For a development validation run, follow [docs/testing.md](docs/testing.md).
Installed commands include `pdb-prepare-wizard`, `pdb-wizard-workflow`,
`pdb-wizard-prepare-docking`, `pdb-wizard-dock`, and
`post-docking-analysis`. The repository examples prefer `python main.py`.

## Add workflow-specific tools

Install and validate only the components you intend to use:

- RDKit/Meeko for the chemistry path;
- Open Babel or AutoDockTools for preparation backends that name them;
- GNINA, Vina, Smina, and/or AutoDock4 plus AutoGrid4 for docking;
- PLIP/ProLIF and optional LigPlot+, PoseView, PyMOL, or py3Dmol for selected
  analysis/visualization stages;
- SSH/synchronization utilities and Slurm or HTCondor clients/workers for remote
  deployment.

Confirm executable versions and licensing from their upstream projects. See
[DEPENDENCIES.md](DEPENDENCIES.md) and the
[HPC guide](HPC_DEPLOYMENT_GUIDE.md).

## Platform notes

The CI definition targets Windows and Ubuntu on Python 3.10/3.12, but a CI
matrix does not guarantee every optional backend on each host. Bash preparation
and scheduler scripts require a Unix shell (native Linux, WSL, or a remote
worker). CUDA paths require a compatible driver/runtime. A dated Windows run
recorded an RDKit DLL blocked by Application Control; see
[testing notes](docs/testing.md#windows-limitations).

## First project

```bash
python main.py workflow init --project-dir docking_project --layout-profile canonical
python main.py workflow status --project-dir docking_project
```

Omitting `--layout-profile canonical` selects the `docking_legacy`
compatibility default.

Continue with the [CLI recipes](USAGE.md). Do not treat successful installation
as validation of a preparation protocol, docking engine, or scientific result.
