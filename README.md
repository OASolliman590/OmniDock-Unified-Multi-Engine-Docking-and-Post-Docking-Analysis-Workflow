# OmniDock

OmniDock is a Python 3.10+ workflow for preparing receptor/ligand projects,
running GNINA, AutoDock Vina, Smina, or AutoDock4, and performing provenance-
aware post-docking analysis. The Python distribution and installed commands
retain the legacy `pdb-prepare-wizard` / `pdb-wizard-*` names for compatibility.

## Start here

```bash
python -m pip install -e '.[chemistry,dev]'
python main.py --help
python main.py workflow init --project-dir docking_project --layout-profile canonical
python main.py workflow status --project-dir docking_project
```

The explicit layout flag selects the canonical project topology. Omitting it
selects the `docking_legacy` compatibility default.

The canonical documentation index is [docs/README.md](docs/README.md). See the
[installation guide](INSTALLATION_GUIDE.md), [CLI recipes](USAGE.md), and
[architecture guide](docs/architecture.md) for the supported path.

## Workflow

1. Initialize or stage a canonical project with `--layout-profile canonical`.
2. Prepare receptors and ligands and review the preparation/QC records.
3. Materialize `pairlist.csv` and the project manifest.
4. Dry-run, then execute docking locally or create a scheduler deployment.
5. Analyze one engine or compare compatible, direction-aware normalized scores.
6. Inspect provenance, coverage, `not_evaluable` rows, and generated reports.

Use `python main.py workflow interactive` for a guided shell. Installed legacy
entrypoints remain available, including `pdb-wizard-workflow` (declared as
`workflow.cli:main`), but examples use `python main.py` so the repository
checkout and installed package expose the same command tree.

## External tools

Requirements depend on the selected path:

- preparation: RDKit/Meeko and, for selected backends, Open Babel or
  AutoDockTools;
- docking: the chosen GNINA, Vina, Smina, AutoDock4/AutoGrid4 executable or
  container;
- interactions and visuals: optional PLIP, ProLIF, LigPlot+, PoseView, PyMOL,
  or py3Dmol integrations;
- remote execution: a Unix worker plus the selected Slurm or HTCondor tools,
  SSH, and file synchronization utilities.

Installing a Python extra does not install external binaries, GPU drivers, or a
scheduler. See [DEPENDENCIES.md](DEPENDENCIES.md).

## Scientific and operational limits

Docking scores are model outputs, not measured affinities. Energy-like scores
are lower-is-better; GNINA CNNscore/CNNaffinity and OmniDock consensus scores
are higher-is-better. Unlike scoring functions are normalized within their own
engine/function and target/site before comparison.

Redocking RMSD is measured in the native receptor coordinate frame and requires
an authoritative chemical graph and atom mapping. Missing mappings, mismatched
graphs, missing poses, or missing frame provenance are `not_evaluable`, not
passes. Missing comparisons reduce coverage. Preparation still requires human
decisions about assemblies, protonation, metals/cofactors, missing residues,
boxes, and unsupported chemistry.

Regression tests enforce software invariants; they do not establish enrichment,
affinity accuracy, clinical usefulness, or broad scientific validity. Real
engines, clusters, CUDA, and optional scientific backends must be validated in
the intended environment. See [docs/correctness-update.md](docs/correctness-update.md)
for implemented behavior and explicit limitations.

## Project status and history

This repository includes historical specifications, audits, release snapshots,
and specialized research artifacts. They are evidence and design history, not
current usage instructions or proof of a green build. The documentation index
labels each collection. Contributions should follow [CONTRIBUTING.md](CONTRIBUTING.md)
and [docs/testing.md](docs/testing.md).
