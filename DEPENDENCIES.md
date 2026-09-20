# Dependencies

## Core Python environment

`requirements.txt` supplies the default runtime libraries: NumPy, pandas,
Biopython, Matplotlib, seaborn, openpyxl, questionary, pdb-tools, SciPy,
PyYAML, requests, and scikit-learn. Package metadata requires Python 3.10+.
Builds use setuptools and wheel through `pyproject.toml`.

## Optional Python extras

Declared extras in `setup.py` are:

- `chemistry`: RDKit, Meeko, and gemmi;
- `interactions`: RDKit, PLIP, ProLIF, and MDAnalysis;
- `dev`: pytest, build, and wheel;
- `notebooks`: Jupyter.

Example:

```bash
python -m pip install -e '.[chemistry,interactions,dev]'
```

Optional imports fail or skip their stage when absent; installing an extra does
not install every external executable used by that stage.

## External preparation and docking tools

Choose tools that match the selected commands and engine profile:

- Open Babel and/or AutoDockTools (`prepare_ligand4.py`,
  `prepare_receptor4.py`) for backends that explicitly use them;
- GNINA for GNINA docking and CNN scores;
- AutoDock Vina for Vina docking;
- Smina for Smina docking;
- AutoDock4 and AutoGrid4, including the required parameter file, for AD4.

Engine binaries, containers, GPU drivers, and force-field/parameter assets are
not bundled. Record the actual executable identity and parameters; resume is
disabled when identity cannot be verified.

## Analysis and visualization tools

PLIP and ProLIF are optional interaction backends. LigPlot+, PoseView, PyMOL,
py3Dmol, and other visual tools are optional and stage-specific. Their outputs
do not replace source structures, score provenance, coverage accounting, or
manual inspection.

## Deployment tools

Local generation works through Python. Remote execution additionally needs a
Unix environment, SSH/file synchronization, the selected Slurm or HTCondor
worker commands, and any engine runtime described by the profile. Start from the
public-safe templates in `examples/hpc_profiles/`; never commit credentials or
private paths.

See [installation](INSTALLATION_GUIDE.md), [HPC deployment](HPC_DEPLOYMENT_GUIDE.md),
and [testing](docs/testing.md).
