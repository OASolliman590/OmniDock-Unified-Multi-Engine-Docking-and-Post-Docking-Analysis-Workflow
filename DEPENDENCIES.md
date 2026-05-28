# Dependencies

This file documents required and optional dependencies for OmniDock.

## Python

- Python: `3.9+` (tested on `3.9`)

## Required Python Packages

From `requirements.txt`:

- `numpy>=1.21.0`
- `pandas>=1.3.0`
- `biopython>=1.79`
- `matplotlib>=3.5.0`
- `seaborn>=0.11.0`
- `openpyxl>=3.0.0`
- `questionary>=2.0.1`
- `pdb-tools>=2.5.0`

## Analysis Packages

- `plip>=2.2.0` (interaction analysis)
- `scipy` (statistical utilities)
- `scikit-learn` (clustering/model helpers)

## External Tools

- `openbabel` / `obabel` (structure conversions and chemistry prep)
- `LigPlot+` (optional)
- `PyMOL` (optional)
- `Apptainer`/`Singularity` (HPC container runs, optional)

## Optional Python Packages

- `prolif` (optional interaction maps)
- `py3Dmol` (optional 3D visualization)
- `jupyter>=1.0.0` (optional notebook workflows)

## Environment Setup

Quick one-liner (Conda):

```bash
conda env create -f environment.yml && conda activate pdb-prepare-wizard && pip install -e .
```

Quick one-liner (pip fallback):

```bash
python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt && pip install -e .
```

Recommended:

```bash
conda env create -f environment.yml
conda activate pdb-prepare-wizard
pip install -e .
```

Alternative:

```bash
pip install -r requirements.txt
pip install -e .
```

## Runtime Environment Variables

- `DOCKFORGE_CLEAN_INTERACTION_DATASET_ROOT` (optional dataset root for clean interaction pipeline)
- `DOCKFORGE_RESEARCH_BASE_PATH` (optional base path for research adapter)
- `LIGPLUS_ROOT` or `LIGPLUS_HOME` (optional LigPlot root)

## Security Note

- Do not commit local runtime config under `.workflow/` (SSH targets, local paths, state files).
- Commit templates only (for example `examples/hpc_profiles/ssh_systems.template.yaml`).
