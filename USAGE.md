# Usage

This guide provides minimal, practical command sets for OmniDock workflows.

## 1) Setup

### One-shot dependency download/install (Conda recommended)

```bash
conda env create -f environment.yml && conda activate pdb-prepare-wizard && pip install -e .
```

### One-shot dependency download/install (pip-only fallback)

```bash
python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt && pip install -e .
```

Then verify:

```bash
python main.py --help
```

---

Standard step-by-step:

```bash
conda env create -f environment.yml
conda activate pdb-prepare-wizard
pip install -e .
python main.py --help
```

## 2) Interactive Workflow

```bash
python main.py workflow interactive
```

Use this for guided project setup, preparation, docking orchestration, HPC flow, and analysis.

## 3) Non-Interactive Core Commands

### Initialize a project

```bash
python main.py workflow init --project-dir <project_dir> --layout-profile canonical
```

### Build pairlist

```bash
python main.py prep pairlist --project-dir <project_dir>
```

### Run docking

```bash
python main.py dock run --project-dir <project_dir> --engine gnina
```

### Analyze favorite engine

```bash
python main.py analyze favorite-engine \
  --project-dir <project_dir> \
  --favorite-engine gnina
```

### Run clean interaction pipeline

```bash
python main.py analyze clean \
  --project-dir <project_dir> \
  --dataset-root <dataset_root>
```

## 4) Interaction Contract

- Standalone `prolif/ligplot/pandamap` aliases are routed to `analyze clean`.
- Mandatory clean stages are enforced through the clean-contract path.
- Retired legacy targets:
  - `analyze.stage.structure_quality`
  - `analyze.visuals.pymol`

## 5) HPC Profile Templates

- Use template files from `examples/hpc_profiles/`.
- Keep real SSH configs in `.workflow/` (ignored by git).

Template example:

`examples/hpc_profiles/ssh_systems.template.yaml`

## 6) Useful Environment Variables

- `DOCKFORGE_CLEAN_INTERACTION_DATASET_ROOT`
- `DOCKFORGE_RESEARCH_BASE_PATH`
- `LIGPLUS_ROOT` / `LIGPLUS_HOME`

## 7) Smoke Check

```bash
MPLCONFIGDIR=/tmp/mpl PYTHONPATH=".:test" python test/test_dockforge_smoke.py
```
