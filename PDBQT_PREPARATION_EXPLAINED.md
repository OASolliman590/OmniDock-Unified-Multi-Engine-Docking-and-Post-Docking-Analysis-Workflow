# DockForge PDBQT Preparation Explained

## Quick Summary

DockForge has three preparation entry points:
1. `autodock_preparation.py` (Python orchestrator)
2. `prep_autodock_enhanced.sh` (shell runner for orchestrated prep)
3. `docking/preparation/ligand_preparation.py` (ligand-focused preparation module)

All produce docking-ready PDBQT assets.
Preparation is decoupled from PLIP and other interaction-analysis tools.

## Why PDBQT

PDBQT extends PDB by adding:
- atom typing for AutoDock-family workflows
- partial charges used by docking engines

## Current Ligand Preparation Design

Ligand preparation is profile-driven and engine-aware.

### Core enrichment (Open Babel path)

When Open Babel normalization is used, DockForge applies:
- hydrogen addition
- protonation at chosen pH (default `7.4`)
- 3D conformer generation
- minimization
- partial charge assignment (best effort)

### Supported profiles

- `engine_aware_full`
- `openbabel_only`
- `meeko_only`
- `autodocktools_only`
- `openbabel_meeko`
- `openbabel_meeko_autodock`
- `openbabel_autodocktools`

### Engine-aware behavior

`engine_aware_full` resolves automatically:
- with `autodock4` selected: `openbabel_meeko_autodock`
- without `autodock4`: `openbabel_meeko`

Invalid profile/engine combinations are blocked before preparation starts.

## Receptor Preparation Design

Receptor preparation builds PDBQT outputs for docking.
The pipeline prefers stronger tools when available and falls back safely if optional tools are missing.

## CLI Examples

```bash
# Create default config snapshot
python autodock_preparation.py --create-config

# Engine-aware prep for multi-engine run
python autodock_preparation.py \
  --ligands-input ./ligands_raw \
  --receptors-input ./receptors_raw \
  --ligands-output ./ligands_prep \
  --receptors-output ./receptors_prep \
  --ligand-profile engine_aware_full \
  --selected-engines gnina,vina,smina,autodock4

# Open Babel -> Meeko only
python autodock_preparation.py \
  --ligands-input ./ligands_raw \
  --receptors-input ./receptors_raw \
  --ligands-output ./ligands_prep \
  --receptors-output ./receptors_prep \
  --ligand-profile openbabel_meeko \
  --selected-engines gnina,vina,smina
```

## AutoDockTools Resolution

When profile logic requires AutoDockTools, DockForge resolves `prepare_ligand4.py` from:
- explicit CLI path
- environment variables (`AUTODOCKTOOLS_PREPARE_LIGAND4`, `ADT_PREPARE_LIGAND4`)
- PATH lookup
- known bundled/local fallback locations

Optional companion settings:
- `--autodocktools-prepare-receptor4`
- `--autodocktools-python`

## Output Artifacts

Typical outputs:
- prepared ligand PDBQT files
- prepared receptor PDBQT files
- preparation summary report
- ligand validation report for malformed/invalid PDBQT outputs

## Troubleshooting

### Dependency missing

If preparation reports missing dependencies, verify:
- `obabel`
- `jq`
- AutoDockTools script path (if profile needs it)

### Profile/engine incompatibility

If you see compatibility blocking, switch to:
- `engine_aware_full`
- `openbabel_meeko_autodock`
- `openbabel_autodocktools`
- `autodocktools_only`
for workflows that include `autodock4`.
