---
name: omnidock-project-intake
description: Intake and configure multi-engine docking projects in OmniDock/DockForge, verifying pairlists, engine selections, and output directories.
---

# OmniDock Project Intake Skill

## Trigger Conditions
Activate when initializing an OmniDock project (`workflow init`, `prep project`), setting up pairlists, or selecting docking engines (GNINA, Vina, Smina, AD4).

## Explicit Inputs and Outputs
- **Inputs:** Project root directory, receptor PDB files, ligand files, engine selections, box coordinates.
- **Outputs:** Materialized project layout (`0-Input/` to `7-Reports/`), `pairlist.csv`, `.workflow/state.json`.

## Refusal-to-Guess Rules
- **NEVER** guess docking box dimensions; require pocket coordinates or co-crystal reference.
- **NEVER** overwrite an existing project without explicit confirmation.

## Verification Commands
- `python main.py prep project --help`
- `python main.py workflow status --project-dir <path>`
