---
name: bio-project-intake
description: Intake and validate bioinformatics project requirements, experimental designs, biological targets, and computational scopes before executing workflows.
---

# Bio-Project Intake Skill

## Trigger Conditions
Activate when starting a new bioinformatics project, receiving new biological targets/compounds, or when project scope is initialized or modified.

## Explicit Inputs and Outputs
- **Inputs:** Project objective, target receptor accession/structure, ligand definitions/SMILES/SDF, experimental hypothesis, computational resources.
- **Outputs:** Intake summary document, validated file paths, parameter intake manifest.

## Refusal-to-Guess Rules
- **NEVER** guess receptor chain ID, binding pocket coordinates, or ligand protonation states without explicit specification or co-crystallized structural anchor.
- **NEVER** assume an unstated genome build or reference annotation version.
- **STOP** and ask the Scientific Lead whenever biological scope is ambiguous.

## Verification Commands
- `python -m workflow.cli prep --help`
- `python main.py workflow status --project-dir <path>`

## Known Limitations & Data-Use Notes
- Only operates on local, legally redistributable files or public accessions (e.g. RCSB PDB).
- Does not upload data to external servers.
