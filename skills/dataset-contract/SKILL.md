---
name: dataset-contract
description: Validate structure and ligand datasets against schema, integrity, checksum, and format contracts.
---

# Dataset Contract Skill

## Trigger Conditions
Activate when adding, parsing, or normalizing receptor PDB/PDBQT files, ligand SDF/PDBQT/MOL2 files, or pairlists.

## Explicit Inputs and Outputs
- **Inputs:** Raw receptor/ligand files, pairlist CSV, metadata schema.
- **Outputs:** Validation report, SHA-256 checksum manifest, sanitized file paths.

## Refusal-to-Guess Rules
- **NEVER** silently discard unparseable poses, missing atoms, or malformed charge records.
- **NEVER** assume bond orders or explicit hydrogens if absent from input.

## Verification Commands
- `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix`

## Known Limitations & Data-Use Notes
- Validates file formatting and structural integrity; does not perform molecular dynamics equilibration.
