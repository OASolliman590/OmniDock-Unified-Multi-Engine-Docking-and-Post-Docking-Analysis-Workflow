---
name: reference-build-audit
description: Audit structural and genomic reference versions, coordinate systems, and annotation sources.
---

# Reference Build Audit Skill

## Trigger Conditions
Activate when importing reference structures, PDB entries, UniProt sequences, or coordinate grids.

## Explicit Inputs and Outputs
- **Inputs:** PDB accession code, resolution, experimental method, deposition date, coordinate frames.
- **Outputs:** Reference audit record, resolution/clash analysis, binding site coordinates.

## Refusal-to-Guess Rules
- **NEVER** assume an apo structure has the same pocket geometry as a holo structure.
- **NEVER** silently swap reference structures during redocking or RMSD validation.

## Verification Commands
- Inspect PDB REMARK headers and crystallographic resolution.

## Known Limitations & Data-Use Notes
- Audits metadata and coordinates; does not re-refine electron density maps.
