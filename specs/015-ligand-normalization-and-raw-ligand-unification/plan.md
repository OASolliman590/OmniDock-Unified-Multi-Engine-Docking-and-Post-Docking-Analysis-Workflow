# Implementation Plan: Ligand Normalization And Raw-Ligand Unification

## Scientific Anchor

- Meeko documents that ligand inputs must already contain explicit hydrogens
  and 3D coordinates before `mk_prepare_ligand.py` runs.
- The normalization path therefore creates a 3D SDF intermediate first, then
  hands that intermediate to Meeko or the Open Babel fallback.

## Implementation Notes

1. Add a shared Python ligand-normalization helper under
   `docking/preparation/` that:
   - preserves existing 3D coordinates when present,
   - embeds 3D coordinates for 2D ligands,
   - adds explicit hydrogens,
   - writes a normalized SDF before PDBQT generation.
2. Route the active shell preparation flow, the extracted-PDB ligand path, and
   the repair flow through the same helper.
3. Collapse `raw_ligands` and `raw_ligands_sdf` onto one physical project
   directory while keeping compatibility metadata keys mapped to that same
   location.
4. Refresh the user-facing docs/specs so the single-folder raw-ligand contract
   is visible in workflow and preparation guidance.
