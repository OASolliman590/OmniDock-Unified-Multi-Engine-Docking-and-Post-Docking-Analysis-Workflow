# Feature Specification: Ligand Normalization And Raw-Ligand Unification

## Context

The ligand-preparation path must be reliable for Vina-family engines when users
provide `.mol` files or other ligands that only contain 2D coordinates.
Meeko expects explicit hydrogens and 3D coordinates, so the workflow needs a
scientifically grounded normalization step before PDBQT generation. At the same
time, project initialization should stop creating two raw-ligand staging
folders for the same compounds.

## Functional Requirements

- FR-001: Ligand preparation must normalize `.sdf`, `.mol`, `.mol2`, and `.pdb`
  inputs to an explicit-hydrogen 3D SDF intermediate before generating PDBQT.
- FR-002: Existing 3D ligand coordinates should be preserved when present;
  2D ligands must be embedded to 3D before Meeko is invoked.
- FR-003: The shell preparation path, PDB-extracted ligand path, and ligand
  repair path must all use the same normalization contract.
- FR-004: New staged projects must use a single raw-ligand staging directory.
  SDF-normalized copies should be colocated with the original raw ligand files
  instead of being written into a separate sibling folder.
- FR-005: Manifest and workflow state compatibility fields may remain, but they
  must resolve to the same physical raw-ligand directory so older resume and
  repair code does not break.

## Acceptance Criteria

1. A 2D `.mol` ligand can be prepared into a valid PDBQT for Vina-family
   docking.
2. A flat `.sdf` ligand without 3D coordinates can be prepared into a valid
   PDBQT.
3. A PDB-extracted co-crystal ligand still prepares successfully through the
   shared pipeline.
4. `workflow init --layout-profile docking_legacy` no longer creates a separate
   `1-Raw_Ligand_SDF/` directory.
5. `pdb collect` stores source ligands and normalized SDF companions in the
   same raw-ligand folder.
