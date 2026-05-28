# Feature Specification: Docking-First Project Workflow

## Context

The primary workflow is now a docking project lifecycle rather than a generic
collection of unrelated modes. Users should start by creating or resuming a
project, gather proteins and ligands into staged folders, build a pairlist, and
 then materialize docking-ready folders for GNINA, Vina, and Smina.

## Functional Requirements

- FR-001: The canonical external layout for interactive docking projects must
  be the staged `docking_legacy` layout with:
  - `1-Raw_Ligand/`
  - `2-Raw_Protien/`
  - `3-Preparation/`
  - `4-Docking/`
  - `5-Post_Docking_Analysis/`
  - `.workflow/state.json`
- FR-002: PDB collection must feed the staged raw directories directly.
- FR-002a: `1-Raw_Ligand/` must act as the single ligand staging directory and
  may contain both source ligands and normalized SDF companions.
- FR-003: Pairlist generation must be explicit and cocrystal-first.
- FR-004: Docking materialization must preserve GNINA HPC-compatible folders
  inside `4-Docking/` while also provisioning Vina and Smina output roots.
- FR-005: Post-docking analysis commands must accept the staged project root
  and resolve the correct docking and analysis subdirectories automatically.
- FR-006: Workflow initialization must not classify a fresh staged project as a
  legacy GNINA project based only on a nearby parent `pairlist.csv`; legacy
  registration requires real GNINA output signals inside the project.

## Acceptance Criteria

1. A fresh staged project contains the expected raw, preparation, docking, and
   post-docking directories.
2. GNINA-compatible folders exist inside `4-Docking/` when GNINA is enabled.
3. Vina and Smina output roots exist inside `4-Docking/` when those engines
   are enabled.
4. Post-docking helpers can resolve `4-Docking/gnina_out`, `4-Docking/logs`,
   `4-Docking/receptors`, and `4-Docking/pairlist.csv` from the staged project
   root.
5. Creating a new staged project inside a larger research directory that
   already contains a parent-level `pairlist.csv` still produces a staged
   project manifest instead of a false legacy GNINA registration.
