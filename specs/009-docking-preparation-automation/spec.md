# Feature Specification: Docking Preparation Automation

## Context

Docking preparation is now part of a staged docking project, not just a flat
canonical asset bundle. The preparation flow must populate raw protein, raw
ligand, prepared protein, prepared ligand, and docking-ready directories in a
way that stays aligned with GNINA HPC-style projects while still serving Vina
and Smina.

## Goals

1. Build staged docking projects with raw and prepared asset directories.
2. Support project-aware PDB collection that writes source PDBs, cleaned
   receptors, extracted co-crystal ligands, and workbook updates into the
   project tree.
3. Add a dedicated pairlist-building stage before docking project
   materialization.
4. Materialize GNINA, Vina, and Smina-ready docking folders from prepared
   assets and an explicit pairlist.

## Functional Requirements

- FR-001: `python main.py workflow init --layout-profile docking_legacy` must
  scaffold the staged docking layout.
- FR-002: `python main.py pdb collect --project-dir <dir> --pdbs ...` must:
  - save source PDB and cleaned receptor files into `2-Raw_Protien/`
  - save extracted co-crystal ligands into `1-Raw_Ligand/`
  - write SDF-normalized ligand copies beside the source ligand inside `1-Raw_Ligand/` when possible
  - append/update `2-Raw_Protien/multi_pdb_analysis.xlsx`
- FR-003: `python main.py prep pairlist` must build `pair_intent.csv` and
  `4-Docking/pairlist.csv` from prepared assets plus workbook site metadata.
- FR-004: Pairlist generation must support:
  - `cocrystal_only`
  - `cocrystal_plus_all`
  - `curated_cartesian`
  - `curated_per_protein`
- FR-005: Pairlist generation must preserve cocrystal benchmark rows when
  curated or expanded modes are used.
- FR-006: `python main.py prep project` must materialize:
  - `4-Docking/receptors/`
  - `4-Docking/ligands/`
  - `4-Docking/metadata/`
  - selected engine output folders
  from prepared assets and an explicit pairlist. Legacy
  `python main.py prepare-docking ...` must continue to normalize to this
  grouped surface.
- FR-007: Default docking values remain:
  - `site_id = site_1`
  - `size_x = size_y = size_z = 20`

## Acceptance Criteria

1. A fresh `workflow init --layout-profile docking_legacy` creates the staged
   docking directory structure.
2. `python main.py pdb collect --help` exists, and the live fetch/workbook
   update path remains a rerun-backed verification item.
3. `python main.py prep pairlist --help` exists and supports all four pair
   generation modes.
4. `python main.py prep project --project-dir <dir>` can materialize a staged
   docking project using the project defaults and an existing pairlist, while
   legacy `python main.py prepare-docking ...` remains an alias for the grouped
   project-build surface.
5. `python -m py_compile` passes for the preparation and workflow modules.
