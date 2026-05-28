# Feature Specification: Engine-Aware Post-Docking Analysis

## Context

The current analysis stack is strongest for GNINA and largely assumes one
engine-specific score schema. The repository now needs a comparison layer that
understands multiple docking engines and can either compare all engines, focus
on one engine, or continue downstream analysis from a chosen favorite engine.

## Goals

1. Normalize GNINA, Vina, and Smina outputs into one analysis schema.
2. Add an engine-aware analysis mode to the post-docking CLI.
3. Support three analysis behaviors:
   - `single_engine`
   - `comparative_all_engines`
   - `favorite_engine_continue`
4. Emit unified score tables and best-pose tables that fit the same downstream
   score-based pipeline shape regardless of engine.
5. Reuse the existing simplified GNINA pipeline when the favorite engine is
   GNINA and canonical SDF/log folders are present.
6. Provide a structural downstream bridge for Vina and Smina favorites that
   can extract best-pose complex PDBs and continue with the safe subset of the
   simplified downstream stack, including affinity analysis, polypharmacology,
   generic visualizations, and artifact consolidation.

## Out of Scope

- Full parity with every GNINA-only downstream stage for non-GNINA engines in
  this phase, especially GNINA-specific all-pose SDF-native workflows and
  engine-specific interaction stages that depend on optional external tools.
- Replacing the existing simplified pipeline for old GNINA-only projects.

## Functional Requirements

- FR-001: `python -m post_docking_analysis --project-dir <project>` must detect
  engine-aware analysis mode.
- FR-002: Normalized scores must include engine, tag, protein, ligand, site,
  pose, affinity, and score metadata columns.
- FR-002a: The engine-aware pipeline must also emit compatibility exports with
  unified score columns such as `complex_name`, `mode`, and `vina_affinity`
  aliases so existing score-based tooling can consume them.
- FR-003: Comparative mode must emit combined score tables and per-engine best
  pose summaries.
- FR-004: Single-engine mode must emit filtered reports for the selected engine.
- FR-005: Favorite-engine mode must emit filtered reports and attempt GNINA deep
  downstream continuation when GNINA is selected.
- FR-006: Favorite-engine mode must also provide a non-GNINA structural bridge
  for Vina/Smina that extracts best-pose complex PDBs and writes continued
  downstream artifacts under `favorite_engine/deep_analysis/`.
- FR-007: The non-GNINA bridge should reuse the safe subset of the simplified
  downstream pipeline where possible.
- FR-007a: The non-GNINA bridge must run PDB-based RMSD analysis from extracted
  complex PDBs, including global best-pose and per-protein best-pose scopes
  whenever comparable structures exist.
- FR-007b: The non-GNINA bridge should also emit per-complex all-pose RMSD
  outputs when multiple comparable poses exist for the same tag.
- FR-008: Bridge stage summaries must distinguish explicit outcomes such as
  `completed`, `disabled`, `missing_dependency`, `missing_configuration`, and
  partial external-render cases instead of collapsing them into generic success.

## Acceptance Criteria

1. The post-docking CLI accepts `--project-dir`, `--analysis-mode`,
   `--engine`, and `--favorite-engine`.
2. Multi-engine analysis writes combined reports under `analysis/reports/`.
3. Multi-engine analysis writes unified score tables under `analysis/raw_data/`.
3. Favorite-engine mode writes filtered reports under
   `analysis/favorite_engine/`.
4. GNINA favorite-engine continuation can invoke the simplified pipeline against
   the canonical project layout.
5. Vina/Smina favorite-engine continuation writes best-pose complex PDBs,
   simplified affinity/polypharmacology outputs, continued visual/report
   artifacts, unified raw score tables, and PDB-based RMSD outputs under
   `analysis/favorite_engine/deep_analysis/`.
6. Bridge summary notes distinguish disabled stages from missing dependencies
   and missing tool configuration.
7. `python -m py_compile` passes for the new engine-aware analysis module.
