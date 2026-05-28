# Feature Specification: Unified Output Topology

**Feature Branch**: `021-unified-output-topology`
**Created**: 2026-04-01
**Status**: Draft
**Input**: User observation that post-docking results are fragmented across `5-Post_Docking_Analysis/`, `5-Analysis/`, `4-PostDocking/`, and `7-Reports/` with no clear entry point.

---

## Summary

Post-docking analysis currently writes results to four separate sibling directories under the project root. A user running a single analysis session must navigate four different folders to find complexes, hit tables, interaction maps, and reports. There is no single folder they can open to review their run.

This spec collapses the fragmented output layout into one canonical analysis root (`5-Analysis/`) with a clear internal hierarchy and a machine-generated `START_HERE.md` at the top level. Intermediate working data stays in `4-PostDocking/` (renamed `4-Working/`). Consolidated cross-run reports remain in `7-Reports/`.

---

## Problem

### Observed Folder Split (current)

```
4-PostDocking/          ← intermediate: filtering, rescoring, raw scores, SQLite
5-Post_Docking_Analysis/ ← session data, complexes, scores, RMSD, interactions,
                           visualizations, reports  ← MAIN OUTPUT (old)
5-Analysis/             ← comparative tables, top_pose, polypharmacology,
                           interaction maps         ← MAIN OUTPUT (new)
7-Reports/              ← consolidated cross-run reports, dashboard export
```

### Root Causes

- **RC-A**: Two pipeline generations (`simplified_pipeline.py` → `5-Post_Docking_Analysis/`, `multi_engine_pipeline.py` → `5-Analysis/`) each chose their own output root. Spec 020 unified the pipelines but did not reconcile the output roots.
- **RC-B**: Spec 016 introduced `5-Analysis/` as the canonical root for comparative and top-pose artifacts but did not migrate per-session data away from `5-Post_Docking_Analysis/`.
- **RC-C**: No top-level navigation document tells users which folder to open first.
- **RC-D**: `4-PostDocking/` name implies it is a peer of `5-Analysis/` for post-docking outputs, not intermediate scratch space.

---

## Goals

1. One canonical output root for all human-facing post-docking results.
2. Clear separation between final outputs (`5-Analysis/`) and working/intermediate data (`4-Working/`).
3. A generated `START_HERE.md` at the project root listing what was produced and where.
4. No duplicate or orphaned output paths after a run.
5. Backward-compatible symlinks for `5-Post_Docking_Analysis/` during a transition window.

## Non-Goals

- Changing what is computed (this is output routing only).
- Restructuring the docking output folders `0-Input` through `3-Docking`.
- Merging `7-Reports/` into `5-Analysis/` (cross-run consolidation stays separate).
- Changing session subdirectory naming inside `5-Analysis/`.

---

## User Stories

### US1: Single Folder to Open
As a scientist, I want to open one folder after a run and find all my results so I do not need to search across four directories.

### US2: Clear Intermediates vs. Outputs
As an operator, I want intermediate pipeline files (filtering manifests, rescoring shortlists, SQLite cache) clearly separated from the final analysis artifacts so I can archive or delete them independently.

### US3: Navigation Document
As a new user, I want a generated file at the root of the analysis folder that tells me what was produced, where each artifact type lives, and which session is the most recent.

### US4: Safe Transition
As an existing user with prior runs already written to `5-Post_Docking_Analysis/`, I want those results to remain accessible (via symlink or explicit path note) so nothing is silently lost.

---

## Functional Requirements

### Output Root Unification

- **FR-001**: All human-facing post-docking outputs MUST be written under `5-Analysis/` as the single canonical root.
- **FR-002**: `5-Post_Docking_Analysis/` MUST NOT be created by any new pipeline run after this spec is implemented.
- **FR-003**: During a configurable transition window, the pipeline MUST create a compatibility symlink `5-Post_Docking_Analysis → 5-Analysis` if the old directory does not already exist, and emit a deprecation notice in the run log.
- **FR-004**: If `5-Post_Docking_Analysis/` already exists from a prior run, the pipeline MUST NOT overwrite or delete it; it MUST log a notice directing the user to `5-Analysis/` for new results.

### Canonical Folder Hierarchy

- **FR-005**: `5-Analysis/` MUST use the following internal structure:

```
5-Analysis/
├── START_HERE.md                        ← generated navigation guide
├── LATEST_SESSION -> sessions/<id>/     ← symlink to most recent session
├── sessions/
│   └── <session_id>/                    ← one folder per analysis run
│       ├── scores/                      ← raw and normalized score CSVs
│       ├── best_poses/                  ← best-pose CSV + PDBQT files
│       ├── complexes/                   ← merged receptor+ligand PDB files (MD-ready)
│       ├── rmsd_analysis/               ← per-complex RMSD matrices and clusters
│       ├── interactions/
│       │   ├── prolif/                  ← frequency tables, barcode, network
│       │   ├── pandamap/                ← 2D interaction maps
│       │   ├── poseview/                ← PoseView diagrams
│       │   └── pymol/                   ← PyMOL scripts and sessions
│       ├── visualizations/
│       │   ├── 2d/                      ← all 2D visual outputs (mirrored)
│       │   └── 3d/                      ← all 3D visual outputs (mirrored)
│       └── reports/                     ← per-session HTML/PDF/MD reports
├── comparative/                         ← cross-session hit ranking tables
├── top_pose_ligand_performance/         ← top-pose atlas artifacts (spec 019)
├── polypharmacology/                    ← multi-target scoring tables
├── structure_quality/                   ← validation gate outputs
└── raw_data/                            ← manifests, indexes, inventory CSVs
    ├── START_HERE.md                    ← (duplicate of root for quick copy)
    ├── visualization_manifest.csv
    ├── raw_artifact_inventory.csv
    └── outputs_index.csv                ← run tracking index
```

- **FR-006**: The `complexes/` subdirectory under each session MUST be the canonical source for MD-ready complex PDB files.
- **FR-007**: `comparative/`, `top_pose_ligand_performance/`, and `polypharmacology/` MUST be written at the `5-Analysis/` root level (not inside a session folder) because they aggregate across sessions.

### Working Data Separation

- **FR-008**: `4-PostDocking/` MUST be renamed `4-Working/` to clearly signal intermediate/scratch status.
- **FR-009**: `4-Working/` MUST contain: `filtering/`, `rescoring/`, `rerun_manifests/`, `scores/` (raw pre-filtering), `storage/` (SQLite cache).
- **FR-010**: Nothing in `4-Working/` is required for scientific review; it MAY be deleted after a run without loss of reportable results.

### Navigation Document

- **FR-011**: The pipeline MUST generate `5-Analysis/START_HERE.md` at the end of every run (consolidation step).
- **FR-012**: `START_HERE.md` MUST include:
  - Run timestamp and session ID
  - Table of key artifacts with relative paths and one-line descriptions
  - Explicit path to `complexes/` for MD users
  - Explicit path to `top_pose_ligand_performance/` for hit-ranking users
  - Explicit path to `7-Reports/` for consolidated cross-run users
  - Warning if `5-Post_Docking_Analysis/` still exists with prior data

### Cross-Run Reports

- **FR-013**: `7-Reports/` remains the canonical location for consolidated cross-run artifacts (`consolidated_run_summary.csv`, dashboard export contract, per-protein figures/tables).
- **FR-014**: `START_HERE.md` MUST link to `7-Reports/` for cross-run comparisons.

---

## Migration Contract

| Old Path | New Path | Action |
|---|---|---|
| `5-Post_Docking_Analysis/analysis/sessions/` | `5-Analysis/sessions/` | Write here directly |
| `5-Post_Docking_Analysis/complexes/` | `5-Analysis/sessions/<id>/complexes/` | Write here directly |
| `5-Post_Docking_Analysis/best_poses/` | `5-Analysis/sessions/<id>/best_poses/` | Write here directly |
| `5-Post_Docking_Analysis/scores/` | `5-Analysis/sessions/<id>/scores/` | Write here directly |
| `5-Post_Docking_Analysis/rmsd_analysis/` | `5-Analysis/sessions/<id>/rmsd_analysis/` | Write here directly |
| `5-Post_Docking_Analysis/interactions/` | `5-Analysis/sessions/<id>/interactions/` | Write here directly |
| `5-Post_Docking_Analysis/visualizations/` | `5-Analysis/sessions/<id>/visualizations/` | Write here directly |
| `5-Post_Docking_Analysis/reports/` | `5-Analysis/sessions/<id>/reports/` | Write here directly |
| `5-Analysis/comparative/` | `5-Analysis/comparative/` | No change |
| `5-Analysis/top_pose_ligand_performance/` | `5-Analysis/top_pose_ligand_performance/` | No change |
| `5-Analysis/interactions/` | `5-Analysis/sessions/<id>/interactions/` | Moved into session |
| `4-PostDocking/` | `4-Working/` | Rename |
| `7-Reports/` | `7-Reports/` | No change |

---

## Acceptance Criteria

1. After a new run, zero files are written to `5-Post_Docking_Analysis/`.
2. `5-Analysis/START_HERE.md` exists and contains paths to: complexes, top-pose atlas, hit tables, and `7-Reports/`.
3. `5-Analysis/LATEST_SESSION` symlink resolves to the most recently completed session folder.
4. `4-Working/` is created instead of `4-PostDocking/` on new runs.
5. A compatibility symlink `5-Post_Docking_Analysis → 5-Analysis` is present when the old path did not exist before the run.
6. Existing `5-Post_Docking_Analysis/` from prior runs is untouched after a new run.
7. `complexes/` under the session folder contains MD-ready PDB files accessible without navigating into any sub-subfolder.
8. Running `grep -r "5-Post_Docking_Analysis" post_docking_analysis/` returns only the compatibility shim and deprecation notice — no active write paths.

---

## Risks

- Existing downstream scripts or user notebooks that hardcode `5-Post_Docking_Analysis/` paths will break. Compatibility symlink mitigates this during transition.
- Session-level vs. root-level distinction for `comparative/` and `polypharmacology/` requires pipeline to know whether it is writing a per-session artifact or an aggregate artifact.

---

## Implementation Phases

### Phase 1: Output Root Switch
Redirect all write paths from `5-Post_Docking_Analysis/` to `5-Analysis/sessions/<session_id>/`. Add compatibility symlink logic (FR-002, FR-003, FR-004).

### Phase 2: Hierarchy Enforcement
Enforce canonical subfolder names inside each session. Move `interactions/`, `visualizations/` writes into session scope (FR-005, FR-006, FR-007).

### Phase 3: Working Data Rename
Rename `4-PostDocking/` to `4-Working/` across all write paths (FR-008, FR-009, FR-010).

### Phase 4: Navigation Document
Generate `START_HERE.md` in consolidation step (FR-011, FR-012).

### Phase 5: Verification
Smoke tests for all acceptance criteria. Grep assertion that no active write paths reference `5-Post_Docking_Analysis/`.
