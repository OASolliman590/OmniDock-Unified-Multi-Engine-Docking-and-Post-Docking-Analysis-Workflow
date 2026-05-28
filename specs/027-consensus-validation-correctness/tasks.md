# Tasks: Consensus and Validation Correctness

**Input**: `/specs/027-consensus-validation-correctness/spec.md`
**Scope**: 3 bug fixes, no new artifacts, no new public API surface beyond small optional parameters.

## Format: `[ID] [P?] Description`
- **[P]**: parallelizable (no dependency on other open tasks)

## Fix 1 — Geometric Consensus Cross-Format Soft Alignment

- [x] T001 Add `soft_alignment_used` and `soft_alignment_rmsd` fields to the `pair_row` dict schema in `_pairwise_rmsd()`. Default both to `False` / `None`.
- [x] T002 Add `pose_formats: Optional[Dict[str, str]]` parameter to `_pairwise_rmsd()`. When `atom_count_mismatch` occurs and formats of engine_a/engine_b differ AND `|count_a - count_b| <= 2`: centroid-sort both coordinate arrays, truncate to `min(count_a, count_b)`, compute Kabsch RMSD, record as `soft_alignment_used=True`. Increment `valid_pairs` and `pass_pairs` based on soft RMSD.
- [x] T003 Collect each engine's pose file suffix from the `pose_file` column in `compute_geometric_consensus()`. Build `pose_formats` dict and pass to `_pairwise_rmsd()`.
- [x] T004 Update `geometric_agreement` condition (line ~332): `True` when all pairwise entries have `rmsd <= cutoff OR (soft_alignment_used AND soft_rmsd <= cutoff)` and no hard failures (missing file / count mismatch > 2).

## Fix 2 — Redocking Validation SDF Parser

- [x] T005 [P] In `redocking_validation.py`, import `_extract_sdf_pose` from `post_docking_analysis.geometric_consensus`.
- [x] T006 [P] Add `pose_index: int = 1` parameter to `_parse_pdb_like_heavy_atoms()`. When suffix is `.sdf` or `.mol`, call `_extract_sdf_pose(file_path, pose_index)` and return its `(coords, elements, error)` directly instead of `_parse_sdf_heavy_atoms(lines)`.
- [x] T007 [P] Remove `_parse_sdf_heavy_atoms()` from `redocking_validation.py` (now superseded by shared implementation). Update `_compute_pose_rmsd()` to pass `pose_index=1` to `_parse_pdb_like_heavy_atoms()`.
- [x] T008 [P] Add comment to `_extract_sdf_pose` in `geometric_consensus.py` marking it as semi-public (used by `redocking_validation.py`).

## Fix 3 — `favorite_engine_continue` Engine Membership Guard

- [x] T009 [P] In `_legacy_multi_engine_pipeline_impl.py`, after the `if not favorite_engine` check at line ~1842, add validation: if `self.engines_in_scope` is non-empty and `favorite_engine not in self.engines_in_scope`, raise `ValueError` with a clear message listing valid engines.

## Tests

- [x] T010 Unit test Fix 1 — soft alignment path: build two fake pose files (`.sdf` and `.pdbqt`) for the same ligand with 44 vs 43 heavy atoms at nearly identical coordinates. Call `compute_geometric_consensus()`. Assert `geometric_agreement=True` and `soft_alignment_used=True` in the pairwise JSON.
- [x] T011 [P] Unit test Fix 1 — hard failure path: same setup but with 44 vs 30 heavy atoms (count delta > 2). Assert `geometric_agreement=False` and `soft_alignment_used=False`.
- [x] T012 [P] Unit test Fix 2 — multi-pose GNINA SDF: write a 2-pose SDF (88 atom lines split by `$$$$`) to a temp file. Call `_parse_pdb_like_heavy_atoms(path, pose_index=1)`. Assert returned coord count equals single-pose heavy atom count (44), not 88.
- [x] T013 [P] Unit test Fix 2 — redocking end-to-end: run `run_redocking_validation()` with a best_by_engine row pointing to a multi-pose GNINA SDF and a reference PDBQT. Assert `redocking_classification` is `pass`, `warn`, or `fail` (not `not_evaluable`).
- [x] T014 [P] Unit test Fix 3: instantiate `MultiEngineAnalysisPipeline` with `engines_in_scope=["vina"]` and `favorite_engine="gnina"`. Call `run()` in `favorite_engine_continue` mode. Assert `ValueError` is raised with `gnina` and `vina` in the message.
- [x] T015 [P] Unit test Fix 3 — empty scope does not raise: same setup but `engines_in_scope=[]`. Assert `ValueError` is NOT raised at the guard (may fail downstream for other reasons, but the membership guard must not fire).
