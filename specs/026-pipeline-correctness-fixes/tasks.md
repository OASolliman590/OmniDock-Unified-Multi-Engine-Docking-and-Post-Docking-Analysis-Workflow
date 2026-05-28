# Tasks: Pipeline Correctness Fixes

**Input**: `/specs/026-pipeline-correctness-fixes/spec.md`
**Scope**: 5 bug categories, 7 discrete code changes, no new artifacts, no new API surface beyond one optional parameter.

## Format: `[ID] [P?] Description`
- **[P]**: parallelizable (no dependency on other open tasks)

## Fix 1 — Top-Pose Selector Direction

- [x] T001 Add `_CONSENSUS_SCORE_ASCENDING_MODES` frozenset and `_consensus_score_ascending(mode: str) -> bool` helper to `top_pose_selector.py`. Rank-pct modes (`strict_consensus`, `weighted_hybrid`, `""`) return `True`; agreement-based modes (`dockbox_geometric`, `favorite_guardrails`) return `False`.
- [x] T002 Add `consensus_mode: str = ""` parameter to `build_top_pose_atlas()`. Pass it through to sort calls. Replace the three hardcoded `ascending=[..., True, ...]` entries for `_sort_primary` with `_consensus_score_ascending(consensus_mode)`.
- [x] T003 Fix `best_consensus_score` summary aggregation: after the main `.agg()` call, recompute `best_consensus_score` using `min` when `_consensus_score_ascending(consensus_mode)` else `max`.
- [x] T004 Propagate `consensus_mode=self.consensus_mode` from `_dag_compute_top_pose_atlas_node()` in `_legacy_multi_engine_pipeline_impl.py` to the `build_top_pose_atlas()` call.

## Fix 2 — Reference Baselines Wiring

- [x] T005 Add `"reference_baselines"` path key to `_dag_artifact_paths()` in `_legacy_multi_engine_pipeline_impl.py`, pointing to `<validation_gate_dir>/reference_baselines.csv`.
- [x] T006 In `_dag_compute_classified_hits_node()`, load `paths["reference_baselines"]` if it exists and pass as `reference_baselines=` to `classify_hits_target_aware()`. Pass empty DataFrame if file absent.

## Fix 3 — Engine Detection Cache Override

- [x] T007 [P] In `engine_detector.py` cache-hit early-return block (lines ~199–202), apply `override_engine` to the cached dict before returning: set `detection_override`, `override_reason`, `routed_to_engine`, `routing_decision` fields when override is non-empty.

## Fix 4 — RMSD Scope Normalizer

- [x] T008 [P] In `_legacy_simplified_pipeline_impl.py`, define `_VALID_RMSD_SCOPES = ("per_complex", "per_protein", "global")`. Update `_normalize_rmsd_scopes()` to accept all valid tokens (guard against `_VALID_RMSD_SCOPES` instead of `_DEFAULT_RMSD_SCOPES`). Update `"all"` / `"*"` expansion to return `_VALID_RMSD_SCOPES`.

## Fix 5 — `strict_consensus` Single-Engine Warning

- [x] T009 [P] In `consensus.py`, when `mode == "strict_consensus"` and `single_engine=True`, emit `logger.warning(...)` stating that strict consensus filtering is skipped in single-engine mode. Add one-line note to function docstring.

## Tests

- [x] T010 Unit test Fix 1: call `build_top_pose_atlas()` with `consensus_mode="dockbox_geometric"` and two rows where row A has higher `consensus_score` and worse `affinity_kcal_mol` than row B. Assert row A is selected (highest consensus_score wins).
- [x] T011 [P] Unit test Fix 1 inverse: same setup with `consensus_mode="weighted_hybrid"`. Assert row with lower `consensus_score` is selected.
- [x] T012 [P] Unit test Fix 2: create a `reference_baselines.csv` on disk with one reference row; run `_dag_compute_classified_hits_node()` via the DAG; assert the loaded baselines were passed (check classified output contains `reference_anchor`-derived class if policy is set).
- [x] T013 [P] Unit test Fix 3: call `detect_engines()` with a pre-populated cached manifest and `override_engine="vina"`. Assert `routed_to_engine == "vina"` and `detection_override == True` without `redetect=True`.
- [x] T014 [P] Unit test Fix 4: call `_normalize_rmsd_scopes("per_protein,global")` and assert result is `("per_protein", "global")`. Call with `"all"` and assert result is `("per_complex", "per_protein", "global")`.
- [x] T015 [P] Unit test Fix 5: call `compute_consensus_ranking()` with `mode="strict_consensus"` and a single-engine frame; assert a warning was logged and all rows are retained (no filtering).
