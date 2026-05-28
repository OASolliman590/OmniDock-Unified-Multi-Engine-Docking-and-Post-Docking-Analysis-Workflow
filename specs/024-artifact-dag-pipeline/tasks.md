# Tasks: Artifact DAG Pipeline

**Input**: `/specs/024-artifact-dag-pipeline/spec.md`
**Scope**: Replace 16-step linear pipeline with a declarative artifact DAG — parallel branches, hash-based caching, per-node failure isolation, single/multi-engine strategy at consensus node only.

## Format: `[ID] [P?] Description`
- **[P]**: parallelizable

## Phase 1: ArtifactGraph Executor

- [x] T001 Create `post_docking_analysis/artifact_graph.py` with `ArtifactNode` dataclass (name, inputs, outputs, compute, optional).
- [x] T002 Implement `ArtifactGraph.register()` and `ArtifactGraph.request(artifact, force=False)`.
- [x] T003 Implement topological sort → tier grouping (nodes in same tier have no inter-dependencies).
- [x] T004 Implement tier-based `ThreadPoolExecutor` execution; advance to next tier when current tier completes.
- [x] T005 Implement per-node status tracking: `pending | running | completed | completed_with_warnings | failed | skipped_optional | blocked_by_failure | cache_hit`.
- [x] T006 Implement hash-based cache: `cache_key = sha256(sorted input file mtimes)`; store in `4-Working/dag_cache.json`.
- [x] T007 Implement `_is_cached()`: return True if all outputs exist and cache_key matches stored key.
- [x] T008 Implement optional node contract: `optional=True` failure marks node `skipped_optional`; dependents not blocked.

## Phase 2: Node Registration

- [x] T009 Register N-001 `raw_scores`: inputs=engine outputs + pairlist, outputs=`raw_scores.csv`.
- [x] T010 Register N-002 `normalized_scores`: inputs=`raw_scores.csv`, outputs=`normalized_scores.csv`.
- [x] T011 Register N-003 `consensus_ranked`: inputs=`normalized_scores.csv` + `engine_scope_config`, outputs=`consensus_ranked.csv` + `engine_agreement.csv`.
- [x] T012 Register validation gate as side-input to N-003 (not a blocking dependency).
- [x] T013 Register N-004 `classified_hits`: inputs=`consensus_ranked.csv` + `validation_gate.json`, outputs=`classified_hits.csv`.
- [x] T014 Register `top_pose_atlas` node: inputs=`classified_hits.csv`, outputs=`top_pose_per_ligand_global.csv` + `ligand_performance_summary.csv`.
- [x] T015 Register N-005 `complexes`: inputs=`classified_hits.csv` + `top_pose_atlas` + pose files, outputs=`complexes/`.
- [x] T016 Register N-006 interaction nodes as `optional=True`: prolif, pandamap, poseview, pymol; all input=`complexes/`.
- [x] T017 Register parallel analysis nodes: `polypharmacology`, `comparative`, `biology_correlation`; all input=`classified_hits.csv`.
- [x] T018 Register N-007 `reports`: inputs=all N-006 nodes + parallel analysis nodes + `classified_hits.csv`.
- [x] T019 Wire existing pipeline functions as `compute` callables in each node — no logic changes.

## Phase 3: Consensus Strategy Injection

- [x] T020 Implement strategy loader: read `engine_scope_config` → select `multi | single/gnina | single/vina | single/smina` strategy.
- [x] T021 Implement multi-engine strategy: DockBox geometric consensus + cross-engine rank aggregation (spec 017).
- [x] T022 [P] Implement GNINA solo strategy: `cnn_affinity` primary, `cnn_confidence` annotation.
- [x] T023 [P] Implement Vina solo strategy: `vina_affinity` primary, `rmsd_lb` pose diversity.
- [x] T024 [P] Implement Smina solo strategy: `vina_affinity` primary, scoring provenance.
- [x] T025 Assert `consensus_ranked.csv` schema is identical for all four strategies (column names, types).

## Phase 4: Parallel Interaction Nodes

- [x] T026 Convert prolif node to per-complex sub-tasks; run via inner `ThreadPoolExecutor`.
- [x] T027 [P] Convert pandamap node to per-complex sub-tasks; run via inner `ThreadPoolExecutor`.
- [x] T028 [P] Convert poseview node to per-complex sub-tasks; skip if binary absent.
- [x] T029 [P] Convert pymol node to per-complex sub-tasks; mark `optional=True`; skip if binary absent.
- [x] T030 Per-complex failure in any interaction node: log warning, continue remaining complexes, set node status `completed_with_warnings`.

## Phase 5: Scope System

- [x] T031 Implement `SCOPE_MAP` dict mapping scope names to artifact names.
- [x] T032 Wire `--scope` CLI flag to `graph.request(SCOPE_MAP[scope])`.
- [x] T033 Validate all scope values produce correct subgraph (no extra nodes executed).
- [x] T034 Add `--force` flag to bypass cache and re-run all nodes in requested subgraph.

## Phase 6: Execution Report + Cache

- [x] T035 Emit `dag_execution_report.json` to `4-Working/metadata/` after every run.
- [x] T036 Include per-node status, duration, cache_hit in report.
- [x] T037 Compute and include `parallel_time_saved_s` (sum of parallel branch durations minus wall clock).
- [x] T038 Update `dag_cache.json` with new cache keys for completed nodes.

## Phase 7: Migration + Tests

- [x] T039 Replace step-sequence loop in `unified_pipeline.py` with `graph.request(scope_artifact)`.
- [x] T040 Remove hard-coded step ordering; confirm no references to step numbers in pipeline logic.
- [x] T041 Smoke test: full 3-engine run; assert `parallel_time_saved_s > 0` in execution report.
- [x] T042 Smoke test: re-run unchanged → all nodes except `reports` show `cache_hit: true`.
- [x] T043 Smoke test: change normalization → N-001 cache_hit, N-002 miss, N-003 through N-007 re-executed.
- [x] T044 Smoke test: ProLIF failure for one complex → pandamap, polypharmacology, reports all complete normally.
- [x] T045 Smoke test: `--scope comparison_only` → no pose files or complexes written.
- [x] T046 Smoke test: GNINA solo → `consensus_ranked.csv` `primary_score_name = cnn_affinity`; same schema as multi-engine.
- [x] T047 Smoke test: `--scope report_only` with all upstream cached → completes in < 5 seconds.
- [x] T048 Smoke test: `dag_execution_report.json` present and valid JSON after every run.
