# Tasks: Top Pose Ligand Performance Atlas

**Input**: `/specs/019-top-pose-ligand-performance/spec.md`, `/specs/019-top-pose-ligand-performance/plan.md`  
**Prerequisites**: comparative allscore outputs available

## Format: `[ID] [P?] Description`

- **[P]**: parallelizable

## Phase 1: Selection Core

- [x] T001 Create `post_docking_analysis/top_pose_selector.py` with policy-driven selectors (`best_affinity`, `best_consensus`, `hybrid`).
- [x] T002 [P] Implement deterministic tie-break utility (consensus score, affinity, lexical tag).
- [x] T003 [P] Implement per-protein ligand top-pose reducer.
- [x] T004 Implement global ligand top-pose reducer.
- [x] T005 Add confidence metadata fields (`engine_support_count`, `selection_policy`, `tie_break_rule`).

## Phase 2: Pipeline Integration

- [x] T006 Integrate selector into `MultiEngineAnalysisPipeline._write_comparative_reports`.
- [x] T007 Write session-local outputs under `output_dir/top_pose_ligand_performance/`.
- [x] T008 Write canonical mirrored outputs under `post_docking_root(project)/5-Analysis/top_pose_ligand_performance/`.
- [x] T009 Generate `top_pose_selection_manifest.json` with selection policy, counts, timestamps, and provenance context.
- [x] T010 Register top-pose outputs in run tracking index (`run_tracking/outputs_index.csv`).

## Phase 3: Scope + UX Hooks

- [x] T011 Add `top_pose_only` analysis scope in `post_docking_analysis/multi_engine_pipeline.py` and CLI validators.
- [x] T012 [P] Wire scope prompt/flag in `workflow/cli.py` and `workflow/interactive.py`.
- [x] T013 Ensure scope skips heavy visual/interaction stages while still generating top-pose artifacts.

## Phase 4: Validation

- [x] T014 Add unit tests for selector determinism and tie-break behavior.
- [x] T015 Add smoke test asserting one global row per ligand and one per `(ligand, protein)` pair.
- [x] T016 Add smoke test asserting top-pose files are present in run tracking index.
- [x] T017 Update docs (`POST_DOCKING_ANALYSIS_GUIDE.md`) with top-pose atlas usage and output locations.

## Exit Criteria

- [x] E001 All top-pose artifacts are produced in both session and canonical mirrors.
- [x] E002 Results are deterministic across reruns with unchanged inputs.
- [x] E003 Run tracking includes top-pose artifacts and step status.
- [x] E004 Scope-only execution (`top_pose_only`) works without full visualization pipeline.
