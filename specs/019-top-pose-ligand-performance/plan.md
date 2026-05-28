# Implementation Plan: Top Pose Ligand Performance Atlas

**Branch**: `019-top-pose-ligand-performance` | **Date**: 2026-03-31 | **Spec**: `/specs/019-top-pose-ligand-performance/spec.md`  
**Input**: Feature specification from `/specs/019-top-pose-ligand-performance/spec.md`

## Summary

Implement a first-class ligand-centric top-pose layer on top of existing comparative outputs. The plan reuses existing best-pose and consensus logic, then adds a dedicated, deterministic export contract with confidence annotations and run-tracking integration.

## Current-State Assessment

Already present:
- Best pose extraction per tag and per group in `multi_engine_pipeline.py`.
- Best-per-ligand summaries in simplified/hierarchical analysis paths.
- Consensus and class annotations (`docking_quality_class`, `qc_status`, `admet_status`).

Missing:
- One canonical ligand-focused artifact set in a stable location.
- Explicit selection policy contract and tie-break documentation.
- Scope-level entrypoint (`top_pose_only`) and dedicated manifest.

## Architecture Approach

1. **Selection Engine (new utility layer)**
- Module: `post_docking_analysis/top_pose_selector.py`
- Inputs: normalized comparative score table + consensus table
- Outputs: per-protein and global top-pose tables + summary metrics

2. **Pipeline Integration**
- Integrate into `MultiEngineAnalysisPipeline._write_comparative_reports(...)`
- Execute after consensus generation so class/QC/biology fields are available.
- Register files in run-tracking outputs index.

3. **Canonical Output Routing**
- Write to: `post_docking_root(project)/5-Analysis/top_pose_ligand_performance/` (canonical mirror)
- Mirror to session output (`output_dir/top_pose_ligand_performance/`) for immediate run context.

4. **Scope Support**
- Add `top_pose_only` to analysis scopes for fast reruns.
- Skip heavy visualization/interactions when this scope is selected.

## Data Contract

### Input Tables
- `best_pose_per_tag_by_engine.csv`
- `consensus_ranked_hits_with_classes.csv` (or biology-augmented variant)

### Output Tables
- `top_pose_per_ligand_per_protein.csv`
- `top_pose_per_ligand_global.csv`
- `ligand_performance_summary.csv`
- `top_pose_selection_manifest.json`

### Key Columns
- Identity: `engine`, `tag`, `protein`, `ligand`, `site_id`, `pose`, `pose_file`
- Scoring: `affinity_kcal_mol`, `consensus_score`, `normalized_affinity_score`
- Interpretation: `docking_quality_class`, `qc_status`, `admet_status`, `bio_*`
- Confidence: `engine_support_count`, `selection_policy`, `tie_break_rule`

## Execution Phases

### Phase 1: Selection Core
- Implement policy-driven selectors:
  - `best_affinity`
  - `best_consensus`
  - `hybrid`
- Implement deterministic tie-break rules.

### Phase 2: Pipeline + Scope Integration
- Add `top_pose_only` scope handling.
- Invoke selector in comparative report generation.
- Persist canonical + session-local outputs.

### Phase 3: Tracking + Validation
- Add manifest generation (`top_pose_selection_manifest.json`).
- Add smoke tests for deterministic output and row uniqueness.
- Validate run tracking index inclusion.

## Testing Strategy

1. Unit tests for selector tie-break behavior.
2. Smoke test over mini comparative dataset:
- one global row per ligand,
- one row per `(ligand, protein)` where scores exist,
- stable ordering across reruns.
3. Integration check that run tracking lists all top-pose output files.

## Risks and Mitigations

- **Risk**: Conflicting definitions of “top” across teams.  
  **Mitigation**: explicit `selection_policy` in manifest + output column.

- **Risk**: Sparse proteins skew global ligand score.  
  **Mitigation**: default `best_target` aggregation, add optional `mean_top_k` later.

- **Risk**: Compatibility with old analysis sessions.  
  **Mitigation**: additive outputs only; no destructive schema change.
