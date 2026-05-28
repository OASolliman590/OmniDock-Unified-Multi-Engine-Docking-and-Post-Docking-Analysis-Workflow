# Feature Specification: Top Pose Ligand Performance Atlas

**Feature Branch**: `019-top-pose-ligand-performance`  
**Created**: 2026-03-31  
**Status**: Draft  
**Input**: User request to make top-performing poses per ligand a first-class, queryable, and report-ready analysis product.

## Summary

Current pipeline logic already computes best poses in multiple places (`best_poses.csv`, per-engine best-per-tag, and best-per-ligand summary tables), but the behavior is spread across modules and outputs are not unified into a single ligand-centric decision artifact.

This feature adds a dedicated **Top Pose Ligand Performance Atlas** that standardizes how the best pose for each ligand is selected, compared across proteins/engines, confidence-scored, and exported.

## Problem

- “Best pose” exists, but not as one authoritative ligand-centric endpoint.
- Per-ligand top performance is not consistently linked to consensus confidence, QC status, ADMET status, and external biology signals.
- Users cannot quickly answer: “For ligand X, what is the best pose overall, best per target, and how stable is that decision?”

## Goals

1. Create a canonical ligand-centric top-pose artifact with reproducible rules.
2. Support two views:
- best pose per ligand per protein
- global best pose per ligand across all proteins
3. Attach explainability fields (engine agreement, rank stability, QC/ADMET flags, optional biology evidence).
4. Write outputs to predictable canonical folders and include them in run tracking indexes.

## Non-Goals

- Replacing existing allscore consensus implementation.
- Changing docking engine scoring formulas.
- Introducing new rescoring models (handled in separate rescoring features).

## User Stories

### US1: Ligand-Centric Selection
As a scientist, I want one table that shows each ligand’s best pose per target and globally so I can prioritize compounds quickly.

### US2: Evidence-Aware Ranking
As a reviewer, I want top-pose rows to include confidence/QC/ADMET/biology context so ranking decisions are explainable.

### US3: Stable Re-Runs
As an operator, I want top-pose outputs to be deterministic and resumable, with run metadata proving what changed.

## Functional Requirements

- **FR-001**: System MUST compute `top_pose_per_ligand_per_protein` from normalized comparative scores.
- **FR-002**: System MUST compute `top_pose_per_ligand_global` across proteins using configured aggregation.
- **FR-003**: Selection policy MUST be configurable (`best_affinity`, `best_consensus`, `hybrid`).
- **FR-004**: Tie-breaking MUST be deterministic (consensus score, then affinity, then tag lexical).
- **FR-005**: Output rows MUST include pose identity fields: engine, tag, protein, ligand, site_id, pose index, pose file path.
- **FR-006**: Output rows MUST include interpretability fields: docking quality class, qc_status, admet_status.
- **FR-007**: When available, output rows MUST include biology correlation fields (`bio_*`) and mapping status.
- **FR-008**: Output rows MUST include confidence metadata (`engine_support_count`, optional bootstrap/rank-stability metrics when available).
- **FR-009**: System MUST export canonical files under `5-Analysis/top_pose_ligand_performance/`.
- **FR-010**: System MUST emit both CSV and JSON summary artifacts.
- **FR-011**: System MUST register these outputs in run tracking (`run_tracking/outputs_index.csv`).
- **FR-012**: System MUST support scope-only execution (`top_pose_only`) without rerunning full downstream visualization.

## Output Contract

Primary artifacts:
- `5-Analysis/top_pose_ligand_performance/top_pose_per_ligand_per_protein.csv`
- `5-Analysis/top_pose_ligand_performance/top_pose_per_ligand_global.csv`
- `5-Analysis/top_pose_ligand_performance/ligand_performance_summary.csv`
- `5-Analysis/top_pose_ligand_performance/top_pose_selection_manifest.json`

Optional artifacts:
- `5-Analysis/top_pose_ligand_performance/top_pose_confidence_metrics.csv`

## Acceptance Criteria

1. For a multi-engine project, each ligand appears exactly once in global top-pose table.
2. For each `(ligand, protein)` pair with valid scores, exactly one top pose is selected.
3. Rerunning with identical inputs produces identical top-pose rows and order.
4. Output files appear in canonical location and are listed by run tracking index.
5. CLI/interactive can run this scope without requiring full visualization stages.

## Risks / Open Questions

- Whether global ligand aggregation should default to `best target` or `mean of top-k targets`.
- Whether confidence metrics should be required when bootstrap inputs are sparse.
- Potential confusion between “best by affinity” vs “best by consensus score”; requires explicit metadata in output.

## Recommended Defaults

- Selection policy: `hybrid` (consensus first, affinity second).
- Global aggregation: `best_target` (scientifically conservative for hit triage).
- Include QC/ADMET gates as annotations, not hard exclusions by default.
