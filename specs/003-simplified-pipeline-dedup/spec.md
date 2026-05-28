# Feature Specification: Simplified Pipeline Deduplication

## Context

The simplified post-docking pipeline currently repeats the same logic in multiple stages for:
- selecting best poses from `scores_df`
- checking and listing generated complex PDB files

This duplication increases maintenance cost and creates drift risk when behavior is updated in one stage but not another.

## Goal

Unify duplicated helper logic in `post_docking_analysis/simplified_pipeline.py` without changing runtime behavior.

## Scope

- Add one helper for best-pose selection from scores.
- Add one helper for discovering generated complex PDB files.
- Refactor existing stages to use these helpers.

## Out of Scope

- Functional changes to scoring, visualization generation, or external tool invocation.
- CLI changes.
- New dependencies.

## Acceptance Criteria

1. The following stages use centralized helper methods instead of duplicated inline logic:
   - RMSD analysis
   - reports generation
   - py3Dmol visualization stage
   - ProLIF stage
   - LigPlot stage
   - PandaMap stage
2. Pipeline behavior remains unchanged for successful and missing-input paths.
3. `python -m py_compile` passes for modified files.
4. A simplified CLI smoke run completes successfully in the `pdb-prepare-wizard` conda env.
