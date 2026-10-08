# Decision needed — divergence between this branch and `main`

**Status:** open; Scientific Lead decision required. No merge or rebase has been performed.

## Facts (2026-10-08)

- This branch and `main` diverge at `cfe3305`.
- `main` has two commits this branch lacks:
  - `7fc130a` "Correct scientific preparation, docking execution, and post-docking validation" (90 files, about +7.5k/−2.4k lines). It adds an `audit/` package, strict receptor/ligand preparation, a job contract, reference chemistry, workflow locking, and changes to consensus, geometric consensus, redocking validation, pose selection and the artifact graph.
  - `8e36027`: dependency declarations.
- This branch has Spec 031 (accepted methods) and Spec 032 (MD-input export) work that `main` lacks.
- A trial merge (`git merge-tree`, no write to any branch) reports content conflicts in 15 files, including `post_docking_analysis/consensus.py`, `geometric_consensus.py`, `redocking_validation.py`, `top_pose_selector.py`, `multi_engine_pipeline_impl.py`, `artifact_graph.py`, `core_pipeline.py` and `docking/project_layout.py`.

## Why this is a scientific decision

The two lines implement overlapping scientific methods differently, for example RMSD/redocking, consensus and pose selection. Resolving the conflicts means choosing which method is authoritative, so it cannot be treated as a mechanical merge.

## Spec 033's bounded use of `main`

Spec 033 ports only the strict receptor-preparation module from `7fc130a` (R2a), because it directly fixes P2 and does not overlap Spec 031/032 methods. Everything else in `7fc130a` stays unported until this decision is made.

## Options for the Scientific Lead

1. Keep this branch authoritative and cherry-pick selected `main` changes under new specs.
2. Make `main` authoritative and re-apply Spec 031/032 on top under review.
3. Perform a full reconciliation spec that compares each conflicting method side by side.
