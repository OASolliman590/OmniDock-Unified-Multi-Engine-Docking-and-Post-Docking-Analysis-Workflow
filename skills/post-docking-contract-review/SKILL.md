---
name: post-docking-contract-review
description: Review post-docking analysis stages, mandatory vs optional stage execution, score normalization, and consensus ranking.
---

# Post-Docking Contract Review Skill

## Trigger Conditions
Activate when executing post-docking analysis, consensus calculation, or report generation.

## Explicit Inputs and Outputs
- **Inputs:** Raw docking poses and logs, `post_docking_analysis/config/schema.yaml`.
- **Outputs:** Master summary CSV, consensus ranking tables, stage status records.

## Refusal-to-Guess Rules
- **NEVER** treat an optional tool skip (e.g. PoseView, PLIP) as an error, and never treat an execution error as a clean skip.
- **NEVER** produce a report that appears complete when mandatory stages failed.

## Verification Commands
- `python -m post_docking_analysis.cli --help`
