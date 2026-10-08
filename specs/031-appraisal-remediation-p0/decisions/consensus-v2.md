# Gate 2 Decision — Versioned Consensus Policy

**Decision ID:** `consensus_rank_geometry_qc_v2`  
**Status:** Approved  
**Approved by:** Scientific Lead  
**Approval date:** 2026-08-21  
**Implements:** SR-002 / T014

## Approved Ranking Policy

1. Consensus v2 ranks ligands only within the same biological target and uses equal-weight per-engine rank percentiles. Raw scores from different engines or targets are never treated as a common scale.
2. The complete-case ligand set for a target contains ligands with a valid primary score from every engine in that target's eligible engine set. Per-engine percentiles and the composite are calculated on that same ligand set.
3. For `N > 1`, engine percentile is `p_e = (N - average_rank_e) / (N - 1)`, where rank 1 is best. For `N = 1`, the engine percentile is `1.0`, but the result is annotated `single_engine` and is not called multi-engine consensus.
4. Composite score is `mean(p_e)` with equal engine weights. Higher is better. No missing-component renormalization or neutral imputation is allowed.
5. If an engine is absent for the entire target, it is disclosed and excluded from that target's eligible engine set. If an eligible engine is missing for one ligand, that ligand is `consensus_incomplete` and receives no v2 composite score.

## Engine Metrics and Pose Selection

- **Vina/Smina:** lowest affinity selects and ranks; more negative is better. Ties use lower pose index, then stable tag.
- **AutoDock4:** lowest reported binding/free energy selects and ranks; more negative is better. Ties use lower pose index, then stable tag.
- **GNINA:** highest `cnn_score` selects the representative pose. Ties use higher `cnn_affinity`, then lower Vina affinity, then lower pose index, then stable tag. The selected pose's `cnn_affinity` ranks compounds; higher is better.
- Missing required GNINA CNN metrics are disclosed. Consensus v2 does not silently substitute Vina affinity for a missing CNN metric; legacy modes retain their documented behavior.
- Equal raw-score ties use average rank before percentile conversion.

## Geometry Policy

Geometry is a QC signal and tie-breaker, not a proportional score component.

- When two composite scores tie and both geometry results are comparable, verified agreement precedes verified disagreement.
- If either geometry result is `not_comparable`, geometry is not used to break that tie; stable ligand/tag identity provides deterministic ordering.
- `not_comparable`, missing geometry, and proven disagreement remain distinct states.
- No value such as `0.5` is imputed for missing geometry.

## Scope and Legacy Boundary

- Preserve `dockbox_geometric` and existing rank modes under their legacy names for historical reproduction.
- New outputs use `consensus_rank_geometry_qc_v2` and explicitly state that higher composite score is better.
- Cross-target tables may display per-target normalized results but cannot apply raw cross-target score comparisons or use unrelated target rank inversions as validation criteria.
- Known-active/decoy validation, engine-specific weights, and CNN integration beyond the approved GNINA policy require separate scientific approval.

## Required Provenance

Record formula ID/version, score direction, engine metrics and directions, target and complete-case ligand universe, engines requested/eligible/included/missing, equal weights, raw and normalized components, pose-selection/tie rules, missing-data policy, geometry status, RMSD method/cutoff/version, exclusions, and final deterministic tie resolution.

## Acceptance Cases

- Synthetic within-target data identify the expected winner with all engine directions represented.
- Permuting input rows does not change pose selection, percentiles, or order.
- One missing ligand/engine score yields `consensus_incomplete`, not a renormalized score.
- Missing geometry never changes the numeric composite.
- Proven geometry can break an exact composite tie only when both results are comparable.
- Single-engine output is annotated and never represented as multi-engine agreement.
