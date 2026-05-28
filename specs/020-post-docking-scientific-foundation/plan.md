# Implementation Plan: Post-Docking Scientific Foundation

**Branch**: `020-post-docking-scientific-foundation` | **Date**: 2026-03-31 | **Spec**: `/specs/020-post-docking-scientific-foundation/spec.md`  
**Input**: Scientific-correctness spec from `/specs/020-post-docking-scientific-foundation/spec.md`

## Summary

This plan hardens post-docking analysis so outputs are scientifically defensible and reproducible.  
The critical track is correctness-first: fix ranking direction, mapping integrity, statistical validity, and structure integrity before adding higher-level consensus and classification refinements.

## Current State Snapshot

Completed foundations already in code:
- `unified_pipeline.py` exists and execution paths are routed through unified flow.
- `structure_processor.py` retired; structure validation now uses `complex_validation.py`.
- Consensus normalization direction fixed for min-max and z-score.
- Pairlist matching is strict/equality-based with unmatched and empty log reporting.
- Correlation analyzer enforces minimum N, rank normalization, and BH/FDR q-values.
- Single-engine QC handling avoids universal low-agreement downgrades.
- Atomic QC downgrade behavior is in place.
- Deterministic tie-break logic exists in pose selection.
- Dashboard/report newline contract fixed.
- Validation gate artifacts + policy fallback (`reference_anchor -> target_percentile`) wired in comparative outputs.
- RMSD defaults are now per-complex only for interactive/CLI pipeline paths.

Open gaps to close:
- Full retirement contract for legacy entrypoints (`simplified_pipeline.py`, `multi_engine_pipeline.py`) is not complete yet.
- Scientific gate behavior needs broader regression coverage across comparative modes.
- Some legacy modules still contain fallback behavior that should be normalized under the unified contract.
- Final docs and migration notes for users need completion.

## Technical Scope

Primary modules:
- `post_docking_analysis/unified_pipeline.py`
- `post_docking_analysis/multi_engine_pipeline.py`
- `post_docking_analysis/consensus.py`
- `post_docking_analysis/generate_scores_csv.py`
- `post_docking_analysis/correlation_analyzer.py`
- `post_docking_analysis/pose_extractor.py`
- `post_docking_analysis/redocking_validation.py`
- `post_docking_analysis/report_generator.py`
- `workflow/interactive.py`
- `workflow/cli.py`

Test harness:
- `test/test_dockforge_smoke.py`
- `test_pipeline.py`

## Phase Plan

### Phase 1: Pipeline Unification

Goal:
- Make unified pipeline the only user-facing post-docking entrypoint.

Work:
- Keep wrappers for compatibility but enforce delegation-only behavior.
- Remove any divergent logic paths between simplified/multi-engine legacy scripts.
- Maintain output-path compatibility.

Exit criteria:
- One orchestrator contract for single-engine and multi-engine runs.
- Legacy modules no longer contain primary analysis logic.

### Phase 2: Critical Correctness Fixes

Goal:
- Eliminate known bugs that corrupt scientific meaning.

Work:
- Confirm and lock normalization direction behavior across all methods.
- Keep strict pairlist mapping and unmatched/empty log artifacts.
- Keep deterministic pose tie-break behavior and residue identity preservation.
- Keep report formatting/data-contract correctness.

Exit criteria:
- Deterministic and directionally correct ranking outputs in synthetic checks.

### Phase 3: Validation Gate

Goal:
- Prevent unsafe reference-anchored classification.

Work:
- Emit and persist `validation_gate_status.json`.
- Enforce policy fallback when gate disallows reference anchors.
- Surface gate state in consolidated run context and canonical mirrors.

Exit criteria:
- `reference_anchor` never applies when gate is not validated.

### Phase 4: Scoring Correctness

Goal:
- Ensure all ranking and classes are target-aware and statistically meaningful.

Work:
- Keep per-target normalization for cross-engine comparison.
- Keep N<3 suppression (`classification_basis=insufficient_n`).
- Maintain single-engine-safe class/QC behavior.

Exit criteria:
- No false confidence for sparse targets or single-engine runs.

### Phase 5: Statistical Rigor

Goal:
- Make correlation claims publication-safe.

Work:
- Keep minimum-N guardrails.
- Keep BH/FDR correction and explicit q-values.
- Keep low-power flags for marginal sample sizes.

Exit criteria:
- Correlation tables always include corrected significance context.

### Phase 6: Classification Calibration

Goal:
- Maintain defensible hit classes under reference or percentile policies.

Work:
- Retain validated reference-anchor behavior.
- Preserve fallback behavior and explicit classification basis labels.
- Keep docking/QC/ADMET channels separate and auditable.

Exit criteria:
- Hit class semantics are explicit and traceable in all reports.

### Phase 7: DockBox-Style Consensus Hardening

Goal:
- Keep geometric agreement as first-class consensus evidence.

Work:
- Preserve `dockbox_geometric` default.
- Keep weighted/strict/favorite modes as explicit alternatives.
- Maintain explainability artifacts for mode, normalization, and class distributions.

Exit criteria:
- Consensus mode and rationale are explicit in artifacts and summaries.

## Migration Strategy

1. Preserve current output contracts while shifting orchestration to unified path.
2. Keep compatibility wrappers in place during transition.
3. Add regression gates for every fixed scientific defect before deleting legacy paths.
4. Publish migration notes after wrappers become pass-through-only.

## Validation Strategy

Required commands:
- `python -m compileall workflow/interactive.py workflow/cli.py post_docking_analysis/multi_engine_pipeline.py post_docking_analysis/redocking_validation.py post_docking_analysis/cli.py post_docking_analysis/simplified_cli.py post_docking_analysis/simplified_pipeline.py`
- `python test_pipeline.py`
- `python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix`

Targeted smoke contracts:
- normalization direction + single-engine QC
- strict pairlist matching
- minimum-N + BH/FDR correction
- pose parser/tie-break determinism
- validation gate + effective policy fallback
- run tracking/report contract integrity

## Risks and Mitigations

- Risk: Legacy and unified paths drift again.
  Mitigation: Force wrappers to delegate; add smoke assertion for active entrypoint path.

- Risk: Users perceive stricter stats as “worse” outcomes.
  Mitigation: Keep raw p-values alongside q-values and label low-power contexts.

- Risk: Gate fallback surprises users expecting reference anchors.
  Mitigation: Write explicit gate reason and effective policy into summaries.

- Risk: Existing downstream scripts expect prior file placements.
  Mitigation: Keep canonical mirrors and compatibility outputs during transition.
