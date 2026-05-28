# Tasks 017: DockBox-Style Consensus and Post-Docking Maturation

## Core
- [x] Add `post_docking_analysis/consensus.py` with consensus/rescoring utilities.
- [x] Thread consensus settings through comparative analysis execution path.
- [x] Generate consensus artifacts and explainability outputs.
- [x] Apply consensus shortlist gating to exhaustive rerun manifest generation.

## UX
- [x] Add interactive prompts for consensus mode, rescoring scope, and rescoring top-N.
- [x] Clarify stage-level engine vs favorite-engine continuation prompts.
- [x] Add explicit `Exit wizard` action to analysis control panel.
- [x] Add CLI flags for consensus/rescoring in both workflow CLI and post-docking CLI.

## Output Maturation
- [x] Materialize canonical analysis view (`analysis/raw`, `processed`, `figures/2d`, `figures/3d`, `reports`) for comparative runs.
- [x] Trigger artifact consolidation after successful stage-target runs.

## Verification
- [x] Run workflow CLI smoke checks for comparative analysis with each consensus mode.
- [x] Run post-docking CLI smoke checks for rescoring scope variants.
- [x] Validate canonical analysis layout on a full-session run in maturation workspace.
