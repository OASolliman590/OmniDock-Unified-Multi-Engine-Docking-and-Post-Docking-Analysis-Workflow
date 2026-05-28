# Spec 017: DockBox-Style Consensus and Post-Docking Maturation

## Goal
Strengthen scientific hit selection by adding a consensus decision layer to multi-engine analysis and mature post-docking outputs into a predictable, review-friendly structure.

## User Outcomes
- Comparative analysis supports three consensus policies:
  - `weighted_hybrid`
  - `strict_consensus`
  - `favorite_guardrails`
- Rescoring shortlist selection supports both:
  - `top_n_per_protein`
  - `top_n_global`
- Comparative runs emit explicit explainability and candidate artifacts:
  - `consensus_inputs_all_engines.csv`
  - `consensus_ranked_hits.csv`
  - `rescoring_candidates.csv`
  - `consensus_explainability.json`
- Post-docking outputs expose a canonical analysis view:
  - `analysis/raw`
  - `analysis/processed`
  - `analysis/figures/2d`
  - `analysis/figures/3d`
  - `analysis/reports`

## Scope
- In-scope:
  - Comparative analysis ranking logic and rerun-candidate gating.
  - Interactive/CLI configuration surfaces for consensus and rescoring.
  - Output consolidation for comparative analysis and stage-level runs.
- Out-of-scope:
  - Replacing docking engines or changing scoring formulas produced by engines.
  - Forced migration of legacy historical analysis outputs.

## Acceptance Criteria
- Comparative analysis runs with any of the 3 consensus modes without code changes.
- Rescoring candidates are generated for both scopes with configurable top-N.
- Exhaustive rerun manifests carry consensus metadata and are filtered by consensus shortlist when available.
- Canonical analysis view is materialized and populated after successful analysis runs.
