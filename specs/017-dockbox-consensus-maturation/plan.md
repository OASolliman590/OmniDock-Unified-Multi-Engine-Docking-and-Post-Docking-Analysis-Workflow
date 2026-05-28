# Plan 017: DockBox-Style Consensus and Post-Docking Maturation

## Design Decisions
- Default policy:
  - Consensus: `weighted_hybrid`
  - Rescoring scope: `top_n_per_protein`
- Consensus is built from best-pose-per-engine rows (`engine`, `tag`), not raw all-pose rows.
- Rerun promotion remains engine-aware but is gated by consensus shortlist when available.
- Output maturation keeps existing files and adds a canonical analysis index/view for faster navigation.

## Implementation Strategy
1. Add consensus utility module with:
   - mode normalization
   - consensus ranking table generation
   - rescoring candidate selection
   - explainability payload generation
2. Integrate consensus/rescoring into multi-engine comparative reporting and rerun manifest generation.
3. Extend interactive and CLI surfaces to collect/pass consensus configuration.
4. Consolidate outputs into canonical `analysis/{raw,processed,figures,reports}` view.
5. Ensure stage-target analysis paths call output consolidation when successful.

## Validation Strategy
- Static validation via import/compile checks.
- CLI dry validation:
  - comparative analysis arg parsing with new options.
- Runtime smoke:
  - comparative mode creates consensus/rescoring artifacts.
  - canonical analysis view directories are created and populated.
