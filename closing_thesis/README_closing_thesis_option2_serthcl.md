# Closing Thesis Similarity Analysis (Option 2 with SertHCl)

This bundle includes:
- Primary enzyme-context comparisons vs biology and docking comparators
- SertHCl anchor comparisons for the same query compounds and enzyme contexts

Query compounds: MS1, MS2, MS4, MS6, PABA, MSPS
Enzyme contexts: EGFR, VEGFR-II, CDK-2, TrkA, COX-II, COX-I

Coverage summary:
{
  "expected_plan_rows": 108,
  "actual_plan_rows": 108,
  "expected_metric_rows": 324,
  "actual_metric_rows": 324,
  "actual_pair_mean_rows": 108,
  "actual_mcs_rows": 108,
  "missing_metric_rows": 0,
  "missing_pair_mean_rows": 0,
  "missing_mcs_rows": 0
}

## Bundle layout

Manuscripts (moved from the repository root; bytes unchanged):

- `closing_thesis/manuscripts/Sertraline_Docking_Methods_Numbered_With_Citation_Placeholders.docx`
- `closing_thesis/manuscripts/Sertraline_Docking_Methods_Omni_DockForge.docx`
- `closing_thesis/manuscripts/Sertraline_Docking_Methods_Thesis_Chapter.docx`

Figures live in `closing_thesis/figures/` (see `figures/FIGURE_INDEX.md`), including
the relocated `multiscale_overlapped_experimental_vs_triple_engine.{png,svg}` pair.

## Provenance limitation

This README previously pointed at `run_manifest_option2_serthcl.json` for full
provenance. That file is **not tracked** (`.gitignore` already excludes
`closing_thesis/run_manifest*.json` and `closing_thesis/*manifest*.json`) and is
not present at the fixed review base
`8d4434d3a33eb83b1e12cad82944b02c83270e47`. Coverage numbers, figures, and
manuscript claims in this bundle therefore require independent provenance
verification before reuse as scientific evidence.
