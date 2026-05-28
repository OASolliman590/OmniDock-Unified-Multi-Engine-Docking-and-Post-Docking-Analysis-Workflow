# Spec 025: Visualization Report Suite

## Problem

The post-docking analysis pipeline produces no cohesive analytical figure suite. Figures are fragmented across three legacy code paths (`_legacy_multi_engine_pipeline_impl.py`, `_legacy_simplified_pipeline_impl.py`, `hierarchical_analyzer.py`), many are duplicated under different names, several key publication-ready figures from the old simplified pipeline are absent from the new DAG pipeline, and no figures adapt to single-engine mode.

## Goal

Define and implement a canonical **visualization report suite** as a single `visualizations` DAG node — 14 figures in 6 semantic groups — that works correctly for both multi-engine and single-engine runs, reads exclusively from already-computed DAG artifacts, and requires no new scientific computation.

---

## Inputs (all from upstream DAG nodes)

| Artifact | Source node |
|---|---|
| `classified_hits.csv` | `classified_hits` |
| `consensus_ranked.csv` | `consensus_ranked` |
| `engine_agreement.csv` | `consensus_ranked` |
| `normalized_scores.csv` | `normalized_scores` |
| `engine_scope_config.json` | written at graph build time |
| `validation_gate.json` | `validation_gate` |
| `top_pose_per_ligand_global.csv` | `top_pose_atlas` |

## Output

All figures written to `5-Analysis/sessions/<id>/visualizations/`, organized into subdirectories by group.

```
visualizations/
  overview/
    score_distribution.png
  ranking/
    consensus_leaderboard.png
    affinity_per_ligand.png
  multi_target/
    full_affinity_heatmap.png
    series_affinity_heatmap.png
    top_ligand_profiles.png
    polypharmacology_delta_heatmap.png      ← only when refs exist
    hit_class_matrix.png
  per_protein/
    affinity_by_protein.png
    best_ligand_per_protein.png
    hit_class_distribution_per_protein.png
  engine_agreement/                          ← multi-engine only
    engine_rank_correlation.png
    engine_agreement_per_complex.png
    per_protein_cross_engine/               ← existing batch plots
    pose_diversity.png                      ← single/vina only
    cnn_confidence_distribution.png         ← single/gnina only
  validation/                               ← only when validation gate passed
    rmsd_distribution.png
    reference_vs_docked.png
    rmsd_per_complex.png
  figures_manifest.json
```

---

## Figure Definitions

### Group 1 — Score Overview

**`score_distribution.png`**
- Histogram of all best-pose affinities from `normalized_scores.csv`
- Blue bars = favorable (≤ 0 kcal/mol), red bars = unfavorable (> 0)
- Vertical lines: mean (red dashed), median (green dashed), zero (orange dotted)
- Annotation box: unfavorable count and percentage
- Multi-engine: one panel per engine (subplots); Single-engine: single panel, title shows engine name and primary score label
- Score column: `affinity_kcal_mol` (multi-engine) or `primary_ranking_score` (single-engine)

---

### Group 2 — Ligand Ranking

**`consensus_leaderboard.png`**
- Horizontal bar chart, one bar per ligand
- Sorted descending by polypharmacology score = number of proteins where ligand is Moderate or better
- Bar color: gradient from green (high score) to grey (low score)
- Annotation on each bar: `targets=N, better_ref=M` (M omitted when no references exist)
- Multi-engine: score from `classified_hits.csv` `docking_quality_class`
- Single-engine: same, using `primary_ranking_score`-derived class

**`affinity_per_ligand.png`**
- Vertical bar chart, one bar per ligand, sorted ascending by best affinity
- Bars colored by hit class: Strong=green, Moderate=yellow-green, Weak=orange, Inactive=red
- Horizontal threshold line at -7 kcal/mol
- X-axis: ligand short name; Y-axis: best affinity (kcal/mol)

---

### Group 3 — Multi-Target Coverage

**`full_affinity_heatmap.png`**
- Protein × ALL ligands matrix
- Left section: reference ligands — only the protein-matched reference cell is filled; unmatched cells are white (diagonal pattern)
- Right section: series ligands — all cells filled
- Cell values: best affinity annotated; color scale: diverging green (≤ −10) → white (0) → red (> 0)
- Proteins on y-axis sorted by mean series affinity; ligands on x-axis: references first, then series sorted alphabetically
- Both sections separated by a vertical spacer column
- Multi-engine: uses `consensus_score` if available, falls back to `affinity_kcal_mol`
- Single-engine: uses `primary_ranking_score`

**`series_affinity_heatmap.png`**
- Series ligands only (references excluded)
- Ligands on y-axis sorted descending by mean affinity across all proteins (best overall binder at top)
- Proteins on x-axis
- Same color scale as full heatmap
- Cell values annotated

**`top_ligand_profiles.png`**
- Line plot: x-axis = protein targets, y-axis = best affinity (kcal/mol)
- One line per top-N ligand (N = min(10, total ligands))
- Selection: top N by polypharmacology score (same as leaderboard); tie-break by mean affinity
- Red dashed threshold line at −7 kcal/mol; pink dotted line at 0
- Legend inside plot; markers at each protein point
- Proteins sorted by mean affinity across top-N ligands (worst target on right)

**`polypharmacology_delta_heatmap.png`** *(only when validation gate has references)*
- Ligand × protein matrix of (series_affinity − reference_affinity)
- Negative = series binds better than reference at that target
- Color: red = positive delta (worse than ref), green = negative (better than ref)
- Cell values annotated
- Only proteins with a matched reference are included

**`hit_class_matrix.png`**
- Ligand × protein matrix, cell color = hit class (Strong=dark green, Moderate=yellow-green, Weak=orange, Inactive=light red, No data=white)
- Bottom row: **Selectivity Index** = count of proteins where ligand is Moderate or better; shown as integer with background color scaled 0→max
- Ligands sorted descending by selectivity index; proteins sorted descending by strong-binder count
- No numeric annotations in cells (class color is self-explanatory); selectivity index row shows number

---

### Group 4 — Per-Protein

**`affinity_by_protein.png`**
- Box + strip plot: one box per protein, all ligand best affinities as dots
- Best ligand for each protein annotated with gold star (⭐) and name
- Red dashed line at −7 kcal/mol; orange dotted line at 0
- Multi-engine: one color per engine, boxes side by side per protein
- Single-engine: single color

**`best_ligand_per_protein.png`**
- Horizontal bar chart, one bar per protein
- Bar = best affinity achieved at that protein (most negative)
- Bar labeled with ligand name
- Sorted ascending by best affinity (best binder at top)

**`hit_class_distribution_per_protein.png`**
- Stacked bar chart, one bar per protein
- Segments: Strong (green) / Moderate (yellow) / Weak (red) / Inactive/Uncategorized (grey)
- Y-axis: count of ligands per class
- Bars sorted descending by strong-binder count

---

### Group 5 — Engine Agreement

**Multi-engine only:**

**`engine_rank_correlation.png`**
- Square correlation matrix: engines × engines
- Cell = Spearman rank correlation of per-ligand best affinity between engine pair
- Color: blue (1.0) → white (0) → red (−1.0)
- Values annotated; diagonal = 1.00 (black)

**`engine_agreement_per_complex.png`**
- Bar chart: one bar per ligand
- Bar height = fraction of engines classifying as Moderate or better
- Bars colored: green ≥ 0.67, yellow ≥ 0.33, red < 0.33
- Sorted descending by agreement fraction

**`per_protein_cross_engine/` batch**
- Existing cross-engine per-protein batch plots retained as-is

**Single-engine only (replaces engine agreement group):**

**`pose_diversity.png`** *(Vina solo only)*
- Histogram of `rmsd_lb` values from best-pose table
- Annotates collapse warning if all rmsd_lb == 0.0

**`cnn_confidence_distribution.png`** *(GNINA solo only)*
- Stacked bar chart: count of ligands per cnn_confidence level (high/moderate/low)
- One bar per protein target

---

### Group 6 — Validation *(only when `validation_gate.json` status = validated)*

**`rmsd_distribution.png`**
- Histogram of RMSD values from `redocking_validation.csv`
- Colored: green = RMSD < 2 Å (successful redock), orange = 2–4 Å, red = > 4 Å
- Vertical line at 2 Å threshold; annotation: pass rate percentage

**`reference_vs_docked.png`**
- Scatter plot: x = reference ligand affinity, y = best docked affinity for same ligand
- One point per protein–reference pair
- Points labeled with protein name
- Diagonal line (y = x); region below diagonal = docking outperforms reference

**`rmsd_per_complex.png`**
- Scatter: x = best docked affinity, y = redocking RMSD
- One point per evaluated complex
- Points colored by RMSD threshold (< 2 Å = green, 2–4 = orange, > 4 = red)
- Vertical line at −7 kcal/mol threshold

---

## DAG Integration

### New node: `visualizations`

```
inputs:  classified_hits.csv, consensus_ranked.csv, engine_agreement.csv,
         normalized_scores.csv, engine_scope_config.json, validation_gate.json,
         top_pose_per_ligand_global.csv
outputs: visualizations/          (directory sentinel file: figures_manifest.json)
optional: True                    (figure failure must not block reports)
```

### Updated node: `reports`

Add `visualizations/` as an input. Report generator reads `figures_manifest.json` to embed figure paths in `START_HERE.md`.

### Dependency graph change

```
Before: classified_hits → reports
After:  classified_hits → visualizations → reports
                        ↑
        (also reads consensus_ranked, normalized_scores, validation_gate, top_pose_atlas)
```

---

## figures_manifest.json Schema

```json
{
  "generated_at": "<iso-timestamp>",
  "engine_mode": "multi_engine | single_engine",
  "engine": "<name or null>",
  "figures": [
    {
      "group": "overview | ranking | multi_target | per_protein | engine_agreement | validation",
      "name": "<figure filename>",
      "path": "<absolute path>",
      "title": "<human readable title>",
      "conditional": false,
      "generated": true
    }
  ],
  "skipped": [
    { "name": "...", "reason": "no_references | single_engine | validation_not_passed" }
  ]
}
```

---

## Engine-Aware Behavior Summary

| Figure | Multi-engine | Single/GNINA | Single/Vina | Single/Smina |
|---|---|---|---|---|
| `score_distribution` | per-engine subplots | single panel, cnn_affinity | single panel, vina_affinity | single panel, vina_affinity |
| `consensus_leaderboard` | consensus_score | primary_ranking_score | same | same |
| `full_affinity_heatmap` | consensus_score | primary_ranking_score | same | same |
| `series_affinity_heatmap` | consensus_score | primary_ranking_score | same | same |
| `top_ligand_profiles` | consensus_score | primary_ranking_score | same | same |
| `polypharmacology_delta_heatmap` | ✓ if refs exist | ✓ if refs exist | ✓ if refs exist | ✓ if refs exist |
| `hit_class_matrix` | ✓ | ✓ | ✓ | ✓ |
| `affinity_by_protein` | per-engine boxes | single-color box | same | same |
| `engine_rank_correlation` | ✓ | ✗ | ✗ | ✗ |
| `engine_agreement_per_complex` | ✓ | ✗ | ✗ | ✗ |
| `per_protein_cross_engine/` | ✓ | ✗ | ✗ | ✗ |
| `cnn_confidence_distribution` | ✗ | ✓ | ✗ | ✗ |
| `pose_diversity` | ✗ | ✗ | ✓ | ✗ |
| `rmsd_distribution` | ✓ if validated | ✓ if validated | ✓ if validated | ✓ if validated |
| `reference_vs_docked` | ✓ if validated | ✓ if validated | same | same |
| `rmsd_per_complex` | ✓ if validated | ✓ if validated | same | same |

---

## Implementation Module

New file: `post_docking_analysis/visualization_suite.py`

### Public API

```python
def generate_visualization_suite(
    *,
    classified_hits_csv: Path,
    consensus_ranked_csv: Path,
    engine_agreement_csv: Path,
    normalized_scores_csv: Path,
    engine_scope_config: dict,
    validation_gate: dict,
    top_pose_csv: Path,
    output_dir: Path,
    dpi: int = 300,
) -> dict:
    """
    Generate all figures. Returns figures_manifest dict.
    Per-figure failures are caught and logged; they set generated=False in manifest.
    Never raises — always returns manifest even if all figures fail.
    """
```

### Internal structure

```
generate_visualization_suite()
  ├── _plot_score_distribution()
  ├── _plot_consensus_leaderboard()
  ├── _plot_affinity_per_ligand()
  ├── _plot_full_affinity_heatmap()
  ├── _plot_series_affinity_heatmap()
  ├── _plot_top_ligand_profiles()
  ├── _plot_polypharmacology_delta_heatmap()     # guarded
  ├── _plot_hit_class_matrix()
  ├── _plot_affinity_by_protein()
  ├── _plot_best_ligand_per_protein()
  ├── _plot_hit_class_distribution_per_protein()
  ├── _plot_engine_rank_correlation()            # guarded: multi-engine
  ├── _plot_engine_agreement_per_complex()       # guarded: multi-engine
  ├── _plot_per_protein_cross_engine_batches()   # guarded: multi-engine
  ├── _plot_cnn_confidence_distribution()        # guarded: gnina solo
  ├── _plot_pose_diversity()                     # guarded: vina solo
  ├── _plot_rmsd_distribution()                  # guarded: validation passed
  ├── _plot_reference_vs_docked()                # guarded: validation passed
  └── _plot_rmsd_per_complex()                   # guarded: validation passed
```

Each `_plot_*` function:
- Accepts a single `PlotContext` dataclass (pre-loaded dataframes + config)
- Writes one file to the output directory
- Returns `FigureResult(name, path, generated, error)`
- Uses `matplotlib.use("Agg")` — no display

---

## Style Constants

All figures share a consistent style defined in `_SUITE_STYLE`:

```python
_SUITE_STYLE = {
    "dpi": 300,
    "figsize_default": (12, 7),
    "color_strong": "#2ecc71",      # green
    "color_moderate": "#f1c40f",    # yellow
    "color_weak": "#e67e22",        # orange
    "color_inactive": "#e74c3c",    # red
    "color_favorable": "#5b9bd5",   # blue (histogram)
    "color_unfavorable": "#e74c3c", # red (histogram)
    "threshold_strong": -7.0,       # kcal/mol
    "threshold_zero": 0.0,
    "font_annotation": 7,
    "font_axis": 9,
    "font_title": 11,
}
```

---

## Exit Criteria

- [ ] `figures_manifest.json` written on every run; contains `generated: true/false` per figure
- [ ] Multi-engine run on 3-engine project: all 14 base figures + delta heatmap present
- [ ] GNINA solo run: engine agreement group absent; `cnn_confidence_distribution.png` present
- [ ] Vina solo run: `pose_diversity.png` present; engine agreement absent
- [ ] Validation-passed run: all 3 validation figures present
- [ ] Validation-absent run: validation group absent; remaining figures unaffected
- [ ] Single ProLIF failure does not prevent any visualization figure from generating
- [ ] All figures render at 300 dpi with no font overlap on 25-ligand × 10-protein study
- [ ] `reports` node embeds figure paths from `figures_manifest.json` in `START_HERE.md`
- [ ] `generate_visualization_suite()` never raises; failed figures set `generated: false`
