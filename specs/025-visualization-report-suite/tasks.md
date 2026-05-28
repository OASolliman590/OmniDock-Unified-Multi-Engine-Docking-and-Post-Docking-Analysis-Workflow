# Tasks: Visualization Report Suite

**Input**: `/specs/025-visualization-report-suite/spec.md`
**Scope**: Canonical visualization suite as a single DAG node — 14+ figures in 6 groups, engine-aware, no-raise contract.

## Format: `[ID] [P?] Description`
- **[P]**: parallelizable

## Phase 1: Module Scaffold + Style

- [x] T001 Create `post_docking_analysis/visualization_suite.py` with `PlotContext` dataclass, `FigureResult` dataclass, `_SUITE_STYLE` constants.
- [x] T002 Implement `generate_visualization_suite()` entry point: load all input CSVs into `PlotContext`, dispatch to per-figure functions, catch per-figure exceptions, return `figures_manifest` dict.
- [x] T003 Implement `_write_figures_manifest()`: write `figures_manifest.json` to output dir.
- [x] T004 Implement engine-mode detection from `engine_scope_config.json`: expose `is_multi_engine`, `engine_name`, `primary_score_col`, `primary_score_label` on `PlotContext`.
- [x] T005 Implement reference detection: scan `classified_hits.csv` for rows where `ligand_type == "reference"` (or `tag` contains `Reference_`); set `PlotContext.has_references`.

## Phase 2: Group 1 — Score Overview

- [x] T006 Implement `_plot_score_distribution()`: histogram of best-pose affinities; blue=favorable, red=unfavorable; mean/median/zero lines; unfavorable annotation box; multi-engine = per-engine subplots.

## Phase 3: Group 2 — Ligand Ranking

- [x] T007 Implement `_plot_consensus_leaderboard()`: horizontal bar chart sorted by polypharmacology score (count of proteins where Moderate or better); `targets=N, better_ref=M` annotation; engine-aware score column.
- [x] T008 [P] Implement `_plot_affinity_per_ligand()`: sorted vertical bar, colored by hit class; threshold line at −7 kcal/mol.

## Phase 4: Group 3 — Multi-Target Coverage

- [x] T009 Implement `_plot_full_affinity_heatmap()`: protein × ALL ligands; reference diagonal pattern on left section; series section on right; diverging color scale; value annotations; vertical spacer column between sections.
- [x] T010 [P] Implement `_plot_series_affinity_heatmap()`: series ligands only; sorted by mean affinity descending; same color scale; annotations.
- [x] T011 [P] Implement `_plot_top_ligand_profiles()`: line plot top-N ligands across proteins; threshold line at −7; proteins sorted by mean affinity of top-N; markers per point.
- [x] T012 [P] Implement `_plot_polypharmacology_delta_heatmap()`: guarded by `PlotContext.has_references`; Δ = series − reference affinity per protein; red=worse, green=better; annotated values.
- [x] T013 [P] Implement `_plot_hit_class_matrix()`: ligand × protein colored by class; selectivity index row at bottom (count Moderate-or-better per ligand); ligands sorted by selectivity index desc; proteins sorted by strong-binder count desc.

## Phase 5: Group 4 — Per-Protein

- [x] T014 Implement `_plot_affinity_by_protein()`: box + strip plot per protein; best ligand annotated with star; threshold lines; multi-engine = per-engine side-by-side boxes.
- [x] T015 [P] Implement `_plot_best_ligand_per_protein()`: horizontal bar per protein; bar labeled with best ligand name; sorted ascending by best affinity.
- [x] T016 [P] Implement `_plot_hit_class_distribution_per_protein()`: stacked bar per protein; Strong/Moderate/Weak/Uncategorized segments; sorted by strong-binder count desc.

## Phase 6: Group 5 — Engine Agreement

- [x] T017 Implement `_plot_engine_rank_correlation()`: guarded `is_multi_engine`; Spearman rank correlation matrix; blue-white-red diverging; annotated values.
- [x] T018 [P] Implement `_plot_engine_agreement_per_complex()`: guarded `is_multi_engine`; bar per ligand; fraction of engines agreeing Moderate-or-better; green/yellow/red coloring.
- [x] T019 [P] Implement `_plot_per_protein_cross_engine_batches()`: guarded `is_multi_engine`; reuse existing `_plot_cross_engine_per_protein_batches` logic; write to `per_protein_cross_engine/` subdir.
- [x] T020 [P] Implement `_plot_cnn_confidence_distribution()`: guarded `engine_name == "gnina"`; stacked bar per protein; high/moderate/low confidence segments.
- [x] T021 [P] Implement `_plot_pose_diversity()`: guarded `engine_name == "vina"`; histogram of `rmsd_lb`; collapse warning annotation if all == 0.0.

## Phase 7: Group 6 — Validation

- [x] T022 Implement `_plot_rmsd_distribution()`: guarded `validation_gate.status == "validated"`; histogram of RMSD values; green/orange/red zones at 2 Å and 4 Å; pass-rate annotation.
- [x] T023 [P] Implement `_plot_reference_vs_docked()`: guarded by validated gate; scatter reference vs. best docked affinity; points labeled by protein; diagonal y=x line.
- [x] T024 [P] Implement `_plot_rmsd_per_complex()`: guarded; scatter best affinity vs. RMSD; colored by threshold; vertical line at −7 kcal/mol.

## Phase 8: DAG Integration

- [x] T025 Register `visualizations` node in `build_artifact_graph()`: inputs = classified_hits + consensus_ranked + engine_agreement + normalized_scores + engine_scope_config + validation_gate + top_pose_atlas; output = `visualizations/figures_manifest.json`; `optional=True`.
- [x] T026 Implement `_dag_compute_visualizations_node()` on `MultiEngineAnalysisPipeline`: load paths, call `generate_visualization_suite()`, return manifest summary dict with figure counts per group.
- [x] T027 Add `visualizations/` as input to `reports` node in `build_artifact_graph()`.
- [x] T028 Update `_dag_compute_reports_node()`: read `figures_manifest.json`; embed figure paths into `START_HERE.md` report section.

## Phase 9: Tests

- [x] T029 Smoke test: multi-engine run → all 14 base figures present; `figures_manifest.json` valid JSON.
- [x] T030 [P] Smoke test: GNINA solo → `cnn_confidence_distribution.png` present; engine agreement figures absent.
- [x] T031 [P] Smoke test: Vina solo → `pose_diversity.png` present; engine agreement figures absent.
- [x] T032 [P] Smoke test: no references → `polypharmacology_delta_heatmap.png` absent; `skipped` list contains it with `reason: no_references`.
- [x] T033 [P] Smoke test: validation gate not passed → validation group absent; remaining 11 figures present.
- [x] T034 Smoke test: inject one figure function that raises → `generate_visualization_suite()` does not raise; that figure shows `generated: false` in manifest; remaining figures generated normally.
- [x] T035 Smoke test: `reports` node output `START_HERE.md` contains at least one figure path from manifest.
