# Spec 028: Visualization Corrections and Scientific Guardrails

## Problem

The current post-docking visualization suite has three classes of defects:

1. Readability defects: crowded categorical axes, overlapping best-hit annotations, and protein labels that omit the requested `Protein name (PDB)` context.
2. Reference/context defects: co-crystal ligands are not consistently labeled or grouped across ligand-centric figures and heatmaps.
3. Scientific-defensibility defects: several figures can render misleading outputs when agreement, pose-diversity, or validation inputs are weak, collapsed, or partially missing.

## Goal

Improve the visualization suite so generated figures are readable, explicitly identify co-crystal ligands, and degrade safely when the underlying metric is not defensible.

## Scope

- `post_docking_analysis/visualization_suite.py`
- `post_docking_analysis/_legacy_multi_engine_pipeline_impl.py`
- `test/test_dockforge_smoke.py`

## Required Changes

### Shared labeling/layout
- Use `Protein name (PDB)` labels in visualization outputs when a protein mapping is available.
- Add a shared ligand label helper that appends `[Ref]` to co-crystal/reference ligands.
- Increase figure size dynamically for crowded categorical plots and wrap long tick labels.

### Figure-specific corrections
- Reduce annotation overlap in `Affinity by Protein`.
- Mark reference ligands in `Affinity Per Ligand`, `Best Ligand Per Protein`, and `Consensus Leaderboard`.
- Group reference ligands ahead of series ligands in full/per-engine/cross-engine heatmaps.
- Rework `Engine Agreement Per Ligand` to report a ligand-level full-support fraction instead of a visually misleading average.
- Restrict `Engine Rank Correlation` to matched complexes with complete multi-engine coverage.
- Treat all-zero/constant Vina `rmsd_lb` as a diagnostic placeholder instead of a normal diversity histogram.
- Improve validation figures so missing baselines / missing affinity joins emit explicit placeholders rather than empty-looking plots.
- Clarify that positive affinities in cross-engine matched-pair distributions are unfavorable, not silently normal.

## Constraints

- Keep existing figure filenames stable.
- No new DAG nodes or new artifact types.
- Figure generation must remain no-raise and continue writing the manifest even when a figure degrades to a placeholder.

## Non-Goals

- Rewriting upstream scoring or consensus math beyond the minimum needed to make plotted metrics defensible.
- Pixel-perfect styling redesign of the whole suite.
