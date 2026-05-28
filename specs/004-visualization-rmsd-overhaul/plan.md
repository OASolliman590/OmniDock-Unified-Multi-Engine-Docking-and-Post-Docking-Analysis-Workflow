# Implementation Plan: Visualization + RMSD Overhaul Closeout

## Current Status

The core implementation from this spec is already present in the current
post-docking stack:

- global protein label formatting is implemented in `protein_naming.py` and
  consumed by both simplified and hierarchical analysis paths
- the simplified CLI exposes `--prompt-protein-names` and `--enable-poseview`
- affinity overview plots now use the current semantics:
  - `affinity_distribution.png` marks the unfavorable region (`> 0`)
  - `top_performers.png` is one best ligand per protein
- RMSD output is split into:
  - `rmsd_analysis/per_complex_all_poses/`
  - `rmsd_analysis/per_protein_best_poses/`
  - `rmsd_analysis/global_best_poses/`
- interaction outputs are routed under `interactions/` for PandaMap, ProLIF,
  LigPlot, and optional PoseView
- user-facing docs are now aligned to the current CLI and visualization
  contract in:
  - `POST_DOCKING_ANALYSIS_GUIDE.md`
  - `post_docking_analysis/README.md`
  - `post_docking_analysis/SIMPLIFIED_IMPLEMENTATION_STATUS.md`

## Verification Completed

1. `python -m py_compile post_docking_analysis/protein_naming.py post_docking_analysis/simplified_pipeline.py post_docking_analysis/hierarchical_analyzer.py post_docking_analysis/publication_pandamap.py post_docking_analysis/prolif_interaction_maps.py post_docking_analysis/ligplot_integration.py post_docking_analysis/poseview_integration.py post_docking_analysis/enhanced_rmsd_analyzer.py`
2. `python -m post_docking_analysis.simplified_cli --help`
3. Code-path audit confirmed:
   - label formatting
   - per-protein top performer behavior
   - three-scope RMSD tree
   - canonical `interactions/` routing

## Remaining Gap

This spec is still open only because a fresh simplified-pipeline smoke run has
not yet been recorded against the documented output tree after the closeout
alignment.

## Remaining Validation Path

1. Run one simplified pipeline smoke command on a representative GNINA project
2. Inspect the resulting tree against the documented output contract
3. Record final dependency caveats for optional tools such as PandaMap,
   LigPlot+, and PoseView
