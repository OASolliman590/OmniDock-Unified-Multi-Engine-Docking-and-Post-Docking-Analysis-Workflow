# Implementation Plan: Simplified Pipeline Deduplication

1. Add helper methods in `SimplifiedPostDockingPipeline`:
   - `_best_poses_from_scores()`
   - `_get_complex_pdb_files()`
2. Replace duplicated inline logic in:
   - `_analyze_rmsd`
   - `_generate_reports`
   - `_extract_poses`
   - `_generate_py3dmol_visualizations`
   - `_generate_prolif_interaction_maps`
   - `_generate_ligplot_diagrams`
   - `_generate_pandamap_analysis`
3. Validate with:
   - `python -m py_compile ...`
   - simplified CLI smoke run (`--no-rmsd --no-visualizations`) in `pdb-prepare-wizard`.
