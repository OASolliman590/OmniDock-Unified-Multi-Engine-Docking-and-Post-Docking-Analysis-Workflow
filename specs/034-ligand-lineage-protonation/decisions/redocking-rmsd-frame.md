# Decision needed — redocking RMSD frame (in place vs Kabsch-superposed)

**Status:** decided 2026-10-10 by the Scientific Lead: option 1. Redocking validation reports in-place heavy-atom RMSD as the primary value, with Kabsch-superposed RMSD as secondary. Inter-engine geometric consensus keeps Kabsch. Both values use the same Spec 031 symmetry-aware mapping. **Implemented (T002a)** in `post_docking_analysis/atom_mapping.py` (`compare_graph_poses_in_place`, method id `in_place_v1`; Kabsch secondary as `kabsch_secondary`) and `post_docking_analysis/redocking_validation.py` (primary `redocking_rmsd_angstrom` is in-place; thresholds unchanged at 2.0 / 3.5 Å; new columns `redocking_rmsd_frame`, `kabsch_secondary_rmsd_angstrom`, `rmsd_secondary_method`). Tests: `test/test_spec034_lineage.py`. Geometric consensus unchanged (Kabsch).

The accepted Spec 031 method (`spec031-atom-mapping-v1`) computes RMSD after Kabsch superposition. Redocking validation conventionally reports **in-place** heavy-atom RMSD: the docked pose and the crystal ligand already share the receptor frame, and superposition can hide a misplaced or rotated pose.

On the Phase 2 1IEP/STI poses, a diagnostic in-place RMSD placed the top pose within the conventional 2 Å threshold and the other poses more than 12 Å away. Superposed values would be smaller and would not separate correctly from misplaced poses.

## Options
1. **(Agent recommendation)** For redocking validation only, report in-place RMSD as the primary value and the Kabsch value as secondary. Inter-engine geometric consensus keeps Kabsch, because there the frames are also shared but pose-shape agreement is the question being asked.
2. Keep Kabsch everywhere (status quo) and document the limitation.
