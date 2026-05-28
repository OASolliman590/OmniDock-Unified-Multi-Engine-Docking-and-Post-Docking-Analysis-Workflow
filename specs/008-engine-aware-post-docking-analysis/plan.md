# Implementation Plan: Engine-Aware Post-Docking Analysis

1. Add `post_docking_analysis/multi_engine_pipeline.py` to:
   - discover canonical multi-engine projects
   - load or build normalized scores
   - emit compatibility exports that preserve a unified score contract
   - emit comparative reports
   - filter to one engine or a favorite engine
   - continue GNINA favorites through the simplified pipeline
   - continue Vina/Smina favorites through a structural downstream bridge
2. Extend `post_docking_analysis/cli.py` with engine-aware arguments.
3. Keep the legacy CLI flow unchanged when a canonical project is not used.
4. Reuse `SimplifiedPostDockingPipeline` as the GNINA downstream bridge and
   `PostDockingAnalysisPipeline` report/quality helpers for Vina/Smina
   continuation.
5. For Vina/Smina, preload simplified-pipeline-compatible score tables and
   best-pose complex PDBs so the safe subset of the simplified stack can run
   without GNINA/SDF-only preprocessing, then run a PDB-based RMSD bridge from
   the extracted complexes.
6. Emit explicit bridge stage-status notes so optional stages are reported as
   `completed`, `disabled`, `missing_dependency`, `missing_configuration`, or
   partial external-render outcomes instead of generic booleans.
7. Validate with:
   - `python -m py_compile`
   - `python -m post_docking_analysis --help`
   - engine-aware dry/smoke invocations on a prepared project, including
     favorite-engine continuation for at least one PDBQT-based engine and
     verification of generated `rmsd_analysis/` outputs
