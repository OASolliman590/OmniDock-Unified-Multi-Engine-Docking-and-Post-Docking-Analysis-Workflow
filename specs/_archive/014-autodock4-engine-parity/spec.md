# Feature Spec: AutoDock4 Engine Parity and Unified Ligand Prep

## Summary
Introduce AutoDock4 (`autogrid4` + `autodock4`) as a first-class docking engine in the unified workflow, with pairlist-first planning and consistent interactive/CLI/HPC deploy surfaces.

## Requirements
1. Engine selection must include `autodock4` anywhere users choose docking panels/engines.
2. Docking run/deploy flows must accept AutoDock4 runtime options and produce executable manifests.
3. AD4 runner must generate per-pair `grid.gpf` and `docking.dpf`, run `autogrid4` then `autodock4`, and write `.dlg` outputs.
4. AD4 scores must be normalized into the same schema used by current multi-engine analysis pipelines.
5. Project materialization must support hybrid ligand layout: shared canonical ligands plus an engine-specific AutoDock4 ligand mirror.
6. Built-in HPC profile templates should include baseline `autodock4` settings.

## Constraints
- Keep pairlist-first contract.
- Preserve existing GNINA/Vina/Smina behavior.
- Avoid introducing a separate legacy AD4-only folder-scan path in this iteration.
