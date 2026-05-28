# Feature Specification: Visualization + RMSD Overhaul

## Context

Current simplified post-docking outputs show three functional issues:

1. Protein naming labels are inconsistent across visualizations and do not
   enforce a global `Display Name (PDB)` format.
2. Affinity visualizations overweight global top rows and can hide
   per-protein ranking behavior.
3. RMSD stage computes one global matrix over mostly non-comparable structures
   and often yields sparse/NaN-heavy metrics.

Additional operational issues:
- OpenBabel receptor conversion emits repeated kekulization warnings.
- Interaction outputs are fragmented across separate root folders.
- ProLIF integration uses legacy API incompatible with installed ProLIF 2.x.
- PoseView exists in codebase but is not wired into simplified pipeline.

## Goals

1. Enforce global protein label format:
   - `ProteinDisplayName (PDBCODE)` when code exists
   - allow optional CLI prompting to override display names.
2. Improve affinity visual logic:
   - per-protein top performer chart (one best ligand per protein)
   - clustered per-protein affinity distribution with best-highlight
   - explicit marking of `vina_affinity > 0` as unfavorable in histogram.
3. Improve RMSD scope:
   - per complex (all poses from each ligand-protein combination)
   - per protein (best poses only)
   - global (best poses across proteins)
4. Reduce OpenBabel warning impact:
   - prefer direct receptor PDBQT -> PDB line conversion before pybel fallback.
5. Unify interaction outputs:
   - save PandaMap / ProLIF / LigPlot / PoseView under `interactions/`
   - keep pose visualization (`py3Dmol`) separate.
6. Make interaction stages resilient:
   - batch limits / graceful skip behavior
   - robust ProLIF 2.x code path.

## Out of Scope

- Replacing PandaMap itself with a different interaction engine.
- New external dependencies beyond current conda env.
- Full redesign of legacy `post_docking_analysis.pipeline` flow.

## Acceptance Criteria

1. A shared protein label helper is used in simplified + hierarchical visual
   outputs and labels include `(PDB)` when available.
2. `top_performers.png` reflects per-protein best performer logic.
3. `affinity_distribution.png` clearly marks unfavorable region (`> 0`).
4. RMSD outputs include three scopes:
   - `rmsd_analysis/per_complex_all_poses/`
   - `rmsd_analysis/per_protein_best_poses/`
   - `rmsd_analysis/global_best_poses/`
5. Simplified pipeline writes interaction tool outputs under:
   - `interactions/pandamap/`
   - `interactions/prolif/`
   - `interactions/ligplot/`
   - `interactions/poseview/` (when enabled)
6. ProLIF stage runs with installed ProLIF 2.x API path and fails gracefully.
7. `python -m py_compile` passes for touched modules.
