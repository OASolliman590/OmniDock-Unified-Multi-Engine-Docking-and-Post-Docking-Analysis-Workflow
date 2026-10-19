# Spec 039 — Conformer-Ensemble Docking

**Status:** specification directed by the Scientific Lead on 2026-10-19 ("add conformer docking where all conformers are docked"). Implementation follows the Spec 035 merge. The method needs validation before it becomes a default.

## Rationale
- AutoDock Vina and its derivatives sample torsions but keep ring geometry fixed (ring puckers, macrocycle shapes). Results therefore depend on the starting conformer for ring-flexible ligands. A single 3D start can also bias a stochastic search.
- **Validity flaw found in the current redocking (fix required).** With `explicit_state`, ligand preparation copies the input heavy-atom coordinates. When the input is the crystal ligand, docking starts from the crystal conformer. That inflates redocking success: it is cognate redocking from the answer. Benchmark-quality redocking starts from conformers generated from the ligand's chemical graph, not from crystal coordinates.

## Requirements

### R1 — Conformer generation (from the chemical graph, never from crystal coordinates)
- Generate conformers from the resolved microspecies (protonation policy and docked state). Use RDKit ETKDGv3 with a recorded random seed, then MMFF94s minimisation and RMSD pruning (heavy atoms). Number of conformers, pruning threshold and energy window are explicit, recorded parameters.
- Input coordinates are discarded for docking starts. When the input came from a crystal structure, provenance records `start_conformer_source: generated_from_graph` and `crystal_coordinates_used_for_start: false`. Crystal coordinates remain the redocking **reference** only.
- Every conformer has the identical microspecies, so the lineage SMILES, net charge and Spec 034 mapping are unchanged across conformers. This is verified per conformer.

### R2 — Dock every conformer
- `--conformer-docking ensemble` (default `single`) docks **every** retained conformer for each pair. Each conformer uses the same box, engine parameters and seeds per the Spec 036 replicate contract. Output names carry `__confNN` (and `__repNN`).
- The cost, which is conformers × replicates × engine runs, is shown at dry-run and recorded. A required `--max-conformers` caps it.

### R3 — Pooling, selection and convergence
- All poses from all conformers and replicates are pooled per (engine, tag) through the shared v2 selector. The selected pose records its conformer, replicate and seed.
- `conformer_convergence`: pairwise in-place RMSD between each conformer's top pose (same method as `pose_reproducibility`), reported as max, median and fraction within 2 Å, plus the spread of top scores across conformers.

### R4 — Box sizing with ensembles
- Default box size remains `rg_scaled_v1` on the prepared conformer. An option `rg_ensemble_percentile_v1` (edge = 2.9 × the chosen percentile of ensemble Rg, default 95th) is available but **experimental** until R5 validates it.

### R5 — Validation before any default change
Before conformer docking or ensemble box sizing becomes a default, benchmark on a public redocking set (Astex Diverse Set, or a PDBbind core subset of at least 30 complexes).

Arms:
- **A:** single generated conformer, `rg_scaled_v1` (the corrected baseline: no crystal start).
- **B:** all conformers docked, `rg_scaled_v1`.
- **C:** all conformers docked, `rg_ensemble_percentile_v1`.
- **Diagnostic only:** single crystal-conformer start, which quantifies the current bias.

Metrics:
- top-1 and top-3 success rate (in-place heavy-atom RMSD ≤ 2.0 Å);
- score-to-RMSD rank correlation;
- runtime;
- results split into ring-flexible and acyclic ligands.

Adoption rule: an arm becomes the default only if its top-1 success rate exceeds arm A with a 95% bootstrap confidence interval that excludes zero improvement, at acceptable runtime. The Scientific Lead decides.

Data access: the benchmark needs RCSB / PDBbind downloads (blocked in the current cloud environment), so it needs network access or a local run.

## Non-goals
No change to the scoring functions, the consensus rules or the Spec 034 lineage contract. No MD.

## Acceptance
- Tests cover: graph-based generation with a recorded seed; no crystal-coordinate start; identical microspecies across conformers; `__confNN` naming; pooling and selection across conformers; the convergence metric; and the cost cap.
- 1IEP/STI is re-run with arm A to quantify the effect of removing the crystal start (bounded run, with approval).
- The R5 benchmark is reported, and the Scientific Lead reviews it.
