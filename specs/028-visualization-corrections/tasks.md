# Tasks: Visualization Corrections and Scientific Guardrails

- [x] T001 Add shared visualization helpers for `Protein name (PDB)` labels, `[Ref]` ligand tags, grouped ligand ordering, and crowded-axis layout.
- [x] T002 Update `Consensus Leaderboard`, `Affinity Per Ligand`, `Affinity by Protein`, and `Best Ligand Per Protein` to use explicit reference labels and reduce overlap.
- [x] T003 Update full/per-engine heatmaps to group co-crystal ligands ahead of series ligands and improve long-axis readability.
- [x] T004 Redefine ligand-level engine agreement to use full-support fraction and tighten engine rank correlation to matched complete rows.
- [x] T005 Replace collapsed/constant Vina pose-diversity histograms with diagnostic placeholders.
- [x] T006 Improve validation joins and placeholder messaging for `Reference vs Docked`, `RMSD Distribution`, and `RMSD Per Complex`.
- [x] T007 Patch legacy cross-engine figures to mark reference ligands and explain positive affinities as unfavorable.
- [x] T008 Add smoke coverage for scientific guardrails, explicit reference labels, and validation-affinity fallback behavior.
