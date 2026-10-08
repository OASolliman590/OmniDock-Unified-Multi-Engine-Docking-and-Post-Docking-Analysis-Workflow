# Spec 033 / Spec 032 Phase 2 — bounded 1IEP/STI Vina run (2026-10-08)

This record has two parts: software evidence (Part A) and scientific evidence (Part B). Numeric docking scores and RMSD values are not in this public record. They are in the Spec 032 T021 review packet given to the Scientific Lead, and stay there until accepted.

- **Code:** `d4974cf` (branch `chore/install-agentic-bioinformatics-stack`), clean working tree.
- **Environment:** Linux, Python 3.10.22. Open Babel 3.2.1, PDB2PQR 3.4.1, Meeko 0.7.1, RDKit 2026.03.6, PLIP 3.0.1. AutoDock Vina is the conda-forge `vina` 1.2.7 build `py314ha160325_0`; the binary reports `f458505-mod`.
- **Raw outputs:** kept outside the repository. Raw manifests contain local absolute paths and are not committed.

## Inputs (public RCSB, `test/fixtures/spec032/phase2_1iep_inputs/`)

| File | SHA-256 |
|---|---|
| `1IEP.pdb` (identical to the Spec 031 fixture receptor) | `577386b4b8bc58a570facba0c8a843386842caef5ed7b60106c685d0f68d4218` |
| `STI_1IEP_A201.sdf` (crystal STI, chain A 201, bond orders) | `3d42bffb25f731026f1bd2501c259d7e6d5dda9dcc5ff9addfb430c435cae1a6` |
| Chain-A protein-only receptor input (derived, no H, 2229 heavy atoms) | `f65ba83470e32be99d2180c52cdb7c2d199a28152083ccc264ef71ac16a867a3` |

## Part A — Software evidence

| Stage | Status | Notes |
|---|---|---|
| `workflow init --layout-profile canonical` | completed | |
| `pdb prepare-protein` (strict module) | completed | PDB2PQR at pH 7.4 (AMBER, `--noopt --nodebump --keep-chain`), then Meeko. 2229 → 2230 heavy atoms. One accepted terminal O (chain A GLN 498). 472 HD, polar-H gate passed. Output `1062c9940d194ba90e778eb7269a5f75fe865bd875159ec93cfd3e9058ce6387`. |
| `pdb prepare-ligand` | completed | `openbabel_then_meeko`. 37 heavy atoms, 2 HD, TORSDOF 7. Output `366f48594df51125ea8edb9fae6d86ffe84b464ef82187896a6fc38aa72e7f07`. |
| `prep pairlist --mode cocrystal_only --default-box-size 20` | completed | 1 pair. Canonical root used; no `4-Docking/` created (Spec 033 R1 verified on a real project). |
| `prep project` | completed | Canonical materialisation. |
| `dock run --dry-run`, then the real run | completed, exit 0 (≈109 s) | Vina: exhaustiveness 8, `num_modes` 5, seed 12345, cpu 2, `--parameter-mode advanced`. Output `ec6985d6d950caf492d034bf1e6a21dd0adc44100790d8a1b1a76a3437aeaf8b`. |
| `analyze comparative` | completed | Single-engine DAG; validation gate `needs_review`. |
| Best-pose table | completed via repository writer driver | `fc124523a5676343d051cab92a4410f34e5a7b925a6f5a77e02da335abe76615`. |
| Spec 031 atom mapping (pose ↔ crystal) | `not_comparable` (`graphs_not_isomorphic`) for all poses | See Part B. |
| `analyze md-inputs` (vina, pH 7.4, `openbabel_predicted`, sdf) | manifest `failed`; row `not_comparable` | G1 completed. G3 `not_comparable` (`missing_pose_to_topology_atom_map`). G2 and G4–G8 not reached. |
| Repository hygiene (Spec 033 R3) | passed | `git status --porcelain` clean after every step; no config or logs in the repo root. |

### Software findings (not fixed)

1. The engine output name and tag carry the receptor file suffix (`1IEP_A_protein.pdbqt_site_1_STI_A_201`, output `….pdbqt.pdbqt`).
2. For a single-engine project, `analyze comparative` / `analyze stage reports` do not write `4-Working/scores/unified/best_pose_per_tag_by_engine.csv`, the input Spec 032 requires. The table had to be produced by calling the repository writer directly.
3. Vina wrote 4 of 5 requested modes. Its default `energy_range` of 3 kcal/mol excludes the fifth, and the pipeline does not record that limit.
4. The preparation steps return exit 1 / `partial` solely because of the informational note "No engines were selected".
5. `add_to_excel_workbook` fails on the dict-valued `binding_site_center_provenance` key (non-fatal).

## Part B — Scientific evidence (status only; values in the T021 packet)

- **Box:** `binding_site_center_v1`, the co-crystal STI chain-A heavy-atom centroid (0.0 Å from the crystal centroid), 20 Å cube. It is 3.79 Å from the historical fixture centre, which was a PLIP residue centroid.
- **Redocking reproduction:** a diagnostic in-place heavy-atom RMSD (element-and-connectivity isomorphism, bond orders ignored, symmetry enumerated, no superposition) places the top-ranked pose well within the conventional 2 Å redocking threshold. This diagnostic is **not** the accepted Spec 031 method. The accepted bond-order-aware mapping returned `not_comparable`, so redocking validation is formally **not established**.
- **Root cause of `not_comparable`:** the pose graph is rebuilt by Open Babel bond perception from PDBQT, and its Kekulé assignment does not match the crystal SDF. Meeko's ligand PDBQT already carries `REMARK SMILES` and `REMARK SMILES IDX`, a bond-order-bearing lineage from PDBQT atoms to a chemical graph. Using it would be a new mapping source and needs a Scientific Lead decision and its own spec. It was not implemented.
- **Ligand protonation for docking:** the Meeko SMILES shows a **neutral** N-methylpiperazine, although the step reports "protonation at pH 7.4 applied". The docked microspecies is therefore neutral imatinib, and it needs Scientific Lead review.
- **Receptor lineage:** the docking receptor is PDB2PQR-protonated. The MD-export receptor-map entry is the unprotonated chain-A PDB, per the approved Spec 032 v1 receptor boundary.
- **MD export:** no completed row exists, so G2 and G4–G8, the predicted formal charge and the charge authority are not established for a real ligand.

Docking scores are hypothesis-generating rankings, not affinity measurements.
