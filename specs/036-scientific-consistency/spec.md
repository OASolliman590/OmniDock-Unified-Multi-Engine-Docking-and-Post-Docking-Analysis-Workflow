# Spec 036 — Scientific Consistency and Validity

**Status:** directed by the Scientific Lead on 2026-10-10. Goal: raise the scientific validity of the pipeline. Source review: an independent read-only review of scientific redundancy, stage order and usability (2026-10-10).

## Decisions (Scientific Lead, 2026-10-10)
- **D1:** One project-wide protonation policy. STI (1IEP) is set to the +1 N-methylpiperazinium microspecies as an `explicit_state`.
- **D2:** Redocking validation **marks** downstream outputs; it does not block them.
- **D3:** Remove legacy consensus modes that do not serve a valid purpose. Consensus v2 becomes the default.
- **D4:** One shared pose selector with a GNINA regression test.
- **D5:** Box size and seed handling follow the literature. Exhaustiveness is unified at the highest level across engines.

## Requirements

### R1 — One protonation policy (D1)
- R1a: The project manifest stores one `protonation_policy` block: receptor pH and force field, ligand policy (`ph_model` / `explicit_state` / `as_input`), ligand pH, and an optional state map. Receptor preparation, ligand preparation and the MD export all read it. Re-entry is unnecessary, and a stage given a conflicting value fails with `protonation_policy_conflict`.
- R1b: The MD export uses the **docked microspecies** by default and does not re-predict it. The ligand's formal charges and hydrogens come from the docked ligand's measured state (`REMARK SMILES` lineage plus preparation provenance). G6 records `charge_authority: docked_state` for `explicit_state` ligands, or `predicted` for `ph_model` ligands. Re-protonation at a different pH needs an explicit override, which is recorded.
- R1c (Spec 034 R2d): The docked microspecies (SMILES and net charge) is carried into `best_pose_per_tag_by_engine.csv`.

### R2 — Redocking marks downstream outputs (D2)
- Each consensus, top-pose, report and MD-export row carries `redocking_validation_status` (`validated` / `failed` / `not_evaluated` / `not_comparable`), its reason, and the reference used. The value is per target, from the redocking result for that receptor and site. Reports show unvalidated targets prominently. Nothing is blocked.

### R3 — Consensus clean-up (D3)
- R3a: Inventory every consensus mode. Remove a mode if it relies on a method forbidden by Spec 031: centroid "soft alignment", atom matching by file order or element order, or mapping-free RMSD. Also remove it if it duplicates v2 without a distinct, documented purpose. Removed modes fail with a clear message naming v2. Silent fallback to legacy is removed.
- R3b: `consensus_rank_geometry_qc_v2` is the default everywhere (CLI, interactive, DAG, execution).
- R3c: Single-engine projects get an explicit single-engine ranking (that engine's native score with v2 pose selection), not a forced legacy hybrid.
- R3d: Retire the legacy `_pairwise_rmsd` in `geometric_consensus.py`. `enhanced_rmsd_analyzer.py` outputs are labelled `diagnostic_not_spec031` (plots only).

### R4 — One pose selector (D4)
- One selector implementing the v2 rules (GNINA: highest `cnn_score`, with the documented v2 tie-breaks; Vina/Smina/AD4: lowest native affinity) is used by every table writer: the best-pose table, the single-engine path, `top_pose_selector` and the MD export input.
- Regression tests: GNINA rows where the minimum `cnn_affinity` is the worst pose must not be selected, and every writer must agree on a mixed-engine fixture.

### R5 — Docking parameters from evidence (D5)
- R5a: **Box size.** The default method is `rg_scaled_v1`: a cubic edge of 2.9 × the radius of gyration of the docked ligand's heavy atoms (Feinstein & Brylinski, J Cheminform 2015, 7:18). The box centre still comes from `binding_site_center_v1`.
  - A fixed size is allowed only as an explicit user value (`box_method: user_fixed`).
  - Every pair records the method, the Rg and the edge length.
  - A warning is raised when the box does not contain all reference-ligand heavy atoms, or when the edge exceeds 30 Å (Vina FAQ guidance).
- R5b: **Seed and replicates.** The seed is required and recorded for every run. The default is 3 independent seeded replicates per pair (seeds recorded). Replicates are pooled for pose selection, and a `pose_reproducibility` metric is reported: in-place RMSD between the replicates' top poses, plus the fraction within 2 Å. A single replicate is allowed and labelled `replicates: 1`.
- R5c: **Exhaustiveness.** Unified at the highest preset level (32) for every engine that has the parameter. Each engine's effective value is recorded.
- R5d: **Vina `energy_range`** is passed explicitly, recorded in the manifest, and reported when fewer poses than requested come back.

### R6 — Correctness and usability fixes (software)
- Single-engine `analyze comparative` writes `best_pose_per_tag_by_engine.csv`.
- Tags and output names use the receptor file stem, without a `.pdbqt` suffix.
- Preparation reports `completed` (exit 0) when the only notes are informational.
- The Excel/pairlist box centre requires a recorded method (`binding_site_center_v1` or `user_explicit`). Values without a method are refused.
- `autodock/prep_autodock.sh` (hard-coded PDB2PQR with flips) is retired or routed to the strict module.
- USAGE documents a single path for one target and one ligand, including `prep project`. The md-input map files are generated automatically from the project.

## Non-goals
No new engines, HPC, uploads or MD. No changes to the Spec 031 mapping rules, the Spec 032 gate semantics or the Spec 033 receptor policy.

## Acceptance
- Tests for R1–R6 pass, and the full suite and the smoke test pass.
- One bounded re-dock of 1IEP/STI on the new defaults:
  - STI +1 as an explicit state;
  - Rg-scaled box;
  - exhaustiveness 32;
  - 3 seeds.
- That re-dock produces a validated redocking mark, a G1–G8-complete MD row with `charge_authority: docked_state` (+1), and a pose-reproducibility value.
- The Scientific Lead accepts the result; the agent does not self-accept.

## References
- Feinstein W.P., Brylinski M. Calculating an optimal box size for ligand docking and virtual screening against experimental and predicted binding pockets. J Cheminform 7:18 (2015). doi:10.1186/s13321-015-0067-5
- AutoDock Vina FAQ (search space size vs exhaustiveness): https://autodock-vina.readthedocs.io/en/latest/faq.html
- Reproducibility of Vina on distributed systems, seed capture: https://pmc.ncbi.nlm.nih.gov/articles/PMC4801993
- Exhaustiveness vs repeated runs (PDBbind): https://sol.sbc.org.br/index.php/bsb/article/view/22861
