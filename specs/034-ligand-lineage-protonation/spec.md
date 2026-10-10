# Spec 034 — Ligand Lineage and Protonation Integrity

**Status:** specification directed by the Scientific Lead on 2026-10-10 (bond-order lineage mapping agreed; protonation and pH entry "designed to fit any case"). Human scientific review pending.
**Trigger:** the Spec 033 Phase 2 run (`specs/033-prep-layout-integrity/evidence/phase2_1iep_run_20261008.md`).

## Problems

- **L1, missing bond-order lineage after docking.** Spec 031 graph mapping and Spec 032 G3 need a bond-order-bearing graph for each docked pose. Vina/Smina poses are PDBQT, so the pose graph is currently rebuilt with Open Babel bond perception. On 1IEP/STI that perception produced a Kekulé assignment different from the crystal SDF, so every pose was `not_comparable` and the MD export stopped at G3. The Meeko ligand PDBQT already carries the lineage needed (`REMARK SMILES` and `REMARK SMILES IDX`, which maps PDBQT atom serials to SMILES atom indices), and Vina copies those remarks into every output model.
- **L2, silent protonation no-op.** Ligand normalisation runs `obabel --gen3d` (which adds explicit hydrogens) and then `obabel -h -p <pH>`. Open Babel does not re-protonate a molecule that already has explicit hydrogens, so the pH step changes nothing, yet provenance records `protonation_applied: true`. Verified on STI: the docked microspecies is neutral. Applied correctly to the hydrogen-free molecule, Open Babel's pH model gives the +2 dication. Imatinib's literature microspecies at pH 7.4 is predominantly +1 (distal N-methylpiperazine), so a single predictor is not safe for every case.
- **L3, receptor pH is a silent default.** PDB2PQR's pH comes from configuration with a default of 7.4. The Scientific Lead requires that the pipeline let the user choose the pH and enter it manually.

## Requirements

### R1 — Bond-order lineage mapping (L1)
- R1a: Read `REMARK SMILES` and `REMARK SMILES IDX` from Meeko-prepared ligand PDBQT and from each docked pose model. Build the pose graph from the SMILES (bond orders, formal charges) with heavy-atom coordinates taken from the pose via the IDX map. Hydrogens follow the same rule as the existing Spec 031 heavy-atom contract.
- R1b: Register this as a new lineage source in the Spec 031 mapping module (`post_docking_analysis/atom_mapping.py`) with method id `meeko_smiles_idx_lineage_v1`. The existing bond-labelled isomorphism, symmetry enumeration and complete-coverage rules are unchanged, and Open Babel perception remains a labelled fallback. Record the lineage source in every mapping result.
- R1c: Spec 032 (Vina/Smina) can derive the pose-to-topology `atom_map` from this mapping when the topology map leaves `atom_map` empty and the pose carries lineage remarks. The derived map is recorded as `atom_map_source: meeko_smiles_idx_lineage_v1`. An explicit `atom_map` still wins, and a missing or ambiguous lineage stays `not_comparable`.
- R1d: Redocking validation reuses this mapping for Vina/Smina poses.

### R2 — Ligand protonation policy for any case (L2)
- R2a: Fix the no-op. When a pH model is requested, remove existing hydrogens before applying `-p <pH>` (or apply it before 3D generation), then add explicit hydrogens.
- R2b: Provenance reports what was **measured** on the output, not that a command ran: net formal charge, charged-atom list and an output microspecies SMILES. `protonation_applied` is true only when that measurement exists.
- R2c: A per-ligand protonation policy, chosen explicitly:
  - `ph_model` (Open Babel pH model, flagged `human_review_required`);
  - `explicit_state` (a user-supplied microspecies SMILES or an approved per-ligand net charge; preparation verifies the prepared ligand matches it and fails otherwise);
  - `as_input` (keep the input's hydrogens and charges, recorded as such).
  The default is `ph_model` with review flagged. The CLI accepts a per-ligand state map file. Interactive mode shows the measured charge and asks for confirmation or an explicit state.
- R2d: The docking-input microspecies (SMILES plus net charge) is carried into `best_pose_per_tag_by_engine.csv` provenance so that Spec 032 G6 can compare it with the MD-export state.

### R3 — Explicit receptor pH entry (L3)
- R3a: `pdb prepare-protein` and the interactive pipeline require the PDB2PQR pH as an explicit user value. No silent default is applied: the CLI fails without `--ph`, and interactive mode prompts for it, showing 7.4 as a suggestion that must be confirmed. The PDB2PQR force field is selectable the same way.
- R3b: Provenance records `ph_source: user_entered` (or `config_file` when it comes from an explicit project configuration value), plus the value.

## Non-goals
- No change to Spec 031 scoring, consensus or reference policy, Redocking RMSD changes from Kabsch-only to in-place primary with Kabsch secondary, per the decision in `decisions/redocking-rmsd-frame.md`. Geometric consensus is unchanged.
- No new docking engines, uploads or MD.
- No automatic pKa prediction beyond Open Babel's existing model. Better predictors are a separate decision.

## Acceptance
1. On the retained Phase 2 poses (1IEP/STI), the lineage mapping returns a complete mapping for every pose, and the Spec 032 row passes G3 and proceeds through G4–G8, or fails for an explicit, scientifically meaningful reason.
2. A focused test shows that the pre-fix `-h -p` sequence leaves STI neutral and the fixed path changes the measured charge, while `protonation_applied` is false whenever no measurement exists.
3. An `explicit_state` test (STI as the +1 N-methylpiperazinium SMILES) prepares a matching ligand, and a mismatching state fails.
4. `pdb prepare-protein` without `--ph` fails with a clear message, and interactive mode prompts for the pH.
5. Full suite plus the smoke test pass, with exact counts recorded. Scientific Lead acceptance is required; the agent does not self-accept.
