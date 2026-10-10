# Spec 037 — Multi-Entity Coherence

**Status:** directed by the Scientific Lead on 2026-10-10 (see `specs/036-scientific-consistency/decisions/lead-decisions-20261010.md`, item 8). The goal is one coherent pipeline for projects with many proteins and many ligands.

## R1 — Layered protonation policy
- R1a: The policy has three levels. A **project default** (receptor pH and force field, ligand policy and ligand pH) sits at the top. **Per-receptor overrides** (pH, force field) and **per-ligand overrides** (policy, pH, explicit state) sit below it. Overrides are stored in the manifest, for example `protonation_policy.receptors.<name>` and `protonation_policy.ligands.<name>`, and set through `workflow protonation-policy --receptor <name> ...` / `--ligand <name> ...`. A ligand state map can carry many ligands.
- R1b: **Prepare once, reuse everywhere.** Each receptor and each ligand is prepared once under its resolved policy and reused in every pair. The same ligand docked into different proteins is the same microspecies, and the same protein screened against many ligands is the same protonated receptor. Preparation provenance records the resolved policy and the level it came from (project, receptor or ligand).
- R1c: **Pair coherence check.** At pairlist and dock time, each pair records the receptor and ligand resolved states. If a `ph_model` ligand was prepared at a pH different from its receptor's pH, the pair is flagged `pair_ph_mismatch`. `explicit_state` ligands carry no pH and are not flagged. Flags are marks, not blocks, and they flow to the best-pose table and the reports.
- R1d: Re-preparation is required when an entity's resolved policy changes (content-hash guard). Stale prepared files are refused.

## R2 — Ligand radius of gyration computed at preparation
- R2a: Ligand preparation computes the heavy-atom radius of gyration of the prepared conformer once (`heavy_atom_unweighted_v1`) and stores it, with its method, in the ligand provenance.
- R2b: The box step (`rg_scaled_v1`) reads Rg from that provenance. When it is missing it computes it, says so, and records `rg_source: computed_at_box_step`. The box records `rg_source` (`ligand_preparation` or `computed_at_box_step`).
- R2c: The Rg definition (unweighted heavy atoms on the prepared conformer) is still an unverified assumption, because the methods section of Feinstein & Brylinski 2015 could not be retrieved. It is recorded as such.

## R3 — Remaining Lead decisions (2026-10-10)
- R3a: Basic mode refuses an explicit `--exhaustiveness` (as it refuses `--num-modes`), and the effective value and source are recorded.
- R3b: `_best_by_tag`, "top overall" and best-per-protein summaries use the v2 consensus rank (or `single_engine_native_v1`). No raw cross-engine affinity comparison remains.
- R3c: The GNINA `summary.txt` label names the actual selection key (`cnn_score`).

## Acceptance
- A synthetic project with 2 receptors × 2 ligands (one `explicit_state` ligand, one `ph_model` ligand, one receptor override) prepares each entity once, reuses them across pairs, and records the policy levels.
- A pH mismatch is flagged on the right pairs only.
- A changed policy forces re-preparation.
- Rg is in the ligand provenance and the box reads it.
- R3 items are tested.
- Full suite and both smoke tests pass.
- Scientific Lead review.
