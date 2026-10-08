# Spec 033 — Preparation and Layout Integrity

**Status:** specification directed by the Scientific Lead on 2026-10-08 ("radical", root-cause fixes, option A for receptor protonation); implementation in progress; human scientific review pending.
**Trigger:** the bounded Spec 032 Phase 2 check (1IEP/STI, Vina) stopped before docking on two defects. Neither defect is specific to that run.

## Problem

### P1 — Layout profile is hard-coded inside workflow steps
`workflow/execution.py` (`run_prepare_pairlist`, plus a sibling step) and several call sites in `workflow/interactive.py` call `ensure_project_layout(root, "docking_legacy")` regardless of the project's recorded layout profile. In a canonical (8-stage) project, `prep pairlist` therefore creates `4-Docking/` and writes `pairlist.csv` at the project root. `dock run` then detects the legacy profile from the presence of `4-Docking/` and reads `4-Docking/pairlist.csv`, which does not exist. As a result, the canonical topology required by AGENTS.md cannot be driven end to end.

### P2 — Receptor preparation silently produces hydrogen-free receptors
`prep_autodock_enhanced.sh` falls back to `obabel -xr` when Meeko is unusable. That path adds no hydrogens, and PDB2PQR is off by default. Neither Meeko nor PDB2PQR is declared in this branch's dependency manifests. A clean install from the declared manifests therefore yields receptor PDBQT files with zero polar hydrogens (`HD`) and reports the step as successful. This removes the explicit hydrogen-bond donor information that Vina-family scoring needs, and it does so without warning.

### P3 — Preparation writes configuration into the caller's working directory
`autodock_preparation.py` (`create_config_file` defaulting to `autodock_config.json`) and the shell script write configuration and log files containing absolute local paths into the current working directory. That can be a public repository root.

### P4 — Branch divergence (record only)
`main` contains commit `7fc130a` ("Correct scientific preparation, docking execution, and post-docking validation"), which this branch never received. A trial merge conflicts in 15 files, including consensus, redocking RMSD, pose selection and artifact graph modules that overlap with approved Spec 031 methods. **Reconciliation is out of scope for Spec 033** and is a Scientific Lead decision (see `decisions/main-divergence.md`).

## Requirements

- **R1:** Every workflow step resolves the layout profile from the project manifest through one resolver. No step may create a layout other than the project's recorded profile. Legacy projects behave exactly as before. A canonical project can run `workflow init → prep pairlist → dock run --dry-run` and finds its own pairlist.
- **R2:** Receptor preparation is strict and fails closed. It never silently substitutes a backend and never silently emits a receptor without polar hydrogens.
  - R2a: Port the strict `docking/preparation/receptor_preparation.py` (with `structure_contract.py` and the required `receptor_quality.py` changes) from `main@7fc130a`. Keep its behavior: refuse multi-model or altloc-ambiguous input, check heavy-atom/coordinate conservation, write provenance JSON (backend, tool version, commands, requested pH, protonation status, hashes) and a preparation log next to the output.
  - R2b: Route `prep_autodock_enhanced.sh` receptor preparation through that module. Remove the silent `obabel -xr` fallback. When a required backend is missing, the receptor is `failed` with reason `skipped_missing_dependency`/`backend_unavailable`, never `completed`.
  - R2c: **Option A (Scientific Lead direction 2026-10-08):** for docking receptors, protonate with PDB2PQR at the explicit configured pH (default 7.4) and force field (default AMBER) before Meeko conversion. `receptor_use_pdb2pqr` defaults to `true`. Disabling it requires an explicit configuration value and is recorded as `protonation_status: template_selected_not_ph_titrated`.
  - R2c decision (option 1, 2026-10-08): PDB2PQR runs with `--noopt --nodebump` (no side-chain flips). The conservation gate accepts only a C-terminal carboxylate O added on the last residue of a chain, at most one per terminus, and records it in provenance. See `decisions/pdb2pqr-heavy-atom-policy.md`.
  - R2d: Polar-hydrogen gate. A receptor PDBQT for a Vina-family engine (vina, smina, qvina, autodock4) must contain at least one `HD` atom, and protein N–H donors must be represented. Otherwise preparation fails with an explicit reason.
  - R2e: Declare `pdb2pqr` and `meeko` in `environment.yml` and the applicable requirements files, with a capability check that reports their versions.
- **R3:** Preparation configuration and logs go into the project (`.meta/` or the project's log directory), never the process working directory. Persisted configuration contains no absolute paths outside the project, other than explicit user-supplied inputs.
- **R4:** Record P4 for Scientific Lead decision. No merge.

## Non-goals

- No change to Spec 031 or Spec 032 approved methods, scoring, consensus, RMSD, reference policy, pose selection or MD-export gates.
- No ligand-preparation method change. Ligand protonation stays as currently implemented and is recorded as a limitation.
- No receptor protonation inside the Spec 032 MD export (unchanged v1 boundary).
- No new docking engines, HPC, uploads or MD.
- No merge or rebase with `main`.

## Known scientific limitations (to disclose, not solve)

- PDB2PQR removes HETATM records (ligands, ions, cofactors, waters). Receptors that need retained cofactors or metals require a separate explicit decision.
- PDB2PQR/PROPKA pKa assignment is a prediction; histidine tautomers and titratable states are not experimentally confirmed.

## Acceptance

1. New focused tests pass for R1–R3, and the full suite plus the smoke test pass. Exact counts are recorded.
2. A canonical project completes `prep pairlist` then `dock run --dry-run` without a legacy folder being created.
3. Preparing the 1IEP chain-A receptor yields a PDBQT with `HD` atoms and a provenance JSON stating `pdb2pqr_applied` at pH 7.4.
4. Removing PDB2PQR or Meeko from `PATH` makes receptor preparation fail with an explicit reason, not succeed.
5. The Scientific Lead reviews and accepts (T-final); the agent does not self-accept.
