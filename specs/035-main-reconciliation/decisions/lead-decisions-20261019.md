# Spec 035 — Scientific Lead decisions (2026-10-19)

The Scientific Lead accepted every agent recommendation in `inventory.md`, section "Decisions":

- **D1** Port order: operational ports first, in small tested commits. **Yes.**
- **D2** AutoDock4 grid points: `main`'s even ceiling (54 for 20 Å at 0.375 Å). Check the AutoGrid convention in the port test.
- **D3** Ligand pH step failure: **refuse** (`main`). The labelled RDKit fallback is removed.
- **D4** Refuse PDB-only ligands; do not re-minimise ligands that already have 3D coordinates. **Yes.**
- **D5** Remove the raw PDB to SDF conversion with inferred bonds in `workflow/execution.py`. **Yes.**
- **D6** Remove the atom-type guessing in `ligand_quality.py`. **Yes.**
- **D7** Pocket volume and druggability: report "Not evaluated" (`main`).
- **D8** Redocking gate: `main`'s stricter rule (coverage 1.0, at least 3 unique experimental complexes, pass rate 0.80).
- **D9 / D10** No raw cross-engine comparisons (consistent with Lead decision 6 and Spec 037 R4b).
- **D11** Inter-engine geometric agreement stays Kabsch. **Yes.**
- **D12** Legacy consensus modes stay removed. **Yes.**
- **D13** Do not port `score_import.py` now.
- **D14** Do not port RCSB reference retrieval now (network policy).
- **D15** QC group (Vina missing RMSD becomes NaN; no false collapse warning; keep record name and chain in exports; unique biology keys). **Yes to all.**
- **D16** Legacy normalisers: missing values become NaN. **Yes.**
- **D17** Clear a reference whose frame differs from the receptor. **Yes.**
- **D18** Basic mode refuses explicit `--exhaustiveness` (Lead decision 5).
- **D19** Do not import `audit/` files. Import only sanitised findings.

**Privacy (2026-10-19):** the audit logs that exposed a local path and host name were removed from `main` in `ec327f4`. Rewriting history to purge them from earlier commits was **not** authorised.

## Follow-up decisions (2026-10-19, after the merge)
1. **DAG write ownership (O6/S9): the analysis DAG is the single owner** of the shared report outputs (for example the top_pose_atlas files). Other report nodes and the legacy `run()` path read those outputs and do not rewrite them. With one owner, port main's output-content fingerprint and its "inputs changed during computation" failure.
2. **Port S18:** the pair reference fields (`reference_frame_id`, `receptor_frame_id`, and `asset_identity`), so that the D17 cross-frame reference guard takes effect.
3. **D16 degenerate groups:** a legacy-normaliser group with one finite value or zero span gives **NaN** (undefined), not 0.0 or 0.5.
