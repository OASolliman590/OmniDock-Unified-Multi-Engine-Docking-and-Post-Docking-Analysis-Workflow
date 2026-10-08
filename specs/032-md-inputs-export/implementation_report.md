# Spec 032 Implementation Report

**Date:** 2026-08-21  
**Engineering status:** completed  
**Scientific review status:** pending  
**Named-consumer status:** `reported_unverified`

## Implemented baseline

- The run requires an explicit pH and records it in configuration and row provenance.
- `openbabel_predicted` is the only v1 protonation policy. Its formal charge is labeled predicted and `human_review_required` unless an approved per-tag charge map is supplied.
- `charmm_gui_cgenff_v1` is restricted to protein receptors. DNA-only or otherwise non-protein receptors fail closed without row outputs.
- Receptor coordinates come from the explicitly mapped prepared PDB or mmCIF lineage. Receptor PDBQT is never accepted as the MD receptor source, and this stage does not alter or protonate the receptor.
- Selection is exact by engine and requested tag from `best_pose_per_tag_by_engine.csv`; the selected one-based pose index and native score are checked.
- GNINA uses the selected SDF record as geometry and topology. Vina/Smina require an explicit connectivity-bearing topology plus a complete atom permutation; PDBQT-only rows are `not_comparable`.
- Each successful row emits `receptor.pdb`, `ligand.mol2`, optional requested ligand fan-out, `system.pdb`, and `provenance.json` under `5-Analysis/md_inputs/{protein_key}/{ligand_key}/{engine}/`.
- G1-G8, source/output SHA-256 values, backend identity/version, charge authority, pH, docking lineage, and force-field status are recorded. Publication is atomic and cache reuse is content-hash guarded.
- `md_inputs` is an explicit optional DAG scope and is not part of ordinary `full` analysis or report generation.

## Verification evidence

| Verification | Result |
|---|---:|
| Spec 032 focused suite with OpenBabel 3.2.1 | 16 passed, 0 skipped |
| Combined Spec 031 + Spec 032 regression with OpenBabel installed | 70 passed |
| Web UI regression | 107 passed, 1 skipped |
| Minimal repository bioinformatics smoke (`--skip-all-engines --skip-prep-matrix`) | passed, exit 0 |
| Python `compileall` | passed |
| `git diff --check` | passed; line-ending warnings only |
| JSON schema parse and manifest/hash assertions | passed in focused suite |
| DAG cache reuse and content invalidation | passed in focused suite |

OpenBabel 3.2.1 was installed in the repository virtual environment from the official CPython 3.12 Windows wheel and pinned in the dependency manifests. The dependency-backed test executed `CorrectForPH(7.4)`, explicit-hydrogen addition, zero-implicit-H validation, formal-charge capture, heavy-coordinate preservation, unique naming, MOL2 graph parse-back, strict `system.pdb` assembly, and end-to-end provenance publication. Predicted states remain `human_review_required`; this engineering verification does not scientifically approve OpenBabel's predicted state for any real ligand.

## Remaining evidence gates

1. **Named-consumer execution:** T020 is satisfied by explicitly recording `reported_unverified`, but no upload or external CHARMM-GUI/CGenFF run was authorized or performed. An OpenBabel round-trip is not substituted for that evidence.
2. **T021 scientific acceptance:** no retained real 2OV5/4ZBE/6LL5 project, authoritative `best_pose_per_tag_by_engine.csv`, prepared receptor map, or corresponding pose files are present in this repository. Those inputs are required to generate the per-row pH, pose/score, predicted or approved charge, G1-G8, and hash bundle for Scientific Lead acceptance.

No force-field parameterization, receptor protonation, solvation, MD system building, remote upload, HPC execution, or biological interpretation is included.
