# Pending Scientific Decisions

**Status**: approved as recommended; retained as the decision matrix  
**Date opened**: 2026-08-21

Implementation must not infer these choices from software defaults.

| Decision | Recommended baseline | Required approval |
|---|---|---|
| Target pH | Require `--ph` on every direct request; use `7.4` only when explicitly supplied or recorded in a project policy. | Approve or replace. |
| Ligand protonation/charge | Allow Open Babel prediction only as `predicted`, require a charge map for `human_approved`, and mark unusual ionizable groups `human_review_required`. | Approve or require charge map for all rows. |
| Receptor source | Require explicit prepared PDB/mmCIF; refuse PDBQT conversion as MD-ready. | Approve and name the receptor source/cleaning policy. |
| Expected receptor class | `protein` for the first profile; nucleic-acid/mixed targets require another explicit profile. | Approve or replace. |
| Receptor protonation | Out of v1; preserve the approved prepared receptor unchanged and record that receptor protonation was not performed. | Approve or define a tool/force-field-specific policy. |
| Consumer profile | Implement only `charmm_gui_cgenff_v1` initially. | Approve or name another profile. |
| Export selection | Require exact `--tags-file` plus `--engine`; do not export all classified hits implicitly. | Approve or define an automatic selection rule. |
| Vina/Smina topology | Require preserved Meeko mapping or explicit SDF/MOL/MOL2 topology plus proven mapping; otherwise `not_comparable`. | Approve; unsafe bond re-perception is not offered. |
| Canonical combined name | Emit `system.pdb`; no duplicate `complex.pdb` in v1. | Approve or request compatibility alias. |
