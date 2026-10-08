# Spec 032 Scientific Baseline Approval

**Decision**: approved  
**Scientific Lead approval date**: 2026-08-21  

The Scientific Lead approved the Spec 032 recommended baseline, including:

- explicit target pH on every direct invocation;
- Open Babel protonation and formal charge recorded as `predicted` and `human_review_required`;
- protein-only `charmm_gui_cgenff_v1` consumer profile;
- unchanged, explicitly mapped prepared receptor PDB/mmCIF input;
- exact tag selection plus engine selection;
- strict ligand topology lineage, with no bond-order inference from arbitrary PDBQT;
- receptor protonation outside v1;
- `system.pdb` as the canonical combined output filename.

This approval authorizes bounded local implementation and verification. It does not authorize CHARMM-GUI uploads, docking or MD execution, HPC submission, cloud resources, or other remote mutations.
