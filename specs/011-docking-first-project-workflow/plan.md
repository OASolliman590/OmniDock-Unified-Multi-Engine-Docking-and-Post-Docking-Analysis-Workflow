# Implementation Plan

1. Add a staged `docking_legacy` layout profile to the project layout helpers.
2. Route workflow init and the interactive shell through that layout profile by
   default.
3. Add project-aware PDB collection that writes source/cleaned proteins and
   extracted ligands into the staged raw folders.
4. Add a dedicated pairlist builder with cocrystal-first logic and curated
   selection modes.
5. Materialize `4-Docking/` from prepared assets and the generated pairlist.
6. Update docking and post-docking path resolution helpers so the staged project
   root works across GNINA, Vina, and Smina.
7. Update dependency manifests and docs for `questionary`.
