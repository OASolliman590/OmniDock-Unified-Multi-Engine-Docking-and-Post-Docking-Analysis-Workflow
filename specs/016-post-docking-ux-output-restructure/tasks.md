# Tasks: Post-Docking UX And Output Restructure

- [x] Add interactive prompt cancellation/back/exit controls in `workflow/interactive.py`.
- [x] Add standalone interactive naming editor stage for protein/ligand display aliases.
- [x] Sessionize interactive multi-target analysis output under `analysis/sessions/<timestamp>/`.
- [x] Add shared naming-directory plumbing in workflow analysis execution paths.
- [x] Fix RMSD scope-signature NameError and add step-progress logging in simplified pipeline.
- [x] Reorganize ProLIF output to include across-pose and best-pose maps with summaries.
- [x] Auto-enable PoseView when running `analyze.interactions.poseview`.
- [x] Harden PyMOL stage success/failure signaling and apply pocket-style script tuning.
- [x] Route PandaMap publication outputs to `maps_2d` / `maps_3d` and keep legacy-path compatibility in consolidation.
- [x] Fix PandaMap filename aliasing to avoid duplicated protein prefixes in exported map names.
- [x] Mirror visual assets into canonical `visualizations/2d` and `visualizations/3d` tool subtrees.
- [x] Add output navigation/index artifacts: `START_HERE.md`, `visualization_manifest`, and `raw_artifact_inventory`.
- [x] Prune empty directories during output consolidation.
- [x] Run targeted verification commands and record closeout notes.
