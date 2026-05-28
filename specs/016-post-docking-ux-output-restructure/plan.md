# Implementation Plan: Post-Docking UX And Output Restructure

## Scope

This implementation focuses on interactive usability, deterministic analysis
progress reporting, and reorganized interaction/visualization outputs without
breaking existing CLI entry points.

## Implementation Notes

1. Update workflow-interactive prompt wrappers to support back/exit semantics
   and controlled cancellation handling.
2. Add a standalone naming-management panel in interactive workflow that writes
   shared protein/ligand override files under post-docking analysis metadata.
3. Route one interactive multi-target analysis launch into a timestamped
   `analysis/sessions/<session_id>` root and reuse that root across selected
   stage targets.
4. Extend analysis dispatch to pass a shared naming directory so stage-only runs
   can reuse previously curated name mappings.
5. Refactor simplified-pipeline run orchestration to emit step progress markers
   and fix RMSD scope-signature generation.
6. Improve ProLIF output organization (`across_poses`, `best_pose_2d`,
   summaries), and auto-enable PoseView for poseview-specific targets.
7. Harden PyMOL visualization execution status handling and enrich interaction
   script generation to emphasize pocket-focused residue labeling.
8. Normalize PandaMap naming/output conventions (`maps_2d`, `maps_3d`) and
   remove duplicated protein-prefix naming artifacts.
9. Promote a canonical visualization tree (`visualizations/2d`, `visualizations/3d`)
   by mirroring renderable assets from tool-specific raw outputs.
10. Generate consolidation indexes (`visualization_manifest`, `raw_artifact_inventory`)
    plus a root `START_HERE.md`, then prune empty directories.
