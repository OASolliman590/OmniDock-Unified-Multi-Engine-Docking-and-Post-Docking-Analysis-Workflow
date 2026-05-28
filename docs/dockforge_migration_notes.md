# Omni-DockForge Migration Notes

This guide maps legacy surfaces to the Omni-DockForge workflow model and clarifies compatibility behavior during the transition.

## Naming migration

- Old platform label: `PDB Prepare Wizard`
- New platform label: `Omni-DockForge`
- New tagline: `End-to-End Docking, Consensus by Design.`
- Old workflow concept: `maturation`
- New workflow concept: `Checkpoint & Revise`

## Command surface mapping

- Legacy top-level `python main.py interactive`
  - Mapped to: `python main.py workflow interactive`
- Legacy top-level `python main.py cli ...`
  - Mapped to: `python main.py pdb run ...`
- Legacy top-level `python main.py batch ...`
  - Mapped to: `python main.py pdb batch ...`
- Legacy alias `workflow clone-maturation`
  - Supported as deprecated alias of `workflow clone-checkpoint`

## Post-docking migration (legacy -> unified)

- Canonical project post-docking commands now resolve to the unified pipeline path:
  - `analyze.comparative`
  - `analyze.favorite_engine`
  - canonical stage targets (`analyze.stage.*`, interactions, visuals) via unified favorite-engine delegation
- Workflow stage targets no longer execute legacy non-canonical fallbacks; they require canonical project context.
- `post_docking_analysis.simplified_cli` remains available as a compatibility wrapper:
  - manifest-backed canonical projects auto-delegate to unified execution
  - legacy-only toggles (`--no-rmsd`, `--no-visualizations`) are ignored in manifest wrapper mode
- User-facing RMSD scope contract is now `per_complex` only:
  - global/per-protein RMSD scope toggles were removed from CLI surfaces

## Layout/profile migration

- Canonical multi-engine profile: `canonical`
  - Uses numbered output structure with unified project manifest.
- Backward-compatible profile: `docking_legacy`
  - Preserves historical tree semantics where needed.

During initialization, existing projects are inspected and either:

1. recognized as canonical, or
2. backfilled with workflow state/manifest metadata without destructive re-layout.

## Interactive behavior migration

- Universal back behavior is now part of all prompt helpers.
- Long-running operations can run in background task mode.
- Timeline and session state are recorded in `.workflow/state.json`.

## Checkpoint & Revise metadata

Each checkpoint clone now writes lineage metadata in workflow state:

- `checkpoint_id`
- `source_project_dir`
- `target_project_dir`
- `layout_profile`
- `marker_file`
- `created_at`

This metadata is designed to support future lineage and reproducibility reports.

## Compatibility guardrails

- Legacy prompts still work where possible, but new commands are preferred.
- Feature flags in workflow state allow staged rollout for new subsystems.
- Deprecated aliases remain available for migration windows and scripts.
- Stage-level delegation emits explicit run notes so users can see when canonical stage targets are executed through unified favorite-engine continuation.
