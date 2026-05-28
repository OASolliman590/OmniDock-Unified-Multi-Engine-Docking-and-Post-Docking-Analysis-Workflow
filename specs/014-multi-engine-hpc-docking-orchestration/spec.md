# Feature Specification: Multi-Engine HPC Docking Orchestration

**Feature Branch**: `014-multi-engine-hpc-docking-orchestration`  
**Created**: 2026-03-09  
**Status**: In Progress  
**Input**: User request to mature docking project preparation for selective pairing, HPC deployment readiness, and cross-engine post-docking compatibility.

## Context

The project already supports staged PDB retrieval, protein and ligand preparation,
project-aware pairlist generation, local multi-engine docking, and comparative
post-docking analysis. The missing operational layer is the bridge from a
prepared project to an HPC-ready multi-engine docking campaign, plus a concrete
way to promote Stage 1 comparative results into targeted exhaustive reruns.

## Goals

1. Add an iterative pair-curation workflow on top of the current pairlist
   contract without breaking downstream analysis.
2. Generate HPC-ready Slurm assets for GNINA, Vina, and Smina from the existing
   project layout.
3. Keep engine outputs normalized and compatible with comparative and
   favorite-engine post-docking analysis.
4. Turn comparative analysis into an actionable rerun-manifest producer for
   exhaustive follow-up docking.
5. Keep site-specific HPC details out of tracked repository files while still
   allowing reproducible deployment generation from local ignored profiles.

## Functional Requirements

- FR-001: `prep pairlist` MUST support reusable pair-curation rounds and round
  freezing into `pairlist.csv` and `pair_intent.csv`.
- FR-002: Pair-curation state MUST persist under the existing project layout and
  be referenced from the project manifest.
- FR-003: `dock deploy` MUST generate deterministic Slurm assets for any subset
  of `gnina`, `vina`, and `smina` without submitting jobs.
- FR-004: Deployment artifacts MUST preserve engine runtime metadata and expected
  outputs in manifest files.
- FR-004a: `dock deploy` MUST support loading site-specific runtime and Slurm
  defaults from a local JSON profile outside tracked public templates.
- FR-004b: Generated Slurm scripts MUST support engine-specific shell preambles
  such as `module load apptainer` and shell bootstrap via `~/.bashrc`.
- FR-005: `analyze comparative` MUST be able to emit an exhaustive rerun
  manifest that selects a target engine and top candidates per protein.
- FR-005a: Comparative rerun promotion MUST support optional pair-level filters
  for engine-win enforcement, minimum affinity advantage, global pair caps, and
  explicit pair allowlists.
- FR-006: The rerun manifest MUST be directly consumable by `dock deploy --mode
  exhaustive --from-rerun-manifest ...`.
- FR-007: Slurm deployment submission MUST happen as one scheduler submission
  per engine rather than one `sbatch` call per docking pair.
- FR-007a: The deployment layer MUST support both `single_job` and
  `slurm_array` submission modes, with cluster profiles allowed to choose the
  default.
- FR-008: `dock deploy` and `dock submit` MUST refuse exhaustive full-pair runs
  by default unless the user provides a rerun manifest or an explicit override.
- FR-009: The interactive workflow MUST expose the comparative-to-exhaustive
  rerun path and surface the same safeguards as the CLI.
- FR-010: Comparative analysis MUST correctly normalize `smina` PDBQT outputs
  that use `REMARK minimizedAffinity` and `REMARK minimizedRMSD`.
- FR-011: Real-project validation MUST record residual data-quality gaps when a
  ligand file produces empty or unparsable engine outputs so downstream
  comparison is not misreported as fully complete.
- FR-012: The GNINA/simplified visualization stage MUST use best-pose
  protein-ligand pairs rather than all raw poses for its global affinity
  distribution figure.
- FR-013: The visualization stage MUST consume benchmark/reference metadata
  from `pair_intent.csv` when available so cocrystal/reference figures work for
  current pairlist modes such as `cocrystal_plus_nonreference`.
- FR-014: The visualization stage MUST support ligand display-name mapping via
  generated mapping/override CSV files, analogous to the existing protein
  naming layer.
- FR-015: Stage-level visualization outputs MUST not generate a duplicate
  `analysis/visualizations/` subtree when the newer top-level visualization
  stage is the canonical output.
- FR-016: Comparative multi-engine analysis MUST emit a dedicated
  cross-engine visualization batch that compares the same protein-ligand pair
  across `gnina`, `vina`, and `smina`, including per-protein batch plots.
- FR-017: Comparative multi-engine analysis MUST render GNINA `cnn_affinity`
  using a visually distinct plot family from classical affinity plots so the
  different score semantics are explicit.
- FR-018: Post-docking analysis MUST support an optional complex-query filter
  that can include or exclude subsets by protein, ligand, site, or tag across
  comparative analysis, favorite-engine continuation, and simplified GNINA
  analysis.
- FR-018a: The complex-query filter MUST be available from the workflow
  interactive shell, workflow CLI, standalone simplified CLI, and direct
  `python -m post_docking_analysis` entry point.
- FR-018b: Query-scoped runs MUST persist the filtered subset as the analysis
  source of truth for that run so downstream reports and plots do not silently
  reload the full unfiltered score set.
- FR-019: The legacy post-docking pipeline MUST reject empty input paths,
  execute a real standard-analysis path for non-GNINA runs, respect configured
  GNINA output directories during pose extraction, and avoid package-relative
  import failures when invoked outside the package directory.
- FR-020: Query-scoped analysis MUST fail explicitly when the requested
  protein/ligand subset has no overlap with usable score tags, instead of
  silently drifting into mismatched pairlist/score state.
- FR-021: PLIP directory analysis MUST honor its advertised `max_workers`
  setting instead of always running sequentially.
- FR-022: The pose-extractor fallback path MUST tolerate irregular PDBQT atom
  spacing by falling back from fixed-column slicing to token parsing.
- FR-023: Remaining post-docking helper edge cases MUST be hardened so that:
  - the legacy standard pipeline can execute binding-affinity analysis without
    missing-symbol failures,
  - PLIP honors configured output formats and tolerates non-numeric
    `max_workers` values,
  - multi-engine best-row selection skips all-NaN affinity groups and pair
    allowlists ignore empty/`nan` rows,
  - PyMOL comparative scenes are affinity-ranked rather than filesystem-order
    dependent.
- FR-024: RMSD analysis MUST support resumable per-scope checkpoints and MUST
  be able to defer the global best-pose scope above a configurable pose-count
  threshold unless the user explicitly forces it.
- FR-025: `ProLIF` and `py3Dmol` MUST remain optional downstream extras.
  Missing installations MUST record skipped stages rather than failing the
  validated baseline pipeline.
- FR-026: Ligand quality control MUST reject invalid AutoDock atom types in
  prepared ligand PDBQT files and MUST repair the known `CG0` / `G0` artifact
  pattern before a repaired ligand is accepted as valid.
- FR-027: The interactive workflows MUST expose their real preparation and
  pair-curation controls, handle `quit` cleanly before PDB-ID validation, allow
  comparative alias prompting, and MUST NOT silently mutate the recorded
  favorite-engine when a user runs a stage-specific analysis with a temporary
  engine selection.
- FR-028: Project metadata MUST support persistent protein and ligand display
  aliases before docking, and pair-intent generation MUST carry those aliases
  forward into `pair_intent.csv` so downstream reporting and visualization can
  use them without requiring analysis-local override files.
- FR-028a: `prep pairlist` MUST generate stable project-level alias templates
  under docking metadata and MAY prompt for alias edits before writing pair
  metadata.
- FR-028b: Pair-intent generation and downstream naming helpers MUST prefer the
  persisted project aliases over filename heuristics when those aliases are
  present.

## Acceptance Criteria

1. A staged project can store and freeze a named pair-curation round.
2. Stage 1 Slurm assets can be generated for one, two, or three engines from
   the same project contract.
3. Public repo templates remain sanitized while a local ignored profile can
   drive a usable BibAlex deployment bundle.
4. Comparative analysis can emit a rerun manifest from normalized multi-engine
   score tables.
5. Stage 2 exhaustive Slurm assets can be generated from that rerun manifest
   without reshaping files by hand.
6. Submitting a deployment bundle launches one scheduler job per engine, not one
   Slurm job per pair.
7. Attempting an exhaustive deployment or submission without a rerun manifest is
   blocked unless the user explicitly overrides it.
8. Comparative analysis on a real staged project includes `gnina`, `vina`, and
   `smina` when all three engines have usable raw outputs.
9. Favorite-engine continuation on a real staged project can enter the deep
   downstream stack, with any remaining long-running stages or data failures
   recorded explicitly in the validation notes.
10. Real-project stage visualization validation produces non-empty benchmark and
    cross-protein figures, exposes ligand naming override files, and avoids a
    duplicate `analysis/visualizations/` output subtree.
11. Comparative analysis on a real staged project writes non-empty
    cross-engine same-pair figures, per-protein cross-engine batch plots, and a
    GNINA CNN-specific visualization set under the comparative output root.
12. A query-scoped post-docking run can limit analysis to explicit
    protein-ligand subsets, and the resulting reports, complexes, and score
    tables only include the requested subset.
13. The interactive workflows expose force-field / pH preparation controls,
    reusable pair-curation round controls, and comparative alias prompting, and
    stage-only engine picks do not overwrite the recorded favorite engine.
13. Large RMSD runs can resume completed per-scope checkpoints and can defer
    the global best-pose scope when the dataset exceeds the configured
    threshold.
14. Missing `ProLIF` or `py3Dmol` does not fail the baseline downstream run; the
    pipeline records a skipped optional stage instead.
15. Repaired prepared ligand PDBQT files no longer contain invalid AutoDock
    atom types such as `CG0` / `G0`.
16. A staged project can persist protein and ligand aliases in project
    metadata, regenerate `pair_intent.csv` with `protein_display_name` /
    `ligand_display_name` columns, and downstream naming helpers resolve those
    aliases automatically.
