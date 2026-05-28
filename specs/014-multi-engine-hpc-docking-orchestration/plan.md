# Implementation Plan: Multi-Engine HPC Docking Orchestration

## Delivered Work

- Added persistent pair-curation rounds under docking metadata.
- Extended `prep pairlist` with iterative round persistence and explicit round
  freezing.
- Added `dock deploy` to render Slurm-ready assets from the current staged
  project layout for GNINA, Vina, and Smina.
- Added secure HPC profile loading so site-specific Slurm accounts, paths, and
  engine defaults can live in ignored local JSON files instead of tracked code.
- Added engine-specific Slurm script preambles so GNINA can bootstrap
  Apptainer/Singularity while Vina and Smina can rely on shell/Conda setup.
- Added richer comparative rerun-promotion controls so exhaustive follow-up can
  be restricted by engine wins, affinity advantage, max promoted pairs, and
  explicit pair allowlists.
- Refactored engine runners so job planning can be reused for both local runs
  and HPC deployment generation.
- Extended comparative multi-engine analysis to emit an exhaustive rerun
  manifest for a selected engine and top-ranked ligands per protein.
- Added the manifest plumbing needed to record the latest pair round,
  deployment root, pair-curation state file, and rerun manifest path.
- Added a public-safe HPC deployment guide and a sanitized BibAlex template,
  plus a local ignored `alex-bibalex` profile for immediate use.
- Verified the current BibAlex cluster state: GNINA binary and Apptainer image
  exist, `vina_dock` currently only provides the Python `vina` package, and no
  `smina` executable was found in the inspected locations.
- Installed and verified user-local Vina and Smina CLIs under `$HOME/.local/bin`
  on the BibAlex cluster, then updated the local ignored deployment profile to
  point directly at those binaries.
- Replaced per-pair Slurm submission with one scheduler bundle per engine and
  added profile-controlled `single_job` vs `slurm_array` submission modes.
- Hardened the deploy and submit commands so exhaustive mode now requires a
  rerun manifest unless the operator passes an explicit override.
- Extended the interactive workflow so comparative analysis can emit rerun
  manifests and the HPC stage can consume them safely.
- Verified the GNINA, AutoDock Vina, and available Smina GitHub documentation so
  the engine wrappers reflect the documented CLI contract instead of memory.
- Added the `cocrystal_plus_nonreference` pairlist mode so each receptor keeps
  its own cocrystal benchmark row while excluding other receptors' reference
  ligands from the expansion pool.
- Fixed `smina` score parsing so comparative analysis accepts real Smina PDBQT
  outputs that emit `REMARK minimizedAffinity` / `REMARK minimizedRMSD` instead
  of `REMARK VINA RESULT`.
- Ran real-project comparative analysis on
  `<dataset-root>/3-Docking`
  and validated unified outputs under
  `5-Post_Docking_Analysis/analysis/analyze_comparative`.
- Ran real-project downstream stage validation for the same project:
  `analyze stage reports` completed and `analyze stage visualizations`
  completed under their own analysis output roots.
- Repaired the simplified post-docking visualization stage so it now:
  - uses best-pose protein-ligand pairs for the affinity distribution figure
    instead of all raw poses,
  - consumes benchmark metadata from `pair_intent.csv` when that richer file is
    present next to `pairlist.csv`,
  - generates ligand naming maps/override templates alongside the existing
    protein naming files,
  - removes the stale `analysis/visualizations` subtree from stage-level runs by
    stopping the older hierarchical visualization writer in the affinity stage.
- Extended the workflow analysis shell and workflow CLI dispatch so protein and
  ligand naming prompts can be requested from the interactive workflow without
  falling back to the standalone simplified CLI.
- Extended comparative multi-engine analysis so it now emits:
  - same-pair cross-engine overview figures,
  - one per-protein cross-engine batch plot set,
  - a GNINA CNN-specific visualization family that separates CNN semantics from
    classical affinity plots,
  - comparative protein/ligand naming map and override files.
- Added optional complex-query filtering across post-docking entry points so
  operators can include or exclude subsets by protein, ligand, site, or tag
  from:
  - workflow interactive analysis,
  - workflow analysis CLI,
  - standalone simplified GNINA analysis,
  - direct `python -m post_docking_analysis` multi-engine analysis.
- Hardened query-scoped simplified GNINA analysis so filtered runs now persist:
  - `filtered_pairlist.csv`,
  - a query-scoped `all_scores.csv`,
  - `complex_query.txt` and `complex_query_summary.txt`,
  which prevents downstream hierarchical and report stages from reloading the
  full unfiltered dataset.
- Hardened the remaining helper edge cases so they now:
  - import `analyze_binding_affinities` into the legacy standard pipeline path
    instead of relying on an undefined symbol,
  - coerce invalid PLIP `max_workers` values to a safe fallback instead of
    aborting the stage,
  - honor PLIP `output_formats` instead of always forcing XML/TXT/PNG command
    flags,
  - skip all-NaN affinity groups during comparative best-row selection,
  - filter blank and literal `nan` rows out of pair allowlists,
  - raise explicitly when promoted rerun tags cannot merge back to the
    canonical pairlist,
  - affinity-rank PyMOL comparison inputs instead of relying on filesystem
    glob order.
- Hardened the legacy post-docking pipeline and adjacent helpers so they now:
  - reject empty `input_dir` instead of silently resolving to the current
    working directory,
  - execute a real standard-analysis path for non-GNINA runs instead of
    returning success without analysis,
  - use package-relative imports for visualization and pose-extraction helpers,
  - respect configured GNINA output directories during best-pose extraction,
  - rank best-binding summaries by affinity rather than file modification time,
  - skip malformed normalization rows with missing affinities instead of
    aborting an entire engine frame,
  - reload the project manifest before writing exhaustive rerun manifests so
    `round_id` stays current.
- Hardened medium-severity analysis helpers so they now:
  - fail fast when a complex-query subset has no surviving overlap between the
    requested canonical pair tags and usable score rows,
  - honor PLIP `max_workers` through a real worker pool instead of a dead
    sequential loop,
  - tolerate irregular PDBQT atom spacing in the pose-extractor fallback path
    by using token parsing when fixed-column slicing is not reliable.
- Added persistent project-level alias metadata so docking projects now keep:
  - `protein_aliases.csv`
  - `ligand_aliases.csv`
  under docking metadata,
  and regenerated `pair_intent.csv` now carries:
  - `pdb_id`
  - `protein_display_name`
  - `ligand_display_name`
  - `cocrystal_ligand_display_name`
  for downstream reporting without analysis-local re-entry.
- Hardened receptor PDB-ID extraction in the docking preparation path so
  prefixed receptor file names such as `VEGFR2_4AG8_cleaned.pdbqt` resolve to
  `4AG8` rather than a leading non-PDB token.
- Hardened the ligand-preparation and project-materialization paths so
  unsupported AutoDock atom types now fail before HPC deployment:
  - `run_enhanced_preparation()` rejects invalid prepared ligand outputs and
    writes a preparation validation report
  - `Prepare docking folders` refuses to materialize an active docking set if
    the selected ligands are not actually AutoDock-ready

## Validation Path

1. `python -m py_compile main.py docking/models.py docking/project_layout.py docking/preparation/project_builder.py docking/preparation/pairlist_builder.py docking/preparation/pair_curation.py docking/runners/base.py docking/deployment.py docking/cli.py workflow/models.py workflow/cli.py workflow/execution.py post_docking_analysis/multi_engine_pipeline.py`
2. CLI help checks for:
   - `python main.py dock --help`
   - `python main.py dock deploy --help`
   - `python main.py prep pairlist --help`
   - `python main.py analyze comparative --help`
3. Synthetic staged-project validation covering:
   - iterative pair-curation round creation
   - round freezing into `pairlist.csv`
   - Stage 1 Slurm deployment generation
   - comparative rerun-manifest emission
   - Stage 2 exhaustive deployment generation from the rerun manifest
4. Profile-aware Slurm rendering validation covering:
   - GNINA GPU deployment with Apptainer preamble and GPU partition defaults
   - GNINA CPU deployment with CPU partition defaults and `--cpu` runtime flag
   - Vina/Smina CPU deployment with direct user-local binary paths
5. Remote BibAlex read-only inspection covering:
   - `apptainer`, `singularity`, and `conda` availability
   - existing GNINA binary/image paths
   - presence or absence of `vina` and `smina` CLIs
6. Comparative rerun validation covering:
   - `--winner-only`
   - `--min-affinity-advantage`
   - `--max-rerun-pairs`
   - `--pair-allowlist`
7. Scheduler-submission validation covering:
   - one scheduler submission per engine
   - visible profile-driven submit mode (`single_job` or `slurm_array`)
8. Exhaustive safety validation covering:
   - deploy rejection without `--from-rerun-manifest`
   - submit rejection when a local exhaustive bundle does not record a rerun
     manifest
9. Real-project comparative validation covering:
   - corrected three-engine ingestion on `3-Docking`
   - engine coverage counts from the generated reports
   - explicit recording of the single residual bad pair
10. Real-project downstream validation covering:
   - favorite-engine continuation enters the deep GNINA pipeline on real data
   - stage-specific reports generation completes on real data
   - stage-specific visualization generation completes on real data
11. Real-project visualization repair validation covering:
   - stage visualization rerun to `analyze_stage_visualizations_aligned`
   - benchmark-aware figure generation from `pair_intent.csv`
   - ligand naming map/override file emission
   - no duplicate `analysis/visualizations` subtree in the new stage output
12. Real-project comparative visualization validation covering:
   - cross-engine same-pair heatmap generation
   - cross-engine disagreement and distribution plots
   - one per-protein batch plot for each detected protein target
   - non-empty GNINA CNN-specific visualization output
13. Query-scoped validation covering:
   - `python -m post_docking_analysis --analysis-mode comparative_all_engines`
     with `--complex-query 'protein=2FVD;ligand=Sorafenib,MS3'`
   - `python -m post_docking_analysis.simplified_cli` with the same query
   - confirmation that the comparative output contains only the selected two
     protein-ligand tags across the three engines
   - confirmation that the simplified GNINA path parses only the two selected
     log files and reports `40` poses / `2` protein-ligand pairs downstream
14. Legacy-pipeline hardening validation covering:
   - empty-input rejection for `PostDockingAnalysisPipeline(input_dir='')`
   - non-GNINA standard-path execution returning failure when no complexes
     exist, rather than a false success
   - synthetic multi-engine normalization proving that a missing
     `vina_affinity` row is skipped while valid rows still normalize
   - synthetic rerun-manifest generation proving `round_id` is re-read from the
     on-disk manifest immediately before writing
15. Medium-severity helper validation covering:
  - synthetic multi-engine query mismatch proving the pipeline raises a clear
    `ValueError` instead of keeping a stale pairlist/score combination
  - synthetic PLIP directory analysis proving `max_workers=3` processes a
    three-file directory successfully through the worker pool
  - parser checks proving the pose-extractor fallback accepts both fixed-width
    and ragged tokenized PDBQT ATOM lines
16. Remaining helper hardening validation covering:
   - a standard legacy binding-affinity call proving the missing-symbol path
     is closed
   - PLIP command construction proving `output_formats=['xml']` emits `-x`
     without forced `-t/-p`
   - invalid PLIP `max_workers` values falling back to `1` instead of raising
   - NaN-only affinity groups being skipped safely during comparative
     best-row selection
   - allowlist CSV rows with blanks/`nan` tokens being ignored
   - exhaustive rerun generation raising clearly when promoted tags cannot
     merge back to the canonical pairlist
17. RMSD resumability validation covering:
   - a synthetic per-scope RMSD rerun proving a cached summary is reused when
     the input signature matches
   - a synthetic global scope proving the run is marked `deferred` when the
     pose count exceeds the configured threshold
   - CLI help checks proving the new RMSD controls are exposed from
     `python -m post_docking_analysis` and `python main.py analyze favorite-engine`
18. Ligand-quality remediation validation covering:
   - repair of `7VKO_ligand_7GI_A_801.pdbqt` in the real `3-Docking` project
   - validation proving both staged copies no longer contain invalid AutoDock
     atom types
   - `ligand_repair_report.json` emission under the real project metadata tree
19. Interactive-workflow hardening validation covering:
   - legacy cleaning menu option `1` prompting for residue names rather than
     treating the literal string `1` as a residue
   - `quit` being accepted before PDB-ID format validation in both single- and
     multi-PDB modes
   - comparative analysis help surfacing protein/ligand prompt flags
   - pair-curation round controls becoming reachable from the guided workflow
   - stage-only analysis no longer updating the recorded favorite engine
20. Alias-metadata validation covering:
   - synthetic staged-project generation of `protein_aliases.csv` and
     `ligand_aliases.csv`
   - regenerated `pair_intent.csv` containing persisted display-name columns
   - downstream naming helpers resolving those persisted aliases without
     requiring analysis-output override CSVs
21. Unsupported-ligand hardening validation covering:
   - a real Pd-containing ligand proving `validate_prepared_ligand_pdbqt()`
     flags `invalid_atom_types`
   - preparation/materialization code paths compiling with the new early-fail
     guards in place

## Expected Deliverable

1. A selective-pairing workflow that can evolve across rounds without changing
   the runtime pairlist schema
2. A deterministic Slurm deployment layer for the three supported engines
3. A secure local-profile mechanism for cluster-specific deployment settings
4. A concrete bridge from comparative analysis into targeted exhaustive reruns

## Real-Project Findings

- Corrected comparative engine coverage on `3-Docking` is:
  - `gnina`: `310/311` complexes
  - `vina`: `311/311` complexes
  - `smina`: `310/311` complexes
- Comparative winners from `best_engine_per_complex.csv` are:
  - `gnina`: `218`
  - `smina`: `66`
  - `vina`: `27`
- The one residual missing tag is:
  - `7VKO_cleaned.pdbqt_site_1_7VKO_ligand_7GI_A_801.pdbqt`
- Root cause of that missing tag is not analysis code. The prepared ligand file
  still contains invalid atom types (`CG0` / `G0`), which leaves GNINA output
  empty and Smina output zero-byte for that pair.
- Real-project favorite-engine continuation for `gnina` is operational through
  complex creation, hierarchical analysis, polypharmacology, per-complex RMSD,
  global-best-pose RMSD, PandaMap generation, and final consolidation.
- Optional environment-backed stages were validated as wiring-correct but
  dependency-limited in the local environment:
  - `ProLIF` stage skips because `ProLIF` is not installed.
  - `py3Dmol` stage skips because `py3Dmol` is not installed.
  - These two stages are now treated as optional extras rather than baseline
    required dependencies.
- RMSD scalability is now controlled explicitly:
  - per-scope summaries persist checkpoint state keyed by an input signature
  - repeated runs reuse completed scope summaries when the input set matches
  - the global best-pose scope now defers automatically when pose count exceeds
    the configured threshold unless the user explicitly forces it
- The upstream source ligand for the one residual missing pair has now been
  repaired:
  - `3-Preparation/4-Prepared_Ligand/7VKO_ligand_7GI_A_801.pdbqt`
  - `4-Docking/ligands/7VKO_ligand_7GI_A_801.pdbqt`
  - both staged copies now validate cleanly and no longer contain `CG0` / `G0`
  - `ligand_repair_report.json` was written under the real project metadata
- The comparative engine counts remain `310/311` for `gnina` and `smina`
  until the corrected `7VKO` pair is rerun and the real project outputs are
  refreshed.
- Interactive workflow hardening is now in place:
  - legacy `interactive_pipeline.py` accepts `quit` cleanly before PDB-ID
    validation and exposes a deterministic cleaning menu path
  - `workflow/interactive.py` now exposes force-field / pH prompts during
    preparation and reusable pair-curation round actions during pairlist
    building
  - comparative analysis can now prompt for protein and ligand aliases, not
    just generate override CSV templates
  - stage-only analysis no longer overwrites the project's favorite-engine
    selection
- Stage-visualization repair on the same real project produced a corrected
  histogram baseline:
  - all raw poses: `6200`, mean affinity `+1.70`, median `-4.65`
  - best pose per protein-ligand pair: `310`, mean affinity `-4.62`, median
    `-6.78`
- The repaired visualization rerun created:
  - benchmark-aware `comparative_redocking.png` with `11` detected cocrystal
    benchmark rows (the `5ZE3` benchmark remains absent because its cocrystal
    ligand selection is still unresolved in source data)
  - ligand naming files `ligand_name_mapping.csv` and
    `ligand_name_overrides.csv`
  - a single top-level `visualizations/` folder with no
    `analysis/visualizations/` duplicate subtree
- Query-scoped validation on the same real project now works end to end:
  - simplified GNINA analysis with
    `protein=2FVD;ligand=Sorafenib,MS3` parsed only `2` log files, generated
    `40` score rows, loaded `1` protein / `2` ligands downstream, and created
    `2` complexes
  - comparative multi-engine analysis with the same query wrote only the two
    selected tags under `/tmp/pdbw_query_compare`, including a `6`-row
    `cross_engine_pair_best_scores.csv` and a single per-protein batch plot
- Legacy hardening checks also passed:
  - empty `input_dir` now fails validation explicitly instead of binding to CWD
  - a Vina-only run on an empty directory now returns `False` from the real
    standard pipeline path instead of returning success with zero analysis
  - synthetic normalization skipped a row with missing `vina_affinity` while
    preserving valid rows
  - synthetic exhaustive rerun manifest generation used the updated on-disk
    `latest_pair_round` value (`round_new`) rather than stale init-time state
- Medium-severity helper checks also passed:
  - a synthetic comparative query mismatch now raises
    `ValueError: Complex query matched canonical pairs, but none of those tags
    had usable engine scores`
  - synthetic PLIP directory analysis with `max_workers=3` returned `3/3`
    successful results
  - the pose-extractor fallback parser now converted both fixed-width and
    ragged tokenized PDBQT ATOM records into valid PDB lines
