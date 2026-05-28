# Tasks: Multi-Engine HPC Docking Orchestration

**Input**: Design documents from `/specs/014-multi-engine-hpc-docking-orchestration/`
**Prerequisites**: `plan.md`, `spec.md`

**Tests**: Use `py_compile`, CLI help checks, and a synthetic staged-project
validation for pair-curation, deployment generation, and rerun-manifest flow.

## Phase 1: Setup

- [x] T001 Create the numbered branch `014-multi-engine-hpc-docking-orchestration`
- [x] T002 Scaffold `spec.md`, `plan.md`, and `tasks.md`

## Phase 2: Pair Curation

- [x] T003 Add persistent pair-curation state paths and manifest fields
- [x] T004 Extend pairlist preparation with iterative rounds and round freezing
- [x] T005 Preserve the downstream `pairlist.csv` and `pair_intent.csv` contract

## Phase 3: HPC Deployment

- [x] T006 Split runner job planning from execution
- [x] T007 Add Slurm deployment generation for GNINA, Vina, and Smina
- [x] T008 Record deployment manifests and runtime metadata in the project state
- [x] T008a Add secure HPC profile loading with ignored local JSON overrides
- [x] T008b Add engine-specific Slurm script preambles for cluster bootstrap

## Phase 4: Comparative Reruns

- [x] T009 Extend comparative analysis to emit an exhaustive rerun manifest
- [x] T010 Add deploy support for rerun-manifest-driven exhaustive jobs
- [x] T010a Add pair-level rerun filters for engine wins, affinity advantage, pair caps, and explicit allowlists

## Phase 5: Validation

- [x] T011 Run `python -m py_compile` on the touched modules
- [x] T012 Run CLI help checks for the new surfaces
- [x] T013 Run a synthetic staged-project validation covering Stage 1 and Stage 2 deployment flow
- [x] T014 Inspect the current BibAlex cluster for GNINA, Conda, Vina, and Smina availability
- [x] T015 Write public-safe deployment/bootstrap guidance and a sanitized profile template
- [x] T016 Install and verify user-local Vina and Smina CLIs on the BibAlex cluster
- [x] T017 Replace per-pair Slurm submission with one scheduler submission bundle per engine
- [x] T018 Add exhaustive-mode safety rails to both deploy and submit flows
- [x] T019 Extend the interactive workflow so comparative analysis can generate rerun manifests and the HPC stage can consume them
- [x] T020 Verify GNINA, AutoDock Vina, and Smina wrapper assumptions against upstream GitHub documentation
- [x] T021 Add the `cocrystal_plus_nonreference` pairlist mode for own-reference-plus-library screening
- [x] T022 Fix comparative `smina` parsing for `minimizedAffinity` / `minimizedRMSD` PDBQT outputs
- [x] T023 Run real-project comparative validation on `3-Docking` and confirm three-engine report generation
- [x] T024 Run real-project downstream report generation on `3-Docking`
- [x] T025 Run real-project downstream visualization generation on `3-Docking`
- [x] T025a Run real-project favorite-engine continuation on `3-Docking` through final consolidation
- [x] T025b Validate ProLIF routing on `3-Docking` and record current environment dependency status
- [x] T025c Repair stage-visualization semantics so benchmark/cross-protein figures align with `pair_intent.csv` metadata and best-pose pair units
- [x] T025d Add ligand naming map/override support and validate the corrected stage visualization output on `3-Docking`
- [x] T025e Expose protein/ligand naming prompts through the workflow interactive analysis shell and workflow analysis CLI dispatch
- [x] T025f Add comparative cross-engine same-pair visualization batches and a GNINA CNN-specific visualization style, then validate them on `3-Docking`
- [x] T025g Add optional complex-query include/exclude filtering across workflow interactive analysis, workflow CLI, simplified GNINA analysis, and comparative multi-engine analysis; persist query-scoped score subsets and validate on `3-Docking`
- [x] T025h Harden the legacy post-docking pipeline and adjacent helpers: fix the dead non-GNINA branch, reject empty input paths, respect configured GNINA directories, replace fragile package-root imports, rank best-binding summaries by affinity, skip malformed normalization rows, and refresh rerun-manifest round IDs before writing
- [x] T025i Harden medium-severity post-docking helpers: fail fast on empty complex-query score overlap, honor PLIP `max_workers`, and replace fragile pose-extractor fixed-column fallback parsing with token-aware parsing
- [x] T025j Close the remaining helper loopholes in the reviewed modules: restore the missing legacy binding-affinity import, validate/coerce PLIP config knobs, skip all-NaN best-row groups, sanitize pair allowlists, guard rerun merge mismatches, and affinity-rank PyMOL comparison inputs
- [x] T025k Harden the interactive workflows: fix legacy `quit` / cleaning-menu control flow, expose comparative alias prompting, and prevent stage-only engine picks from overwriting the recorded favorite engine
- [x] T025l Expose real guided-workflow controls for preparation force-field / pH settings, reusable pair-curation rounds, and clearer downstream-analysis labeling/order
- [x] T025m Add persistent project-level protein/ligand alias metadata, enrich `pair_intent.csv` with display-name columns, and wire those aliases into downstream naming resolution
- [x] T025n Fix receptor PDB-ID extraction in docking preparation so prefixed receptor names still resolve to canonical 4-character PDB codes
- [x] T026 Repair `7VKO_ligand_7GI_A_801.pdbqt` in the real `3-Docking` project and validate that both staged copies no longer contain invalid AutoDock atom types
- [ ] T026b Re-run the repaired `7VKO` GNINA/Smina pair and refresh the real-project comparative coverage counts from `310/311` to `311/311`
- [x] T026c Fail unsupported AutoDock ligands during ligand preparation and project materialization so bad outputs do not survive until HPC preflight
- [x] T027 Reduce favorite-engine global RMSD runtime by splitting it into resumable per-scope stages with configurable global-scope deferral
- [x] T028 Keep `ProLIF` and `py3Dmol` as optional extras, document that contract, and ensure missing installations do not fail the validated baseline workflow
