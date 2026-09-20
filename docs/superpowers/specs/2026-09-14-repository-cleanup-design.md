# OmniDock Repository Cleanup Design

**Status:** approved for implementation by the user-provided cleanup brief
**Fixed review base:** `8d4434d3a33eb83b1e12cad82944b02c83270e47`
**Working branch:** `codex/cleanup-simplify-docs`
**Date:** 2026-09-14

## Goal

Make OmniDock easier to install, navigate, review, and maintain without weakening its scientific safeguards or breaking supported command surfaces. The cleanup must leave a smaller active documentation surface, explicit historical provenance, fewer duplicated orchestration paths, and automated checks for repository/documentation drift.

## Evidence at the fixed base

- 300 tracked files totaling 13,327,496 bytes.
- 125 Python files totaling 2,741,002 bytes and 115 Markdown files totaling 752,032 bytes.
- `closing_thesis/` contains 25 tracked files and 8,655,902 bytes. These are research/manuscript assets, not disposable build output.
- `specs/` contains 84 Markdown files from completed historical implementation cycles, but only one subdirectory is currently named `_archive`.
- `audit/` contains curated reports and executable probes plus five generated verification logs/XML files totaling 184,736 bytes.
- Five relative Markdown links are broken at the base: two references to the absent `CONTRIBUTING.md`, two references to the absent `AUTODOCK_PREPARATION_GUIDE.md`, and one audit link to an untracked probe artifact.
- `4-K_pneumoniae_/5-Post_Docking_Analysis` is a tracked symbolic link whose target is an absolute private workstation path. It is dangling in the isolated checkout and is not referenced by code or documentation.
- Active installation documents disagree on Python 3.8, 3.9, and 3.10; the package contract and CI require Python 3.10+.
- `README.md` is 834 lines and mixes current multi-engine workflow instructions with old PDB Prepare Wizard tutorials, unsupported production-readiness claims, stale module maps, and an incorrect Python 3.8 footer.
- `post_docking_analysis/` contains current README/visualization guidance beside historical redesign, decision-tree, pairlist, implementation-status, and usage documents without a clear archive boundary.
- Largest orchestration functions include `workflow.cli.build_parser` (465 lines), `workflow.interactive._run_analysis_stage` (612 lines), and `workflow.interactive._run_hpc_stage` (524 lines).
- `docking.deployment.generate_condor_deployment` and `generate_slurm_deployment` repeat project initialization, engine job planning, path mapping, payload assembly, and manifest finalization across two roughly 200-line functions.
- `autodock/prep_autodock.sh` and `autodock/prep_ligands_custom.sh` fail Bash parsing, are not packaged or called by the active code, and predate the authoritative-graph/fail-closed preparation contracts. The supported `prep_autodock_enhanced.sh` replacement is packaged and syntax-valid.
- Windows baseline: the exact full test command cannot collect three RDKit-dependent modules because Windows Application Control blocks RDKit's `rdinchi` DLL. The remaining workflow/docking/smoke set records 145 passes and two RDKit-dependent smoke failures; those failures return `rdkit_required_for_chemical_mapping`. This host limitation must not be hidden or converted into skipped assertions.

## Boundaries and parallel-work ownership

The cleanup owns repository hygiene, package/distribution metadata, active documentation and indexes, command/example verification, workflow CLI construction, and non-scientific docking deployment structure.

A separate scientific-integrity task owns behavioral corrections in ligand preparation, protein preparation, and post-docking scientific execution. Until it supplies exact file ownership, this cleanup will not change scientific algorithms in those areas. Shared public docs are owned here, but must preserve the historical scope and limits in `docs/correctness-update.md`; new scientific assessment/method documents can be linked after that task lands.

The benchmark task is paused. This cleanup will not resume it, access NMRBox, submit jobs, or incorporate uncommitted benchmark work.

## Chosen design

### 1. Add a repository-contract module

Create one small checker module with the interface `python scripts/check_repository.py`. It will report broken relative Markdown links, tracked generated verification artifacts, and tracked dangling symbolic links. Tests exercise the checker through temporary repositories/trees rather than mocking its filesystem behavior.

This module creates a reusable seam for local development and CI. The interface stays small while the implementation hides Markdown target resolution and Git tracking details.

### 2. Preserve evidence while removing generated clutter

- Remove the dangling `4-K_pneumoniae_/5-Post_Docking_Analysis` symbolic link because its tracked contents are only an absolute private path and no repository reference consumes it.
- Replace generated audit logs/XML with `audit/verification-history.md`, which records the historical commit, commands, outcomes, limitations, artifact hashes, and how to recover the originals from Git history. Keep curated audit reports, executable probes, and structured probe JSON results.
- Keep research/manuscript data. Move root-level Sertraline manuscripts and the multiscale figure pair under `closing_thesis/` to make ownership explicit; do not delete or rewrite their results. Update the bundle README to state that its referenced run manifest is not tracked and therefore the figures require independent provenance verification.
- Add `specs/README.md` declaring numbered specs historical snapshots. Avoid a large rename-only migration of all spec directories.
- Move obsolete post-docking design/status/usage documents to `post_docking_analysis/docs/_archive/`, preserving content and repairing historical references. Keep the package README concise and point it to the canonical top-level guide.
- Expand ignore rules for canonical runtime roots, generated analysis artifacts, and local configuration without ignoring source fixtures or public templates.
- Retire the two broken, unreferenced legacy AutoDock shell helpers instead of repairing a permissive chemistry path that bypasses current safeguards. Document migration to `prep_autodock_enhanced.sh` and the Python workflow, and syntax-check every remaining tracked shell script.

### 3. Deepen docking deployment

Retain `generate_condor_deployment(...)` and `generate_slurm_deployment(...)` as the public interfaces. Move shared initialization, runtime/job planning, path mapping, and project-manifest finalization behind private helpers. Scheduler-specific rendering remains in the two adapters because Condor and Slurm are real variants at that seam.

The public function signatures and generated path/layout contracts remain compatible. An additive `scheduler` field is required in both top-level deployment manifests so downstream inspection does not need to infer it from profile shape.

### 4. Isolate CLI parser construction

Move parser construction into `workflow/cli_parser.py` with a single exported `build_parser()` interface and focused private group builders. Keep `workflow.cli.build_parser` as a compatibility re-export, and leave dispatch/execution behavior in `workflow/cli.py`.

This separates command declaration from command execution without adding a second parser or changing option names/defaults. Representative parse contracts and existing smoke tests protect compatibility.

### 5. Consolidate active documentation

`README.md` becomes an accurate landing page rather than the exhaustive manual. `docs/README.md` becomes the navigation index. The canonical active guides are:

- `INSTALLATION_GUIDE.md`: Python/dependency/environment installation and verification.
- `USAGE.md`: short end-to-end CLI recipes.
- `PDB_PREPARATION_USAGE.md` and `PDBQT_PREPARATION_EXPLAINED.md`: preparation workflows and scientific decision points.
- `HPC_DEPLOYMENT_GUIDE.md`: generic Slurm and NMRBox/HTCondor use through public-safe templates; no private host/account/path facts.
- `POST_DOCKING_ANALYSIS_GUIDE.md`: current canonical analysis modes, score/RMSD interpretation, output structure, and optional stages.
- `docs/architecture.md`: module map, interfaces, state/artifact flow, and compatibility seams.
- `docs/testing.md`: local gates, RDKit/platform requirements, CI matrix, and clean-wheel checks.
- `CONTRIBUTING.md`: development workflow, standards, tests, and scope rules.
- `docs/dockforge_migration_notes.md`, `docs/correctness-update.md`, and `CHANGELOG.md`: compatibility and historical change records.

Duplicate tutorials in package-local docs become short routing pages or archived historical documents. Past release-candidate reports receive explicit historical headers; their old verification claims are not rewritten as current evidence.

## Scientific invariants that must not regress

- authoritative chemical graphs and atom mappings;
- native receptor-frame RMSD with no ligand superposition that hides displacement;
- explicit score direction and engine/scoring-function-specific comparisons;
- explicit `not_evaluable` statuses and coverage accounting;
- reference source/frame provenance;
- fail-closed execution and QC gates;
- content-verified resume/cache identity;
- atomic, concurrency-safe state changes;
- target-specific limitations and no invented benchmark, affinity, efficacy, or production-readiness claims.

## Verification

Each batch is reviewed from a clean tree and committed only by the orchestrator. Required checks are:

1. focused test-first checks for the batch;
2. `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test -q` in an environment where RDKit can load;
3. on the current Windows host, the full non-RDKit suite plus explicit recording of the RDKit policy blocker;
4. `python -m build` and `python scripts/check_wheel.py` from a clean installation environment;
5. `python scripts/check_repository.py`;
6. shell syntax for tracked shell scripts on Linux/Git Bash;
7. `python -m compileall` and `git diff --check`;
8. Linux and Windows GitHub Actions on the draft pull request.

Final review compares `8d4434d...HEAD` against this specification on separate Standards and Spec axes. The repository currently lacks `docs/agents/issue-tracker.md`; the fixed base and this in-repository specification are the review sources for this work.

## Rejected alternatives

- **Wholesale rewrite:** rejected because it would obscure the scientific corrections and make regression attribution difficult.
- **Delete all historical specs/audits:** rejected because they preserve design and scientific provenance.
- **Treat file size or missing imports as deletion proof:** rejected because CLI, shell, package resources, dynamic tools, and research artifacts require stronger evidence.
- **Rename the distribution immediately:** rejected because `pdb-prepare-wizard` and its console scripts are public compatibility surfaces. Naming can be documented without an unrequested package migration.
- **Refactor scientific post-docking implementation concurrently:** rejected while the scientific-integrity task owns those behavioral files.
