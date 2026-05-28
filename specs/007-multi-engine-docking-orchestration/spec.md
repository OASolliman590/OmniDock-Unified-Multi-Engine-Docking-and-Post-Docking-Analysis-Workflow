# Feature Specification: Multi-Engine Docking Orchestration

## Context

The repository currently prepares proteins and analyzes docking outputs, but the
actual docking execution path is fragmented and GNINA-centric. The new feature
must support a canonical project layout that can drive GNINA through a
Singularity/Apptainer workflow and AutoDock Vina or Smina through conda-based
execution.

## Goals

1. Add a canonical docking project contract with shared receptors, ligands,
   metadata, and per-engine output folders.
2. Add a new top-level docking mode that can execute GNINA, Vina, and Smina
   from the same `pairlist.csv`.
3. Record per-engine commands, runtime settings, and run results in manifest
   files so downstream analysis does not need to rediscover raw engine layouts.
4. Support dry-run validation before executing external binaries.

## Out of Scope

- HPC scheduler submission integration.
- Auto-installing GNINA, Vina, Smina, or container images.
- Replacing the legacy GNINA HPC adapter for old projects.

## Functional Requirements

- FR-001: `python main.py dock run` must accept a canonical or staged docking
  project root and a comma-separated engine list. Legacy
  `python main.py dock --project-dir ...` invocations must continue to
  normalize to the `dock run` surface.
- FR-002: GNINA commands must support containerized execution with optional GPU
  flags and explicit image/binary paths.
- FR-003: Vina and Smina commands must support execution through conda
  environments or direct binaries.
- FR-004: Each engine must write outputs under the active project layout:
  `engines/<engine>/{poses,logs,scores}` for canonical projects, or the staged
  compatibility layout under `4-Docking/` (`gnina_out`, `logs`, `results`,
  `vina_out/`, and `smina_out/`) for `docking_legacy` projects.
- FR-005: Each engine run must emit a `run_manifest.json` capturing commands and
  per-job status.
- FR-006: Successful runs must emit normalized score tables suitable for
  engine-aware analysis.

## Acceptance Criteria

1. `python main.py dock run --help` works, and legacy `python main.py dock --help`
   still exposes the same GNINA, Vina, and Smina runtime options through
   argument normalization.
2. `python main.py dock run --project-dir <project> --dry-run` (or the legacy
   alias form) writes per-engine run manifests without executing binaries.
3. A real or synthetic Vina/Smina output folder can be normalized into
   `normalized_scores.csv` in the engine score directory.
4. GNINA output logs can be converted into normalized scores through the new
   engine runner path.
5. `python -m py_compile` passes for all added docking modules.
