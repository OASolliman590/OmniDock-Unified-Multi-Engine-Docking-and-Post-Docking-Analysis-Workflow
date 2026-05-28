# Feature Specification: Unified Workflow Orchestrator

## Context

The repository now supports PDB preparation, docking preparation, docking
execution, and post-docking analysis across GNINA, Vina, and Smina. The
interactive entrypoint must be docking-first, and the grouped CLI must expose
the same staged workflow directly.

## Goals

1. Make `python main.py workflow interactive` start with create/resume/inspect
   docking project flow.
2. Keep one grouped CLI with `workflow`, `pdb`, `prep`, `dock`, and `analyze`
   command groups.
3. Persist workflow state in `.workflow/state.json`.
4. Support direct access to docking, analysis, and post-docking subfunctions.
5. Preserve legacy aliases such as `interactive`, `cli`, `batch`,
   `prepare-docking`, and `dock`.

## Functional Requirements

- FR-001: `python main.py`, `python main.py interactive`, and
  `python main.py workflow interactive` must route to the workflow shell.
- FR-002: The interactive shell must begin with:
  - create a new docking project
  - resume an existing docking project
  - inspect project status
- FR-003: Creating a new project must ask for working directory, project name,
  and engine panel choice including `All Panels`.
- FR-004: Interactive multi-selection surfaces must use checkbox-style prompts
  rather than numeric-only menus.
- FR-005: The grouped CLI must expose:
  - `python main.py pdb collect ...`
  - `python main.py prep pairlist ...`
  - `python main.py prep project ...`
  - `python main.py workflow init --layout-profile docking_legacy`
- FR-006: Workflow state must record the staged project directories, selected
  engines/panels, and pairlist mode.
- FR-007: The workflow shell must expose direct access to post-docking
  subfunctions including comparative analysis, favorite-engine continuation,
  stage-only analysis, interaction tools, and visualization tools.

## Acceptance Criteria

1. `python main.py --help` shows the grouped command surface.
2. `python main.py workflow interactive --help`,
   `python main.py pdb collect --help`, and
   `python main.py prep pairlist --help` all work.
3. `workflow init --layout-profile docking_legacy` creates the staged project
   structure and workflow state.
4. The interactive shell uses checkbox/select prompts through `questionary`.
5. Running project init, collection, pairlist building, or docking updates
   `.workflow/state.json`.
6. Legacy aliases still resolve to the new grouped CLI.
