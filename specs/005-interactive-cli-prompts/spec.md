# Feature Specification: Interactive Simplified CLI Prompts

## Context

The current simplified post-docking CLI expects required paths via command
arguments and exits if they are missing. Users requested an interactive command
line flow where missing paths are prompted in-session, plus explicit protein
name prompts for each detected PDB target.

## Goals

1. Add interactive CLI prompting for missing required inputs:
   - project root mode (`--project-dir`) or explicit folder mode
   - SDF folder, log folder, receptors folder, output folder
   - optional pairlist file path
2. Keep non-interactive script behavior deterministic:
   - fail with clear errors when required args are missing and no TTY
3. Prompt protein naming per detected PDB/receptor during pipeline setup:
   - one prompt per detected receptor/PDB
   - press Enter to keep current default
4. Keep backward compatibility with existing flags and output behavior.

## Out of Scope

- Replacing current auto-detection logic (`detect_gnina_layout`)
- Changing analysis/visualization computation stages
- Adding GUI prompts

## Acceptance Criteria

1. Running `python -m post_docking_analysis.simplified_cli` in a terminal with
   missing required args enters guided prompts instead of immediate failure.
2. In non-interactive contexts (no TTY), missing required args produce a clear
   error and non-zero exit.
3. During interactive runs, users are prompted for protein names for each
   detected PDB/receptor, with defaults shown.
4. Existing explicit invocation with full args still works unchanged.
5. `python -m compileall` passes for touched files.
