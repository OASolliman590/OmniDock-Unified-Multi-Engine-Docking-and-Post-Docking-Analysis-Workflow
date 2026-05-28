# Implementation Plan: Post-Docking Analysis Remediation & Contract Uniformity

**Branch**: `014-multi-engine-hpc-docking-orchestration` | **Date**: 2026-05-24 | **Spec**: `specs/029-post-docking-remediation/spec.md`
**Input**: Feature specification from `/specs/029-post-docking-remediation/spec.md`

## Summary

Unify post-docking contract behavior across interactive/non-interactive surfaces, remove legacy bypass routes, and bring multi-engine parity to the same hard-fail + observability contract. This implementation starts with P1 correctness: FR-001/FR-002/FR-003 already landed, and FR-004 is the active slice.

## Technical Context

**Language/Version**: Python 3.11
**Primary Dependencies**: pandas, RDKit, MDAnalysis, matplotlib/seaborn (optional analysis stack)
**Storage**: File-system artifacts + JSON/CSV manifests
**Testing**: `python -m py_compile`, repository smoke harness (`test/test_dockforge_smoke.py`)
**Target Platform**: Local Linux/macOS + HPC via Slurm wrappers
**Project Type**: Single Python project

## Current Phase

- Completed: FR-001 (`run_rmsd`/`run_visualizations` enforcement in non-GNINA context)
- Completed: FR-002 (standalone interaction aliases routed to clean contract)
- Completed: FR-003 (legacy `structure_quality` / `pymol` routes blocked with migration guidance)
- In progress: FR-004 (engine parser + detector parity, currently AutoDock4 slice)

## Validation Strategy

1. Compile checks on touched modules.
2. Focused smoke checks for:
   - non-GNINA contract enforcement,
   - legacy-route blocking,
   - non-interactive alias-to-clean routing,
   - engine detection parity including AutoDock4.
3. Follow-up parity expansion for Vina/Smina/AutoDock4 hard-fail branch tests.
