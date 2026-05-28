# Implementation Plan: Pipeline Unification and Deduplication

**Branch**: `001-unify-and-optimize` | **Date**: 2026-02-12 | **Spec**: `specs/002-pipeline-unification/spec.md`  
**Input**: Feature specification from `specs/002-pipeline-unification/spec.md`

## Summary

Refactor duplicated runtime internals by centralizing PandaMap CLI execution into one shared helper, removing unreachable duplicate pipeline flow from `pipeline.py`, and cleaning local duplication/ambiguity in simplified pipeline imports/helpers.

## Technical Context

**Language/Version**: Python 3.9+  
**Primary Dependencies**: pandas, pathlib, subprocess, existing post-docking modules  
**Storage**: File-based outputs under pipeline output directories  
**Testing**: `py_compile` + simplified CLI smoke run  
**Target Platform**: macOS/Conda  
**Project Type**: Single Python package  
**Performance Goals**: No regression in current pipeline runtime behavior  
**Constraints**: Optional tools (PandaMap) may be absent; behavior must remain graceful  
**Scale/Scope**: Internal refactor across `post_docking_analysis` modules

## Constitution Check

- Maintain backward-compatible CLI behavior.
- Keep optional-tool graceful degradation.
- Avoid introducing destructive filesystem behavior.

## Project Structure

### Documentation (this feature)

```text
specs/002-pipeline-unification/
├── spec.md
└── plan.md
```

### Source Code (repository root)

```text
post_docking_analysis/
├── pipeline.py
├── simplified_pipeline.py
├── pandamap_integration.py
├── publication_pandamap.py
└── [new] pandamap_runner.py
```

**Structure Decision**: Add a small shared PandaMap utility module and refactor existing analyzers to use it; clean dead/duplicate flow in pipeline core.

## Phases

1. Add shared PandaMap runner module (`pandamap_runner.py`).
2. Refactor `pandamap_integration.py` and `publication_pandamap.py` to use shared runner.
3. Remove unreachable legacy flow block from `pipeline.py`.
4. Remove import ambiguity / small duplicate patterns in `simplified_pipeline.py`.
5. Validate via compile checks and one smoke run.
