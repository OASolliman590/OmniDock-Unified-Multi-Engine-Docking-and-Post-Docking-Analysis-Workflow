# Feature Specification: Pipeline Unification and Deduplication

**Feature Branch**: `001-unify-and-optimize`  
**Created**: 2026-02-12  
**Status**: Draft  
**Input**: User description: "unify the pipeline, remove redundant or duplicate code"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Consistent PandaMap Runtime (Priority: P1)

As a researcher, I need PandaMap behavior to be consistent across standard and publication analyzers, without duplicated command logic drifting over time.

**Why this priority**: PandaMap is used in multiple analysis paths, and drift causes hard-to-debug differences.

**Independent Test**: Run both PandaMap analyzers and verify they share the same CLI detection/execution behavior.

**Acceptance Scenarios**:

1. **Given** PandaMap analyzers in `pandamap_integration.py` and `publication_pandamap.py`, **When** commands are executed, **Then** both use one shared runner implementation.
2. **Given** PandaMap CLI mode changes, **When** analyzer code runs, **Then** both analyzers detect mode identically.

---

### User Story 2 - Remove Dead Pipeline Flow (Priority: P1)

As a maintainer, I need unreachable and duplicated pipeline flow removed so there is one authoritative execution path.

**Why this priority**: Dead blocks increase maintenance risk and hide real behavior.

**Independent Test**: Run syntax checks and pipeline smoke execution to confirm behavior matches current intended flow after removing dead code.

**Acceptance Scenarios**:

1. **Given** `_generate_all_scores_csv` in `pipeline.py`, **When** code is inspected, **Then** no unreachable legacy pipeline block remains inside the method.
2. **Given** pipeline execution, **When** run in GNINA fast-path, **Then** it still completes with unchanged functional outcome.

---

### User Story 3 - Eliminate Local Duplication and Ambiguity (Priority: P2)

As a maintainer, I need duplicate helpers/import ambiguity removed in simplified pipeline to keep data-labeling behavior clear.

**Why this priority**: Ambiguous duplicate imports and repeated helper logic cause subtle regressions.

**Independent Test**: Run compile checks and smoke invocation; verify naming and output organization still work.

**Acceptance Scenarios**:

1. **Given** simplified pipeline imports, **When** module is loaded, **Then** no duplicate/overwritten utility import exists.
2. **Given** helper logic shared across components, **When** future changes are made, **Then** updates happen in one place.

### Edge Cases

- PandaMap unavailable in environment should still fail gracefully in both analyzers.
- Shared PandaMap runner must support both single-command and subcommand CLI modes.
- Pipeline refactor must not alter existing CLI arguments or output filenames unexpectedly.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST provide a shared PandaMap command runner used by both PandaMap analyzers.
- **FR-002**: System MUST remove unreachable legacy execution block from `post_docking_analysis/pipeline.py`.
- **FR-003**: System MUST preserve current GNINA fast-path behavior after refactor.
- **FR-004**: System MUST remove local duplicate/ambiguous import patterns in simplified pipeline.
- **FR-005**: System MUST pass module compile checks for all changed files.

### Key Entities *(include if feature involves data)*

- **PandaMapRunner**: Shared utility for PandaMap CLI execution and mode detection.
- **PipelineExecutionPath**: Canonical set of runtime steps in `pipeline.py` without dead branches.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: One shared PandaMap runner is referenced by both analyzers.
- **SC-002**: `pipeline.py` contains no unreachable duplicate block in `_generate_all_scores_csv`.
- **SC-003**: All modified files pass `py_compile`.
- **SC-004**: At least one simplified pipeline smoke run completes after refactor.
