# Feature Specification: Engine Scope Filter

**Feature Branch**: `023-engine-scope-filter`
**Created**: 2026-04-01
**Status**: Draft
**Input**: User need to selectively analyze one or a subset of engines in a multi-engine docking project, rather than always running full multi-engine consensus or always running a forced single-engine override.

---

## Summary

Spec 022 (Engine Detection and Solo Routing) handles automatic detection and zero-config routing. This spec handles the complementary case: **the user has a multi-engine project and explicitly wants to choose which engines to include in a specific analysis run**.

This is distinct from the 022 `--engine` override (which forces single-engine on one named engine). The engine scope filter allows:
- **Single-from-multi**: run one engine's solo profile on a project that has multiple valid engines (e.g., "I only want GNINA's CNN-ranked results today")
- **Subset multi-engine**: run 2 of 3 engines through the multi-engine consensus path (e.g., "exclude Smina, compare only GNINA vs Vina")
- **Per-session scope persistence**: each session records which engines were in scope so reruns and comparisons are traceable

The engine scope decision is separate from detection: detection determines what is *available*, scope filter determines what is *used* for this session.

---

## Problem

### Current State
- Detection (spec 022) auto-routes to multi-engine when N ≥ 2.
- `--engine <name>` override exists but is all-or-nothing: forces single-engine mode.
- No way to say "compare GNINA and Vina but not Smina" without manually editing the pairlist or scores CSV.
- No way to re-run analysis with a different engine subset without a full new docking project.

### Concrete User Scenarios
1. Smina ran on all pairs but produced low coverage (< 50%) — user wants GNINA+Vina consensus only, skipping Smina.
2. User wants to generate GNINA solo (CNN-ranked) results for MD candidate selection while also having full 3-engine consensus results in the same project.
3. User wants to quickly compare how Vina-only ranking differs from 3-engine consensus — without re-running docking.
4. HPC job partially failed for one engine — user wants to continue analysis with the engines that completed.

---

## Goals

1. Allow user to specify an engine scope (`--engines gnina,vina`) per analysis session.
2. When scope reduces to one engine, automatically activate that engine's solo profile (spec 022).
3. When scope is a subset ≥ 2, run multi-engine consensus using only the scoped engines.
4. Record scope in session metadata and all output artifacts.
5. Interactive engine picker prompt in interactive mode when multiple engines are available.
6. Session-level scope does not modify the project's detection cache or pairlist.

## Non-Goals

- Changing which engines ran docking (this is post-analysis scoping only).
- Merging or re-weighting scores across differently-scoped sessions.
- Per-ligand engine selection (always per-session).

---

## User Stories

### US1: Exclude Low-Coverage Engine
As a scientist, I want to exclude an engine with poor coverage from consensus analysis so low-quality data does not corrupt my ranking.

### US2: GNINA Solo from Multi-Engine Project
As a medicinal chemist, I want CNN-ranked GNINA solo results from a 3-engine project for MD candidate selection without rerunning docking or creating a new project.

### US3: Pairwise Comparison
As a methods analyst, I want to run GNINA+Vina consensus and GNINA solo as two separate sessions in the same project to compare ranking stability.

### US4: Scope Traceability
As a reviewer, I want session outputs to clearly state which engines were included and excluded so I can reproduce the analysis.

### US5: Partial Failure Recovery
As an HPC operator, I want to run analysis using only the engines that completed successfully, without waiting for or discarding failed engine results.

---

## Functional Requirements

### Engine Scope Selection

- **FR-001**: The pipeline MUST accept an `--engines <comma-separated list>` CLI argument that restricts the session to the named engines.
- **FR-002**: In interactive mode, when the detected engine set has N ≥ 2 engines, the pipeline MUST present an engine picker prompt with checkboxes, defaulting to all valid engines selected.
- **FR-003**: The engine picker MUST display per-engine coverage percentage (from spec 022 detection report) alongside each engine name so users can make informed exclusion decisions.
- **FR-004**: At least one engine MUST remain selected; the pipeline MUST prevent deselecting all engines.
- **FR-005**: The engine scope is session-local: it MUST NOT modify `project_manifest.json["detected_engines"]` or the pairlist.

### Scope-to-Profile Routing

- **FR-006**: When the scoped engine set has exactly one engine, the pipeline MUST activate that engine's solo profile (spec 022 FR-010 through FR-034) for the session, identical to a native single-engine project.
- **FR-007**: When the scoped engine set has N ≥ 2 engines, the pipeline MUST run multi-engine consensus using only the scoped engines. Excluded engines MUST NOT contribute scores, ranks, or agreement fractions.
- **FR-008**: The `engine_support_count` in top-pose atlas MUST reflect the scoped engine count, not the detected engine count. A note `scoped_engine_count: N of M` MUST be included in the manifest.

### Exclusion Reasons

- **FR-009**: When `--engines` or the interactive picker excludes one or more valid engines, the session metadata MUST record each excluded engine with an `exclusion_reason` field. Accepted values: `user_excluded`, `low_coverage`, `parse_failed`, `cli_override`.
- **FR-010**: Excluded engines with valid output files MUST appear in the detection report as `excluded_this_session` (not `absent`) so their data is preserved for future sessions.

### Output Artifacts

- **FR-011**: All output files (hit tables, top-pose atlas, consensus explainability, session reports) MUST include a `engines_in_scope` field listing which engines contributed.
- **FR-012**: `START_HERE.md` MUST display which engines were in scope and which were excluded in the session summary section.
- **FR-013**: If a prior session in the same project used a different engine scope, `START_HERE.md` MUST note the discrepancy so users are aware that sessions are not directly comparable.

### Saved Scope Presets

- **FR-014**: The pipeline MUST support named engine scope presets saved in `project_manifest.json["engine_presets"]`:
  ```json
  {
    "engine_presets": {
      "cnn_only": ["gnina"],
      "vina_family": ["vina", "smina"],
      "full": ["gnina", "vina", "smina"]
    }
  }
  ```
- **FR-015**: `--engine-preset <name>` CLI flag MUST load the named preset and apply it as the session scope.
- **FR-016**: The interactive picker MUST offer a "load preset" option.

---

## Interaction with Spec 022

| Scenario | Spec 022 Behavior | Spec 023 Behavior |
|---|---|---|
| Project has 1 valid engine | Auto-routes to solo | N/A (no selection needed) |
| Project has 3 valid engines, no filter | Auto-routes to multi | All 3 in scope (default) |
| Project has 3 valid engines, `--engines gnina` | N/A | Routes to GNINA solo profile |
| Project has 3 valid engines, `--engines gnina,vina` | N/A | Routes to 2-engine multi-engine |
| Project has 3 valid engines, `--engine gnina` (022 override) | Forces single, ignores detection | Not used; 023 handles this case |

**Note**: `--engine` (spec 022, single value) and `--engines` (spec 023, list) are distinct flags. `--engine` is a hard override that bypasses detection; `--engines` is a filter applied after detection.

---

## Interactive Prompt Design

```
? Engines detected in this project:
  ✔ GNINA   (311/311 pairs, 100% coverage)
  ✔ Smina   (312/312 pairs, 100% coverage)
  ✔ Vina    (312/312 pairs, 100% coverage)
  (Use space to toggle, Enter to confirm)

> After selection (e.g., GNINA + Vina only):
  Routing to: multi-engine consensus [2 engines]
  Excluded: Smina (user_excluded)
```

```
> After selection (e.g., GNINA only):
  Routing to: single-engine [GNINA solo profile — CNN-ranked]
  Excluded: Smina, Vina (user_excluded)
```

---

## Acceptance Criteria

1. `--engines gnina` on a 3-engine project routes to GNINA solo profile; output tables sorted by `cnn_affinity`.
2. `--engines gnina,vina` on a 3-engine project runs 2-engine consensus; Smina scores absent from all output tables; `engines_in_scope: ["gnina","vina"]` in session manifest.
3. Interactive picker shows coverage % per engine before selection.
4. Deselecting all engines shows an error and re-presents the picker.
5. `engine_support_count` in top-pose atlas equals the number of scoped engines, not total detected engines.
6. `START_HERE.md` lists scoped and excluded engines for the session.
7. Two sessions in the same project with different scopes both complete without error; `START_HERE.md` notes the scope difference.
8. `--engine-preset cnn_only` on a project with preset defined routes to GNINA solo.
9. Excluded engines remain in detection report as `excluded_this_session`; their output files are untouched.
10. Session metadata `exclusion_reason` field populated for all excluded engines.

---

## Implementation Phases

### Phase 1: Scope Parsing
Add `--engines` and `--engine-preset` CLI flags; parse and validate against detected engine set (FR-001, FR-014, FR-015).

### Phase 2: Interactive Picker
Add engine picker prompt with coverage display in interactive mode (FR-002, FR-003, FR-004, FR-016).

### Phase 3: Scope-to-Profile Routing
Implement scope size → solo or filtered multi-engine routing; activate correct profile (FR-006, FR-007, FR-008).

### Phase 4: Session Metadata + Exclusion Records
Write `engines_in_scope`, `excluded_engines`, `exclusion_reason` into session manifest and all output artifacts (FR-009, FR-010, FR-011).

### Phase 5: Output Document Updates
Update `START_HERE.md` generator to include scope summary and cross-session scope discrepancy notice (FR-012, FR-013).

### Phase 6: Tests
Smoke tests for all acceptance criteria.
