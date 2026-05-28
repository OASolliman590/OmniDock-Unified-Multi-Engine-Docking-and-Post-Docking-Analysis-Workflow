# Tasks: Engine Scope Filter

**Input**: `/specs/023-engine-scope-filter/spec.md`
**Scope**: Allow selective engine analysis in multi-engine projects via `--engines` flag, interactive picker, and scope presets.

## Format: `[ID] [P?] Description`
- **[P]**: parallelizable

## Phase 1: Scope Parsing

- [x] T001 Add `--engines <comma-list>` CLI flag to unified pipeline and simplified CLI entry points.
- [x] T002 Add `--engine-preset <name>` CLI flag; load preset from `project_manifest.json["engine_presets"]`.
- [x] T003 Validate scoped engines against detected engine set; error if named engine is `absent`.
- [x] T004 Allow `--engines` to specify a single engine and route to that engine's solo profile (no multi-engine).

## Phase 2: Interactive Picker

- [x] T005 Add engine picker checkbox prompt in interactive mode when N ≥ 2 detected engines.
- [x] T006 Display coverage % per engine (from detection report) alongside each engine name in picker.
- [x] T007 Block confirmation when zero engines selected; re-present picker with error message.
- [x] T008 Add "load preset" option to picker; list available presets from `project_manifest.json["engine_presets"]`.

## Phase 3: Scope-to-Profile Routing

- [x] T009 When scoped N=1, activate named engine's solo profile (spec 022 solo profile contracts).
- [x] T010 When scoped N≥2, filter score inputs to scoped engines only before consensus computation.
- [x] T011 Set `engine_support_count` in top-pose atlas to scoped engine count; add `scoped_engine_count: N of M` annotation.
- [x] T012 Ensure excluded engine scores do not appear in any output table (rank correlation, consensus, hit classes).

## Phase 4: Session Metadata + Exclusion Records

- [x] T013 Write `engines_in_scope` and `excluded_engines` lists to session manifest on every run.
- [x] T014 [P] Record `exclusion_reason` per excluded engine (`user_excluded | low_coverage | parse_failed | cli_override`).
- [x] T015 Mark excluded engines as `excluded_this_session` in detection report (not `absent`).
- [x] T016 Preserve excluded engine output files; do not delete or overwrite.

## Phase 5: Output Document Updates

- [x] T017 Update `START_HERE.md` generator: add scope summary section (engines in scope, excluded, exclusion reasons).
- [x] T018 Detect cross-session scope differences in same project; add discrepancy notice to `START_HERE.md`.

## Phase 6: Tests

- [x] T019 Smoke test: `--engines gnina` on 3-engine project → GNINA solo profile, `cnn_affinity` primary sort.
- [x] T020 Smoke test: `--engines gnina,vina` → 2-engine consensus, Smina absent from all outputs.
- [x] T021 Smoke test: interactive picker with 0 engines → error + re-prompt.
- [x] T022 Smoke test: `engine_support_count` equals scoped engine count, not total detected.
- [x] T023 Smoke test: `--engine-preset cnn_only` routes to GNINA solo.
- [x] T024 Smoke test: two sessions with different scopes complete; `START_HERE.md` notes discrepancy.
- [x] T025 Smoke test: excluded engine files untouched after scoped run.
