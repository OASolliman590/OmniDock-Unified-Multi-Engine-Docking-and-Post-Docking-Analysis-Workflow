# Tasks: Engine Detection and Solo Routing

**Input**: `/specs/022-engine-detection-and-solo-routing/spec.md`
**Scope**: Auto-detect docking engines from project outputs, route to single or multi-engine analysis, apply engine-specific solo profiles.

## Format: `[ID] [P?] Description`
- **[P]**: parallelizable

## Phase 1: Detection Core

- [x] T001 Create `post_docking_analysis/engine_detector.py` with `detect_engines(project_dir)` entry point.
- [x] T002 Implement manifest-first detection: read `project_manifest.json` `"engines"` field.
- [x] T003 Implement filesystem fallback: probe `{engine}_out/` directories + log marker matching for GNINA vs Smina/Vina disambiguation.
- [x] T004 Implement validity gate per engine: count pose files, parse one score row, compute coverage %.
- [x] T005 Assign engine status: `valid | partial | absent`; emit coverage warnings for `partial`.
- [x] T006 Write `engine_detection_report.json` to `4-Working/metadata/` on every run.
- [x] T007 Cache detection result in `project_manifest.json["detected_engines"]`; skip re-probe unless `--redetect` passed.

## Phase 2: Routing Integration

- [x] T008 Wire `detect_engines()` into unified pipeline entry point (replace manual engine selection prompt).
- [x] T009 Implement routing: 0 valid → abort with error; 1 valid → `single_engine`; N≥2 → `multi_engine`.
- [x] T010 Implement `--engine <name>` override: force `single_engine` for named engine, record `detection_override: true` in report.
- [x] T011 Keep manual engine selection prompt only as fallback when `no_valid_engines` and user passes `--interactive-fallback`.

## Phase 3: GNINA Solo Profile

- [x] T012 [P] Implement 5-column GNINA log parser (`mode vina_affinity intramol cnn_score cnn_affinity`); handle 3-column fallback when `--cnn_scoring=none`.
- [x] T013 [P] Set `cnn_affinity` as primary ranking score in GNINA solo mode.
- [x] T014 [P] Add `cnn_confidence` field from `cnn_score` thresholds (≥0.5=high, 0.3–0.5=moderate, <0.3=low).
- [x] T015 Compute ligand efficiency as `|cnn_affinity| / heavy_atom_count` in GNINA solo.
- [x] T016 Confirm SDF pose parser (spec 020 fix) is used for complex generation in GNINA solo.

## Phase 4: Smina Solo Profile

- [x] T017 [P] Implement Smina score extractor from `smina_out/scores/normalized_scores.csv`.
- [x] T018 [P] Extract and store Smina scoring function weights from log header into detection report.
- [x] T019 Suppress RMSD clustering steps that depend on `rmsd_lb/rmsd_ub` in Smina solo (mark as `not_available`).
- [x] T020 Add `smina_scoring_function` provenance field to run manifest in Smina solo.
- [x] T021 Emit warning when scoring weights differ across logs in same project (inconsistent custom scoring).

## Phase 5: Vina Solo Profile

- [x] T022 [P] Implement Vina score extractor from `vina_out/scores/normalized_scores.csv`; strip progress-bar lines from log before parsing.
- [x] T023 [P] Use `rmsd_lb` from Vina output for pose diversity reporting.
- [x] T024 Detect and warn on docking collapse: best pose `rmsd_lb = 0.0` across all modes.

## Phase 6: Common Solo Behaviour + Tests

- [x] T025 Disable multi-engine steps in all solo profiles (cross-engine correlation, engine agreement fraction, DockBox consensus).
- [x] T026 Set `engine_support_count=1` and `single_engine_mode: true` in top-pose atlas for all solo runs.
- [x] T027 Confirm `warn_low_engine_agreement` does not fire in solo mode.
- [x] T028 Smoke test: 3-engine project → routes to `multi_engine` without prompt.
- [x] T029 Smoke test: GNINA-only project → routes to `single_engine(gnina)`, output tables sorted by `cnn_affinity`.
- [x] T030 Smoke test: `--engine gnina` override on 3-engine project → `detection_override: true` in report.
- [x] T031 Smoke test: project with 0 valid engine outputs → `no_valid_engines` error with actionable message.
- [x] T032 Smoke test: Smina solo run → no RMSD clustering from `-1.0` values.
- [x] T033 Smoke test: Vina solo run → `rmsd_lb` in pose diversity table, progress-bar not in score table.
