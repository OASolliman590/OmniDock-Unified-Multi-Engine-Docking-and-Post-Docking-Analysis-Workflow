# Feature Specification: Engine Detection and Solo Routing

**Feature Branch**: `022-engine-detection-and-solo-routing`
**Created**: 2026-04-01
**Status**: Draft
**Input**: User observation that post-docking analysis requires manual engine selection; inspection of `3-Docking/` confirmed three distinct output layouts (GNINA SDF, Smina PDBQT/scores, Vina PDBQT/scores) and authoritative engine declarations in `project_manifest.json`.

---

## Summary

Post-docking analysis currently asks users to manually select the docking engine, enter scores paths, and pick between single-engine and multi-engine workflows. This is error-prone and inconsistent — a project already contains all the information needed to determine which engines ran, whether their outputs are valid, and which analysis path to take.

This spec adds an **Engine Detection Engine** that:
1. Reads authoritative engine declarations from `project_manifest.json`.
2. Falls back to filesystem-based detection when no manifest exists.
3. Validates each declared engine's output (file count, parse-ability, coverage).
4. Routes automatically: single-engine analysis when exactly one engine has valid output, multi-engine analysis when N ≥ 2.
5. Defines engine-specific **Solo Profiles** for GNINA, Smina, and Vina that configure score extraction, ranking, confidence signals, and RMSD handling appropriately for each engine when run alone.

---

## Observed Output Layouts (from `3-Docking/`)

| Engine | Pose Format | Score Location | Log Format |
|--------|------------|----------------|-----------|
| GNINA | SDF (`gnina_out/*.sdf`) | Parsed from log (5-col table) | `mode vina_affinity intramol cnn_score cnn_affinity` |
| Smina | PDBQT (`smina_out/poses/*.pdbqt`) | `smina_out/scores/normalized_scores.csv` | `mode affinity rmsd_lb rmsd_ub` (rmsd=-1.0) |
| Vina | PDBQT (`vina_out/poses/*.pdbqt`) | `vina_out/scores/normalized_scores.csv` | `mode affinity rmsd_lb rmsd_ub` (real RMSD) |

**GNINA log distinguishing markers**: `CNNscore`, `cnn_affinity`, 5-column score table, `WARNING: No GPU detected. CNN scoring will be slow.`

**Smina log distinguishing markers**: `smina is based off AutoDock Vina`, weights/terms block (`gauss`, `repulsion`, `hydrophobic`, `non_dir_h_bond`), `rmsd_lb/rmsd_ub = -1.0` in scores CSV.

**Vina log distinguishing markers**: progress bar `0%   10   20   30...`, 4-column table with real RMSD, no CNN lines.

**Manifest source**: `4-Docking/project_manifest.json` → `"engines": ["gnina", "vina", "smina"]`

---

## Goals

1. Zero-click engine detection — no user prompt required in normal use.
2. Correct routing to single vs. multi-engine analysis without user intervention.
3. Engine-specific solo profiles that extract the right scores and apply appropriate confidence signals per engine.
4. A machine-readable detection report for debugging and audit.
5. CLI and interactive `--engine` flag override for expert users.

## Non-Goals

- Adding new docking engines (AutoDock4, rDock, Glide handled in separate specs).
- Changing docking execution (detection is read-only inspection of existing outputs).
- Changing multi-engine consensus logic (covered in spec 017/020).

---

## User Stories

### US1: Zero-Config Start
As a scientist, I want to point the pipeline at my project folder and have it detect which engines ran and route the analysis automatically so I do not need to configure anything manually.

### US2: Correct Solo Profile
As a GNINA user, I want CNN affinity used as my primary ranking score (not just Vina affinity) when running solo so I get the full benefit of GNINA's CNN scoring.

### US3: Valid Output Gating
As an operator, I want engines with zero or corrupt output excluded from detection so a failed engine does not silently corrupt the analysis.

### US4: Detection Transparency
As a reviewer, I want a detection report that shows which engines were found, which were excluded, why, and which path was taken so the analysis is fully reproducible.

### US5: Override
As a power user, I want to override detection with `--engine gnina` to force single-engine analysis even when multiple engines have valid output.

---

## Functional Requirements

### Detection Layer

- **FR-001**: The detection engine MUST first read `project_manifest.json` `"engines"` field as the authoritative engine list.
- **FR-002**: When no `project_manifest.json` exists, the detection engine MUST fall back to filesystem probe: scan for `{engine}_out/` directories (gnina, smina, vina) and `4-Docking/logs/` with engine-identifiable log markers.
- **FR-003**: Filesystem detection MUST use log content markers (not just directory names) to distinguish GNINA from Smina/Vina, since Smina and Vina produce identical directory structures and log table formats differ only in the CNN lines and scoring weights block.
- **FR-004**: Each detected engine MUST pass a validity gate before being included in routing:
  - At least one valid pose file found (`*.sdf` for GNINA, `*.pdbqt` for Smina/Vina)
  - At least one valid score row parseable from scores CSV or log
  - Output coverage ≥ `min_coverage_pct` (default 30%) of pairlist pairs; engines below threshold are flagged as `partial` but still included if ≥ 1 valid file
  - Engines with 0 valid files are marked `absent` and excluded from routing
- **FR-005**: Detection MUST produce a `engine_detection_report.json` artifact in `4-Working/metadata/` containing: detected engines, validation status per engine, coverage counts, routing decision, and detection method (manifest vs. filesystem).
- **FR-006**: The routing decision MUST be: `single_engine` when exactly one engine passes the validity gate; `multi_engine` when N ≥ 2 pass; `no_valid_engines` when none pass (pipeline aborts with actionable error).
- **FR-007**: `--engine <name>` CLI flag MUST override detection and force single-engine mode for the named engine regardless of what other engines have valid output. The detection report MUST record `detection_override: true` and the reason.
- **FR-008**: Detection result MUST be cached in `project_manifest.json` under `"detected_engines"` so repeated runs do not re-probe unless `--redetect` is passed.

### GNINA Solo Profile

- **FR-010**: When GNINA is the sole valid engine, the pipeline MUST use `cnn_affinity` as the primary ranking score.
- **FR-011**: `vina_affinity` MUST be retained as a secondary score and included in all output tables.
- **FR-012**: `cnn_score` (range 0–1, probability of activity) MUST be used as a per-pose confidence signal: `cnn_score ≥ 0.5` → `high_cnn_confidence`; `0.3–0.5` → `moderate`; `< 0.3` → `low`.
- **FR-013**: CNN confidence MUST be included in hit classification output rows as `cnn_confidence` field.
- **FR-014**: Pose files are SDF format. The pipeline MUST use the SDF multi-conformer parser (fixed by spec 020 SDF parser fix) for RMSD and complex generation.
- **FR-015**: GNINA log parser MUST use the 5-column regex: `mode vina_affinity intramol cnn_score cnn_affinity`. The `intramol` column is discarded (ligand internal energy, not a binding score).
- **FR-016**: Hit classification thresholds for GNINA solo MUST use `cnn_affinity` for percentile ranking, not `vina_affinity`.
- **FR-017**: Ligand efficiency MUST be computed as `|cnn_affinity| / heavy_atom_count` for GNINA solo.

### Smina Solo Profile

- **FR-020**: When Smina is the sole valid engine, the pipeline MUST use `vina_affinity` (Smina's primary output) as the ranking score.
- **FR-021**: Smina's custom scoring function weights (the `gauss`, `repulsion`, `hydrophobic`, `non_dir_h_bond`, `num_tors_div` terms logged in the header) MUST be extracted and stored in the detection report for provenance. If weights differ across logs (custom run), a warning MUST be emitted.
- **FR-022**: `rmsd_lb` and `rmsd_ub` values from Smina output are always `-1.0` (Smina does not compute inter-pose RMSD). The pipeline MUST treat these as `not_available` rather than real RMSD values, and MUST suppress any RMSD-based clustering that depends on reported RMSD.
- **FR-023**: Pose files are PDBQT multi-model format. Complex generation MUST use the PDBQT pose parser.
- **FR-024**: Smina solo outputs MUST include a `smina_scoring_function` provenance field in the run manifest capturing which scoring terms were active.

### Vina Solo Profile

- **FR-030**: When Vina is the sole valid engine, the pipeline MUST use `vina_affinity` as the ranking score.
- **FR-031**: Vina's built-in inter-pose RMSD (`rmsd_lb`, `rmsd_ub`) MUST be used for pose diversity reporting. `rmsd_lb` (lower bound) is the preferred metric.
- **FR-032**: Pose files are PDBQT multi-model format. Complex generation MUST use the PDBQT pose parser.
- **FR-033**: Vina solo outputs MUST flag when the best pose has `rmsd_lb = 0.0` for all poses (all modes identical geometry) as a potential docking collapse warning.
- **FR-034**: Because Vina uses a progress-bar log format (`0%   10   20   30...`), the score parser MUST strip the progress bar lines before parsing the score table.

### Common Solo Behaviour

- **FR-040**: All solo profiles MUST run the same 16-step unified pipeline contract (spec 020 FR-000).
- **FR-041**: Solo profiles MUST disable all multi-engine steps (cross-engine rank correlation, DockBox geometric consensus, engine agreement fraction) and replace them with single-engine equivalents: within-engine pose spread, cnn_score confidence (GNINA only), and Vina RMSD pose diversity (Vina only).
- **FR-042**: The `engine_support_count` field in top-pose atlas output MUST be set to 1 for all solo runs and annotated with `single_engine_mode: true`.
- **FR-043**: QC downgrade logic MUST NOT fire the `warn_low_engine_agreement` rule in solo mode (spec 020 FR-013 fix).

---

## Detection Algorithm

```
1. Look for project_manifest.json → read "engines" list
   - If found: engine_source = "manifest"
   - If not: probe filesystem for {gnina,smina,vina}_out/ directories
              engine_source = "filesystem"

2. For each candidate engine:
   a. Count pose files in expected location
   b. Try to parse one score record from scores CSV or log
   c. Compute coverage = valid_pairs / total_pairlist_pairs
   d. Assign status: "valid" | "partial" | "absent"

3. valid_engines = [e for e in candidates if status != "absent"]

4. Route:
   len(valid_engines) == 0 → abort: no_valid_engines
   len(valid_engines) == 1 → single_engine(valid_engines[0])
   len(valid_engines) >= 2 → multi_engine(valid_engines)

5. Write engine_detection_report.json
6. Cache in project_manifest.json["detected_engines"]
```

---

## Detection Report Schema

```json
{
  "detection_method": "manifest | filesystem",
  "detection_override": false,
  "override_reason": null,
  "routing_decision": "single_engine | multi_engine | no_valid_engines",
  "routed_to_engine": "gnina | smina | vina | null",
  "engines": {
    "gnina": {
      "status": "valid | partial | absent",
      "pose_files_found": 311,
      "pose_format": "sdf",
      "valid_score_rows": 3110,
      "coverage_pct": 100.0,
      "score_columns": ["vina_affinity", "cnn_score", "cnn_affinity"],
      "log_marker_matched": "cnn_affinity"
    },
    "smina": {
      "status": "valid",
      "pose_files_found": 312,
      "pose_format": "pdbqt",
      "valid_score_rows": 3120,
      "coverage_pct": 100.0,
      "score_columns": ["vina_affinity"],
      "smina_scoring_weights": {...},
      "rmsd_available": false
    },
    "vina": {
      "status": "valid",
      "pose_files_found": 312,
      "pose_format": "pdbqt",
      "valid_score_rows": 3120,
      "coverage_pct": 100.0,
      "score_columns": ["vina_affinity"],
      "rmsd_available": true
    }
  }
}
```

---

## Engine Output Layout Reference

```
4-Docking/
├── project_manifest.json          ← "engines": ["gnina","vina","smina"]
├── pairlist.csv
├── logs/                          ← GNINA logs ({tag}.log)
├── gnina_out/                     ← GNINA SDF pose files ({tag}.sdf)
├── smina_out/
│   ├── logs/                      ← Smina logs
│   ├── poses/                     ← Smina PDBQT poses
│   └── scores/normalized_scores.csv
└── vina_out/
    ├── logs/                      ← Vina logs
    ├── poses/                     ← Vina PDBQT poses
    └── scores/normalized_scores.csv
```

---

## Acceptance Criteria

1. On a 3-engine project (`gnina`, `smina`, `vina` all valid), pipeline routes to `multi_engine` without user prompt.
2. On a project where only `gnina_out/` has files, pipeline routes to `single_engine(gnina)` without prompt.
3. `engine_detection_report.json` is written on every run with correct `routing_decision`.
4. GNINA solo run uses `cnn_affinity` as primary sort column in all output tables, not `vina_affinity`.
5. Smina solo run does not attempt RMSD clustering using `-1.0` values.
6. Vina solo run uses `rmsd_lb` for pose diversity reporting.
7. `--engine gnina` override forces `single_engine(gnina)` even when smina/vina also have valid output; detection report records `detection_override: true`.
8. An engine with 0 valid pose files is excluded from routing and logged in the detection report as `absent`.
9. Solo runs do not emit `warn_low_engine_agreement`.
10. Detection result is cached in `project_manifest.json` and reused on subsequent runs without re-probing.

---

## Risks

- Projects that ran GNINA without CNN scoring (`--cnn_scoring=none`) produce 3-column logs, not 5-column. Detection must handle both log formats for GNINA.
- Custom Smina runs with non-default scoring terms may produce different weight blocks; the parser must be tolerant of unknown term names.
- A project could have partial outputs for all engines (e.g., an HPC run that timed out). The `partial` status allows routing but must surface coverage warnings to the user.

---

## Implementation Phases

### Phase 1: Detection Core
Implement `post_docking_analysis/engine_detector.py` with manifest-first + filesystem-fallback logic and report generation (FR-001 to FR-008).

### Phase 2: Routing Integration
Wire detection result into the unified pipeline entry point. Replace manual engine selection prompt with detection result; keep manual override prompt only when detection is ambiguous or `--engine` flag is passed (FR-006, FR-007).

### Phase 3: GNINA Solo Profile
Implement GNINA-specific score extraction (5-column parser), CNN confidence field, `cnn_affinity`-first ranking, SDF pose handling (FR-010 to FR-017).

### Phase 4: Smina Solo Profile
Implement Smina-specific score extraction, scoring function provenance capture, RMSD suppression (FR-020 to FR-024).

### Phase 5: Vina Solo Profile
Implement Vina-specific score extraction, progress-bar log stripping, `rmsd_lb` pose diversity, collapse warning (FR-030 to FR-034).

### Phase 6: Common Solo Behaviour + Tests
Enforce solo-mode contracts (FR-040 to FR-043), disable multi-engine steps in solo, smoke tests for all acceptance criteria.
