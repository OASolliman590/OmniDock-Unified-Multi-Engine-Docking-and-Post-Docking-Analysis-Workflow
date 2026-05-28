# Feature Specification: Omni-DockForge: End-to-End Docking, Consensus by Design.

**Feature Branch**: `018-dockforge-omnidock-platform`  
**Created**: 2026-03-30  
**Status**: Draft  
**Input**: User-driven platform refactor request for a spec-first DockForge upgrade covering interactive UX, engine-aware preparation/execution, post-docking promotion, consensus scoring, QC, reproducibility, and extensibility.

## Product Vision

DockForge is the next-generation evolution of the current docking workflow: an omni-docking platform that treats preparation, multi-engine docking, and post-docking analysis as first-class, auditable, and modular workflows.

DockForge must shift from a single-run utility into a reusable project platform with:
- consistent project topology and provenance,
- non-blocking interactive orchestration for long-running tasks,
- engine-aware preparation and execution,
- scientifically defensible cross-engine comparison and consensus,
- optional local and environment-aware runtime backends.

## Personas and Usage Modes

1. **Interactive Discovery User (bench/computational scientist)**
- Uses guided interactive panels.
- Needs visible progress, back navigation, and the ability to keep working while jobs run.

2. **Pipeline Operator (computational chemist / HPC user)**
- Runs multi-engine projects at scale.
- Needs environment selection, strict input validation, predictable folder contracts, and rerunnable runs.

3. **Scientific Reviewer (PI / collaborator / manuscript pipeline)**
- Consumes ranked outputs and supporting evidence.
- Needs explainable consensus, QC traces, confidence metrics, and report-ready artifacts.

4. **Method Developer (platform maintainer)**
- Adds or modifies engines, scoring methods, and analysis plugins.
- Needs clean module boundaries and extension points.

## Workflow Map

1. **Project Entry**
- Create/resume project.
- Choose stage or analysis-only entry.

2. **Input and Preparation**
- Receptor cleanup/QC gating.
- Engine-aware ligand preparation graph (Open Babel, Meeko, AutoDockTools combinations).

3. **Grid and Docking Setup**
- Validate box/search-space.
- Select engines and parameter presets.
- Select execution environment (local CPU/GPU/conda/container/remote adapter).

4. **Docking Execution**
- Run as background jobs with live state updates.
- Keep UI navigable while jobs run.

5. **Post-Docking (first-class)**
- Analyze existing projects without rerunning docking.
- Run modular analysis scopes (comparison-only, rescoring-only, QC-only, report-only).

6. **Consensus, QC, and Classification**
- Build allscore from normalized engine outputs.
- Compute target-aware strong/moderate/weak classes.
- Integrate ADMET and optional external biology correlations.

7. **Reporting and Provenance**
- Emit canonical 0-7 folder hierarchy + `.meta` manifests + session snapshots.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Non-Blocking Interactive Workflow Control (Priority: P1)

As an interactive user, I can start long-running steps, navigate to other panels, monitor step-level progress, and return back from any panel without losing state.

**Why this priority**: Existing blocking behavior and missing back navigation create major workflow friction and user errors.

**Independent Test**: Start a long docking or post-docking task, navigate to two other panels, inspect progress and state timeline, then return and verify the job is still running/completed correctly.

**Acceptance Scenarios**:

1. **Given** a step is running, **When** the user navigates to another panel, **Then** the step continues and shows status `in_progress` with progress updates.
2. **Given** a completed step with output validation, **When** workflow timeline is rendered, **Then** that step is marked `validated` and visually highlighted.
3. **Given** any interactive panel, **When** user chooses `Back`, **Then** the previous panel opens without losing selected project context.

---

### User Story 2 - Engine-Aware Ligand Preparation and Validation (Priority: P1)

As a pipeline operator, I can choose a preparation mode that correctly orchestrates shared and engine-specific ligand preparation steps and rejects invalid combinations before execution.

**Why this priority**: Ligand preparation correctness is foundational for all downstream engines and currently causes failures/confusion.

**Independent Test**: Run each preparation mode with a sample ligand set and verify expected artifacts per engine, including Open Babel enrichment steps and output format validation.

**Acceptance Scenarios**:

1. **Given** `Open Babel -> Meeko -> AutoDock-compatible output`, **When** preparation runs, **Then** Open Babel enrichment executes first, Meeko conversion executes second, and AutoDock4-compatible artifacts are produced with validation report.
2. **Given** `AutoDockTools only`, **When** required script path is missing, **Then** workflow fails fast with actionable diagnostics and does not proceed to docking.
3. **Given** preparation phase menus, **When** user opens options, **Then** PLIP is not present in preparation choices.

---

### User Story 3 - Environment-Aware Docking Execution with Parameter Control (Priority: P1)

As a computational user, I can run docking using a selected execution backend and parameter schema (basic or advanced), with engine-specific overrides validated before launch.

**Why this priority**: Docking runtime and reproducibility depend on explicit environment and parameter control.

**Independent Test**: Configure and run one multi-engine job in local mode and one in environment-backed mode using the same project; verify required file checks, parameter validation, and run manifests.

**Acceptance Scenarios**:

1. **Given** missing required inputs for an engine, **When** user starts docking, **Then** preflight fails and lists missing files by engine.
2. **Given** GNINA with GPU-aware options enabled, **When** parameters are validated, **Then** unsupported combinations are rejected with clear guidance.
3. **Given** advanced mode off, **When** user runs docking, **Then** only safe defaults and preset fields are required.

---

### User Story 4 - Post-Docking as Project-Centric First-Class Entry (Priority: P1)

As an analysis-focused user, I can open a finished project and run analysis-only workflows (comparison, rescoring, QC, reporting) without re-running docking.

**Why this priority**: Post-docking must be central to scientific interpretation, not a hidden downstream add-on.

**Independent Test**: Open an existing completed project and run comparison-only and report-only modes, confirming new outputs are generated under canonical analysis directories.

**Acceptance Scenarios**:

1. **Given** a completed project, **When** user selects `Post-Docking Analysis` entry mode, **Then** analysis modules run without invoking docking steps.
2. **Given** `QC-only` mode, **When** executed, **Then** QC outputs are generated in `4-PostDocking/filtering` and `5-Analysis/structure_quality`.
3. **Given** analysis session restart, **When** reopened, **Then** session state and completed modules are visible.

---

### User Story 5 - Explainable AllScore Consensus, QC, and Biology Correlation (Priority: P2)

As a scientific reviewer, I can inspect normalized per-engine scores, consensus ranking logic, confidence metrics, hit classes, and optional external biology correlations.

**Why this priority**: Scientific credibility requires traceable, target-aware interpretation rather than opaque raw-score ranking.

**Independent Test**: Run consensus generation on a multi-engine project with reference ligand metadata and an external biology file; verify ranking, explainability, correlation outputs, and class labels.

**Acceptance Scenarios**:

1. **Given** scores from multiple engines, **When** allscore is generated, **Then** raw, normalized, and consensus tables are output with method metadata.
2. **Given** optional decoys/reference ligands, **When** benchmarking runs, **Then** EF1%/EF5%/ROC-AUC outputs are produced where data suffices.
3. **Given** external biology file mapping, **When** correlation pipeline runs, **Then** unresolved entities are reported and resolved entities are included in correlation tables.

---

### User Story 6 - Reproducible Checkpoint and Auditability (Priority: P2)

As a platform maintainer, I can freeze a run state and revise from that checkpoint with deterministic manifests and environment snapshots.

**Why this priority**: Reproducibility and recoverability are required for long projects and publication-quality evidence trails.

**Independent Test**: Complete a run, create a checkpoint, revise settings, and verify both checkpoint lineage and updated outputs are auditable in `.meta` and `sessions`.

**Acceptance Scenarios**:

1. **Given** a completed or partial run, **When** user triggers `Checkpoint & Revise`, **Then** a frozen checkpoint manifest is saved and editable revision state is created.
2. **Given** a rerun from checkpoint, **When** run starts, **Then** provenance hash changes only if inputs/configs change.

## Edge Cases

- User selects conflicting preparation path (for example Meeko-only for AutoDock4 without compatible prerequisites).
- User starts background task, exits panel, then resumes after application restart.
- Partial engine failures (one engine success, others fail) within multi-engine run.
- Missing optional tools (Open Babel, Meeko, AutoDockTools, PandaMap, ProLIF, OnionNet2).
- Positive affinity or missing pose ligands in some engines but not all engines.
- Mixed receptor assemblies (protein + DNA/multichain) requiring chain-retention policy.
- External biology file IDs that do not map to ligand/protein naming aliases.
- High-scale projects where CSV-only processing becomes memory-bound.

## Requirements *(mandatory)*

### Functional Requirements

### Interactive UX and Workflow State

- **FR-001**: System MUST rename platform branding to `Omni-DockForge: End-to-End Docking, Consensus by Design.` in interactive headers and key docs.
- **FR-002**: System MUST maintain a workflow state model with statuses: `not_started`, `in_progress`, `completed`, `validated`, `failed`, `skipped`, `needs_review`.
- **FR-003**: Interactive UI MUST display workflow position timeline beneath current panel choices.
- **FR-004**: Completed and validated steps MUST be visually distinct (green for `validated` status).
- **FR-005**: Long-running steps MUST run as background jobs and MUST NOT block navigation.
- **FR-006**: System MUST provide live or near-live progress for active background jobs.
- **FR-007**: Every interactive panel MUST include a universal `Back` action.
- **FR-008**: System MUST allow user to inspect background jobs from any panel.
- **FR-009**: System MUST persist workflow state to session storage for recovery after restart.

### Checkpoint & Revise

- **FR-010**: System MUST replace ambiguous `maturation` label with `Checkpoint & Revise` across UX and docs.
- **FR-011**: `Checkpoint & Revise` MUST freeze current state/config into immutable checkpoint metadata.
- **FR-012**: System MUST support creating a revision branch/state from a checkpoint without mutating frozen artifacts.

### Preparation Pipeline Rules

- **FR-013**: PLIP MUST NOT appear in preparation-stage prompts or execution.
- **FR-014**: Open Babel enrichment stage MUST support hydrogen addition, pH-based protonation, 3D conformer generation, energy minimization, and partial charge calculation where available.
- **FR-015**: Default ligand protonation pH MUST be 7.4 unless user overrides.
- **FR-016**: Preparation choices MUST include: `Open Babel only`, `Meeko only`, `AutoDockTools only`, `Open Babel -> Meeko`, `Open Babel -> Meeko -> AutoDock-compatible output`, `Open Babel -> AutoDockTools`, `Engine-aware full preparation`.
- **FR-017**: System MUST validate preparation-mode compatibility against selected docking engines before execution.
- **FR-018**: Engine-aware full preparation MUST treat Vina/Smina/GNINA as a shared family with optional per-engine overrides.
- **FR-019**: Engine-aware full preparation MUST treat AutoDock4 as a distinct preparation branch.
- **FR-020**: Preparation pipeline MUST emit per-ligand validation reports with failure reasons.

### Docking Execution and Environment Abstraction

- **FR-021**: System MUST expose explicit `Run Docking` workflow with environment selection.
- **FR-022**: Environment model MUST support at least: local CPU, local GPU, conda environment, containerized environment, and remote/HPC adapter hooks.
- **FR-023**: System MUST allow user to define working folder and prerequisites folder for run execution.
- **FR-024**: System MUST run preflight validation of required docking inputs before launch.
- **FR-025**: System MUST support local docking execution without HPC dependency.
- **FR-026**: Job execution layer MUST support both synchronous (debug) and asynchronous (production) modes.

### Docking Parameters and Profiles

- **FR-027**: System MUST provide a common parameter schema for all engines.
- **FR-028**: System MUST provide engine-specific parameter extensions.
- **FR-029**: Common schema MUST include at least: exhaustiveness, number of modes, seed, and search-space/box parameters.
- **FR-030**: GNINA extension schema MUST include CNN/scoring and GPU-relevant controls.
- **FR-031**: System MUST provide `basic` and `advanced` parameter modes.
- **FR-032**: System MUST provide validated default presets for screening and exhaustive modes.

### Post-Docking First-Class Workflows

- **FR-033**: Post-docking MUST be a top-level project entry point independent of rerunning docking.
- **FR-034**: System MUST support modular analysis scopes: comparison-only, rescoring-only, QC-only, report-only, full analysis.
- **FR-035**: System MUST support reopening prior sessions with stage/module completion visibility.

### AllScore, Consensus, and Correlation

- **FR-036**: System MUST ingest raw scores per engine and per pose into a unified score model.
- **FR-037**: System MUST normalize scores before cross-engine comparison.
- **FR-038**: Normalization methods MUST be configurable and recorded in metadata.
- **FR-039**: System MUST generate per-protein and global consensus rankings.
- **FR-040**: System MUST generate ligand-centric and protein-centric ranking views.
- **FR-041**: System MUST compute cross-engine rank correlations per protein and globally.
- **FR-042**: System SHOULD compute bootstrap-based rank stability/confidence intervals where compute budget permits.
- **FR-043**: System MUST emit explainability artifacts describing consensus decisions.

### Hit Classification and QC

- **FR-044**: System MUST implement configurable target-aware strong/moderate/weak hit classes.
- **FR-045**: Hit classification MUST consider normalized score rank bands and QC status.
- **FR-046**: ADMET flagging MUST be separated from docking quality labels and reported independently.
- **FR-047**: System MUST support reference-aware comparisons (reference ligand, known actives, decoy/background distribution, or internal percentile fallback).
- **FR-048**: System MUST reject single universal raw-affinity thresholds across all proteins/engines by default.

### QC, Benchmarking, and Structure Handling

- **FR-049**: System MUST provide upstream receptor QC gating before docking.
- **FR-050**: System MUST provide pre-docking ligand QC/filtering (including problematic ligand flags and ADMET filters).
- **FR-051**: System SHOULD support decoy-based enrichment benchmarking (EF1%, EF5%, ROC-AUC) when reference data exists.
- **FR-052**: System MUST generate and validate protein+pose complex structures for downstream visualization/rescoring.
- **FR-053**: Complex generation MUST preserve chain/residue mapping and report unresolved merges.

### External Biology and Rescoring Plugins

- **FR-054**: System MUST allow user to supply external biological findings files.
- **FR-055**: Supported biology file inputs MUST include CSV and TSV at minimum; JSON support SHOULD be available.
- **FR-056**: System MUST provide entity mapping and unresolved-entity reporting for biology file integration.
- **FR-057**: System MUST produce correlation outputs between biological findings and docking/consensus outputs.
- **FR-058**: System MUST support pluggable rescoring modules.
- **FR-059**: OnionNet2 MUST be available as an optional rescoring plugin, not hardwired as the only option.

### Folder Topology and Reproducibility

- **FR-060**: Project outputs MUST follow numbered canonical structure:
  - `.meta`, `0-Input`, `1-Preparation`, `2-GridBoxes`, `3-Docking`, `4-PostDocking`, `5-Analysis`, `6-Visualizations`, `7-Reports`, `sessions`.
- **FR-061**: `.meta` MUST include frozen run config, manifest with input checksums/tool versions/timestamps, and environment lock snapshot.
- **FR-062**: Per-engine configs/logs/poses MUST be retained under engine-specific directories.
- **FR-063**: System MUST generate rerun manifests in machine-readable formats (CSV + JSON).
- **FR-064**: Sessionized interactive states MUST be durable and auditable.

### Data Backend

- **FR-065**: System MUST preserve CSV exports for interoperability.
- **FR-066**: System SHOULD add a structured query backend (SQLite) for cross-run and cross-engine analysis.
- **FR-067**: If SQLite is enabled, schema MUST track runs, proteins, ligands, poses, scores, QC flags, and correlations.

### Future-Ready Extension Points

- **FR-068**: Architecture MUST expose extension points for ensemble receptor docking and multi-conformation analysis.
- **FR-069**: Architecture MUST expose extension points for automatic rerun execution from rerun manifests.
- **FR-070**: Architecture MUST expose extension points for dashboard-driven project inspection.

### Key Entities *(include if feature involves data)*

- **WorkflowState**: stage-level and step-level statuses, timestamps, validation outcomes, and progress percentages.
- **BackgroundTask**: asynchronous unit of work with task type, owner stage, progress, logs, and cancellation policy.
- **CheckpointRevision**: immutable checkpoint metadata plus mutable revision lineage.
- **PreparationProfile**: selected ligand preparation graph, pH policy, quality flags, and outputs.
- **EngineRunProfile**: selected engines, environment backend, parameter presets, and engine-specific overrides.
- **ExecutionEnvironment**: runtime backend details (local/conda/container/remote), resource constraints, and command templates.
- **ScoreRecord**: per-pose/per-engine raw score, normalized score, source metadata.
- **ConsensusRecord**: per-ligand/per-target consensus score, method, explainability breakdown, confidence indicators.
- **HitClassRecord**: strong/moderate/weak class assignment with contributing evidence and ruleset ID.
- **QCGateRecord**: receptor and ligand QC checks, pass/fail status, metrics, and blockers.
- **ComplexStructureRecord**: pose-to-receptor merge outputs, chain mapping, validation status.
- **BiologyAnnotationRecord**: external biology file mappings and per-entity correlation outputs.
- **RunManifest**: provenance payload with checksums, tool versions, environment snapshot, and run hash.

## Engine Orchestration Rules

### Preparation Graph Policies

| Mode | Shared Steps | Engine-Specific Steps | Required Outputs |
|---|---|---|---|
| Open Babel only | hydrogen/protonation/3D/minimize/charges | none | normalized ligand structures + optional PDBQT |
| Meeko only | none | Meeko conversion/validation | PDBQT (Vina-family compatible) |
| AutoDockTools only | optional format normalization | ADT prepare_ligand4 path | AD4-compatible PDBQT |
| Open Babel -> Meeko | Open Babel enrichment | Meeko convert | Vina/Smina/GNINA PDBQT |
| Open Babel -> AutoDockTools | Open Babel enrichment | ADT convert | AD4-compatible PDBQT |
| Open Babel -> Meeko -> AutoDock-compatible output | Open Babel enrichment | Meeko intermediate + AD4-compatible final conversion/validation | dual-family outputs |
| Engine-aware full preparation | Open Babel enrichment baseline | Vina-family branch + AD4 branch | all selected-engine outputs |

### Output Compatibility Rules

- Vina/Smina/GNINA consume shared-family ligand artifacts unless engine-specific overrides are enabled.
- AutoDock4 consumes AD4 branch artifacts and is validated separately.
- Invalid mode/engine combinations MUST fail before job execution.

## UX Behavior Rules

- Timeline is always visible in interactive workflow panels and includes stage index, status, and last update time.
- Running steps show progress percentage + current substep name.
- `Back` action is mandatory in all nested panels.
- Background task monitor is reachable from all major menus.
- Long-running operations write incremental logs into session folder and UI status stream.

## Output and Project Topology Contract

```text
PROJECT/
├── .meta/
├── 0-Input/
├── 1-Preparation/
├── 2-GridBoxes/
├── 3-Docking/
├── 4-PostDocking/
├── 5-Analysis/
├── 6-Visualizations/
├── 7-Reports/
└── sessions/
```

This structure is mandatory for new projects and migration-target structure for legacy projects.

## Non-Functional Requirements

- Deterministic output paths and naming for reruns.
- Graceful degradation when optional tools are unavailable.
- Preflight-first error handling to reduce wasted runtime.
- Session/state durability against interactive restarts.
- Performance target: interactive status refresh <= 2 seconds for active tasks under normal local workloads.
- Compatibility target: Python 3.9+ and existing conda workflow.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: 100% of interactive panels provide a working `Back` navigation path.
- **SC-002**: Users can start a long-running step and continue navigating to at least two other panels without stopping that step.
- **SC-003**: Engine-aware preparation passes compatibility validation for all supported engine combinations and rejects invalid combinations preflight.
- **SC-004**: Post-docking can be launched directly on an existing project with no docking rerun.
- **SC-005**: Consensus outputs include raw + normalized + consensus tables and explainability metadata for every analyzed ligand-target pair.
- **SC-006**: Hit classification rules are configurable per target and avoid global raw-score thresholds by default.
- **SC-007**: `.meta/run_manifest.json` includes input checksums, tool versions, and environment lock snapshot for every run.
- **SC-008**: CSV exports remain available; SQLite backend (if enabled) reproduces the same ranking outputs as CSV pipeline for smoke datasets.

## Risk Areas

- Asynchronous task orchestration complexity may introduce race conditions in interactive state.
- Preparation graph expansion could increase support burden for toolchain combinations.
- Score normalization choices can bias consensus if defaults are not transparent.
- External biology entity mapping may have high unresolved rates due to naming inconsistencies.
- SQLite adoption may require migration and dual-write consistency controls.

## Unresolved Design Decisions (with recommendation)

1. **Default normalization method for allscore**
- Options: z-score per target/engine, robust rank normalization, quantile normalization.
- **Recommendation**: robust rank normalization default, with optional z-score for sensitivity analyses.

2. **SQLite adoption timing**
- Options: immediate mandatory backend vs optional backend first.
- **Recommendation**: optional in Milestone 6 with dual CSV parity checks, then evaluate default switch.

3. **Background task runtime strategy**
- Options: thread pool, process pool, external queue.
- **Recommendation**: process pool for heavy compute tasks + lightweight thread/event loop status monitor.

4. **Reference-aware hit classification default policy**
- Options: internal percentiles only vs reference/decoy-aware when available.
- **Recommendation**: target-aware tiered policy: reference/decoy-aware when data exists, percentile fallback otherwise.

5. **Checkpoint storage granularity**
- Options: full snapshot copies vs manifest-linked immutable references.
- **Recommendation**: manifest-linked immutable references to avoid storage explosion; copy only small config/state artifacts.

6. **Name for old maturation step replacement**
- Options: `Checkpoint & Revise`, `Freeze & Edit`.
- **Recommendation**: `Checkpoint & Revise` as final product term due to clarity and scientific audit semantics.
