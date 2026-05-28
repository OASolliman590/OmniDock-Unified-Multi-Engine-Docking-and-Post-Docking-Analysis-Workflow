# Feature Specification: Post-Docking Analysis Remediation & Contract Uniformity

**Feature Branch**: `029-post-docking-remediation`
**Created**: 2026-05-23
**Status**: Draft
**Input**: User description: "draft this into a spec it cycle to fix it" — referring to the interactive-pipeline alignment audit plus the still-open items from the original post-docking gap analysis.

## Context

The simplified post-docking pipeline was recently hardened (`simplified_cli.py:334,356,499`; `_legacy_simplified_pipeline_impl.py:156,193,324,403,471,1937,1945,2555,3567,3642,3769`) so that RMSD, visualizations, ProLIF and LigPlot are mandatory and emit `run_tracking/` artifacts. Two legacy entry paths were blocked (`cli.py:338`, `pipeline.py:1366`), and a `python main.py analyze clean` shortcut was added (`workflow/cli.py:107,623,1398`, `workflow/execution.py:1751,1782`).

The interactive surface (`workflow/interactive.py` + `workflow/execution.py`) and several non-GNINA paths still bypass that contract:

- Non-GNINA path forces `run_rmsd=False` (`workflow/execution.py:2482`).
- Per-stage menu items (`analyze.stage.rmsd`, `analyze.interactions.prolif`, `.ligplot`, `.pandamap`, `.poseview`, `analyze.visuals.py3dmol`, `.pymol`) run in isolation, defeating the joint ProLIF+LigPlot requirement (`workflow/interactive.py:73-79`, `workflow/execution.py:2354-2363`).
- `PostDockingAnalysisPipeline` class is still importable and reachable from `analyze.stage.structure_quality` / `analyze.visuals.pymol` (`workflow/execution.py:2366,2371`).
- `run_tracking/` outputs exist but are never surfaced in the interactive completion message (`workflow/interactive.py:2507-2553`).
- Stub modules (`post_docking_analysis/rmsd_analyzer.py`, `post_docking_analysis/affinity_analyzer.py`) still exist; `plugin_manager.py` is marked for removal but still imported by `pipeline.py`; PandaMap is marked broken but still wired.
- Vina/Smina/AutoDock4 have no HPC adapter parallel to `gnina_hpc_adapter.py` and no per-engine parser/hard-fail tests.
- Config surface is fragmented (`config.py`, `config_manager.py`, hardcoded simplified_cli defaults).
- Documentation (`SIMPLIFIED_IMPLEMENTATION_STATUS.md`, `REDESIGN_PROPOSAL.md`, `INTERACTIVE_CLI_DECISION_TREE.md`, `README.md`) is now out of sync with the new contract.

**Intended outcome:** every user-facing surface (simplified CLI, interactive menu, `analyze clean` shortcut, multi-engine path) enforces the same hardened contract, surfaces the same observability, and is backed by the same test suite, for every supported engine — with the dead code, duplicate stubs, and stale docs removed.

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Interactive workflow honors the hardened contract for every engine (Priority: P1)

A researcher launches `python main.py` and walks through the questionary-based interactive workflow on a Vina (or Smina, or AutoDock4, or GNINA) project. Every analysis target they pick — favorite-engine, comparative, or a specific stage — runs with RMSD enforced, ProLIF+LigPlot jointly required, and hard-fails surfaced to the user.

**Why this priority**: This is a correctness regression today. Users believe the new contract protects them; on the interactive path it silently doesn't, producing reports that look complete but are missing mandatory science stages.

**Independent Test**: Run `python main.py` interactively against a fixture project for each of {GNINA, Vina, Smina, AutoDock4}, pick `analyze.comparative`, and assert that `run_tracking/run_manifest.json` shows `enforced == requested` for `run_rmsd`, `run_visualizations`, ProLIF, and LigPlot, regardless of engine and regardless of menu choice.

**Acceptance Scenarios**:
1. **Given** a Vina project, **When** the user selects `analyze.favorite_engine`, **Then** the simplified pipeline runs with `run_rmsd=True` and the manifest records the enforcement.
2. **Given** any engine, **When** the user picks `analyze.interactions.prolif` alone, **Then** the system either also runs LigPlot (joint contract honored) or refuses with a clear message pointing to `analyze.interactions.clean`.
3. **Given** any engine, **When** the user picks `analyze.stage.structure_quality` or `analyze.visuals.pymol`, **Then** the run does not invoke the legacy `PostDockingAnalysisPipeline` stub-bearing class; either the targets are migrated to the hardened path or removed from the menu with migration guidance.

---

### User Story 2 — Vina / Smina / AutoDock4 reach feature parity with GNINA (Priority: P1)

A researcher running any supported engine receives the same post-docking artifacts, the same hard-fail behavior, and the same HPC deployment ergonomics as a GNINA user.

**Why this priority**: The product claims multi-engine support, but `gnina_hpc_adapter.py`, the engine-specific log parsers in `docking_parser.py`/`engine_detector.py`, and the per-engine assertions in the contract are GNINA-tuned. Without parity, the multi-engine pipeline is a façade.

**Independent Test**: Run `python -m post_docking_analysis.simplified_cli --project-dir <vina_fixture>` and `<smina_fixture>` and `<autodock4_fixture>` and assert each produces the same `run_tracking/step_status.json` schema with identical step names and the same enforcement entries.

**Acceptance Scenarios**:
1. **Given** a Vina project with valid logs, **When** the simplified pipeline runs, **Then** affinity, RMSD, ProLIF and LigPlot all complete and the manifest matches the GNINA schema.
2. **Given** a Smina project missing logs, **When** the pipeline runs, **Then** it hard-fails with the same error class GNINA would raise.
3. **Given** an AutoDock4 project, **When** the user runs the HPC deploy/sync/submit subcommands, **Then** an engine-specific adapter exists analogous to `gnina_hpc_adapter.py`.

---

### User Story 3 — Interactive observability and ergonomics match simplified_cli guarantees (Priority: P2)

After any interactive run completes, the user immediately sees the path to `run_tracking/run_manifest.json`, a one-line summary of `stage_contract` (what was requested vs enforced), and which optional features ran or were skipped and why.

**Why this priority**: The hardening already writes these artifacts; they're just not surfaced. Low effort, high signal — fixes the "did anything actually run?" experience.

**Independent Test**: Complete one interactive run and capture stdout; assert it contains the absolute path to `run_manifest.json` and a contract diff line.

**Acceptance Scenarios**:
1. **Given** any completed analysis target, **When** the interactive workflow returns to the menu, **Then** the printed completion block names the manifest path and lists any `enforced != requested` entries.
2. **Given** the favorite-engine continuation flow, **When** the user opts in to interaction analysis, **Then** `analyze.interactions.clean` is selectable from the same prompt (not only from the per-stage submenu).
3. **Given** an optional feature failed (e.g. PoseView REST unreachable), **When** the run completes, **Then** the manifest distinguishes "skipped (disabled)" from "failed (error)" and the completion message mentions any failures.

---

### User Story 4 — Dead and duplicate code eliminated from post_docking_analysis (Priority: P2)

A maintainer running `grep` or `find` sees one canonical implementation per concern, no `_legacy_*` files masquerading as legacy, no stub analyzers, no plugin scaffolding, and no half-removed integrations.

**Why this priority**: The duplication and misleading naming caused the original gap (stubs running where real analyzers were expected). Cleaning it up prevents future regressions of the same shape and shortens onboarding.

**Independent Test**: `grep -rn "_legacy_" post_docking_analysis/` returns nothing; `grep -rn "rmsd_analyzer\|affinity_analyzer\|plugin_manager" .` only matches removed-file references in CHANGELOG/docs.

**Acceptance Scenarios**:
1. **Given** the repo state after this story, **When** a maintainer reads `post_docking_analysis/`, **Then** `_legacy_simplified_pipeline_impl.py` and `_legacy_multi_engine_pipeline_impl.py` are renamed (e.g. `simplified_pipeline_impl.py`, `multi_engine_pipeline_impl.py`) and the wrapper modules either become the canonical files or are deleted.
2. **Given** the cleanup, **When** the test suite runs, **Then** no import errors arise from removing `rmsd_analyzer.py`, `affinity_analyzer.py`, `plugin_manager.py`, `plugins/`.
3. **Given** a PandaMap decision, **When** the spec is closed, **Then** PandaMap is either fully removed (menu, imports, CLI flag, docs) or fully repaired (with a passing test) — not in the half-state it occupies today.

---

### User Story 5 — Configuration is unified across CLIs (Priority: P3)

A user supplies one config file (YAML, matching `post_docking_analysis/config/schema.yaml`) and it produces identical effective behavior across `simplified_cli`, `unified_pipeline`, the `analyze clean` shortcut, and the interactive workflow.

**Why this priority**: Today simplified_cli ignores config files and hardcodes defaults; cli.py requires them; interactive has its own prompt-driven state. Unification reduces "which knob actually applies?" confusion.

**Independent Test**: Apply the same YAML to each entrypoint against the same fixture project and diff the resulting `run_manifest.json` — should be byte-identical except for timestamps.

**Acceptance Scenarios**:
1. **Given** a config file, **When** any entrypoint is invoked, **Then** the same `config_manager` loader resolves it and the same precedence (CLI args > config file > defaults) applies everywhere.
2. **Given** an invalid config, **When** any entrypoint is invoked, **Then** the same validation error is raised with the same message.

---

### User Story 6 — Documentation reflects the current contract (Priority: P3)

A new user reads `README.md`, `post_docking_analysis/README.md`, `SIMPLIFIED_IMPLEMENTATION_STATUS.md`, `REDESIGN_PROPOSAL.md`, and `INTERACTIVE_CLI_DECISION_TREE.md` and finds a single coherent description that matches the live CLI flags and menu options.

**Why this priority**: Doc drift is annoying but not blocking. Doing it last (after code stabilizes) avoids re-writing the same files multiple times.

**Independent Test**: Manually walk each doc's example commands against the current CLIs; every command runs and produces the documented outputs.

**Acceptance Scenarios**:
1. **Given** the updated docs, **When** a reader follows the quickstart, **Then** every command exists and every flag is current.
2. **Given** the contradicting docs (REDESIGN_PROPOSAL marks PandaMap for removal, README still documents it), **When** the spec closes, **Then** the contradiction is resolved in code first, then in docs.

---

### User Story 7 — Regression test suite covers hard-fail branches and per-engine paths (Priority: P3)

Each hard-fail introduced by the recent hardening is backed by a test that fails if the hard-fail is weakened; each supported engine is backed by parser, contract, and interactive smoke tests.

**Why this priority**: Tests lock in the gains. Without them, the next refactor will silently re-introduce the original gaps.

**Independent Test**: `pytest test/` passes; intentionally setting `run_rmsd=False` in `_build_simplified_pipeline` or removing one of the hard-fail raises produces a failing test.

**Acceptance Scenarios**:
1. **Given** the new test suite, **When** RMSD's score-table-missing branch is bypassed, **Then** a dedicated test fails.
2. **Given** the per-engine fixtures, **When** the parser for any of {GNINA, Vina, Smina, AutoDock4} is broken, **Then** a dedicated test fails.
3. **Given** the interactive smoke harness, **When** the non-GNINA path silently sets `run_rmsd=False`, **Then** an interactive smoke test fails.

### Edge Cases

- **Pairlist auto-detection fails silently** beyond the documented 4-parent-level search: must raise or warn explicitly.
- **LigPlus root present but binary missing**: hard-fail should distinguish "tool not installed" from "tool produced zero diagrams".
- **ProLIF dependency present but RDKit conformer parse fails**: must be classified as a hard-fail, not silently skipped.
- **PoseView REST unreachable**: optional, must be recorded in manifest as `skipped (network)`, not `failed`.
- **Mixed-engine project directory** (GNINA + Vina results in the same project): engine_detector must report all engines and the pipeline must enforce the contract per engine, not just on the first detected.
- **User selects multiple per-stage targets in one interactive session**: contract enforcement must apply to the union, not per-target.
- **Legacy `PostDockingAnalysisPipeline` import attempts** from external scripts after removal: import error must include migration pointer to simplified_cli / analyze clean.

## Requirements *(mandatory)*

### Functional Requirements

**Contract enforcement (P1):**
- **FR-001**: System MUST enforce `run_rmsd=True` and `run_visualizations=True` in every code path that constructs `SimplifiedPostDockingPipeline`, including `_prepare_non_gnina_context` (`workflow/execution.py:2482`).
- **FR-002**: System MUST ensure ProLIF and LigPlot are jointly executed (or jointly refused) whenever either is requested via the interactive menu — no isolated `analyze.interactions.prolif`-only or `.ligplot`-only runs.
- **FR-003**: System MUST either migrate `analyze.stage.structure_quality` and `analyze.visuals.pymol` off the legacy `PostDockingAnalysisPipeline` class, or remove them from the interactive menu with a printed migration pointer.

**Engine parity (P1):**
- **FR-004**: System MUST provide log parsers and `engine_detector` coverage for Vina, Smina, and AutoDock4 equivalent to the GNINA implementation, with the same hard-fail surface.
- **FR-005**: System MUST provide HPC adapters for Vina, Smina, and AutoDock4 with the same interface as `gnina_hpc_adapter.py`.
- **FR-006**: System MUST emit the same `run_tracking/` schema regardless of engine.

**Observability (P2):**
- **FR-007**: System MUST print the absolute path of `run_tracking/run_manifest.json` and a one-line `stage_contract` diff at the end of every interactive analysis run.
- **FR-008**: System MUST distinguish in the manifest between `skipped (disabled)`, `skipped (missing dependency)`, and `failed (error: …)` for every optional feature.
- **FR-009**: System MUST surface `analyze.interactions.clean` as an option in the favorite-engine continuation flow, not only in the per-stage submenu.

**Cleanup (P2):**
- **FR-010**: System MUST rename `_legacy_simplified_pipeline_impl.py` and `_legacy_multi_engine_pipeline_impl.py` to non-misleading canonical names and remove or collapse the wrapper modules.
- **FR-011**: System MUST delete `post_docking_analysis/rmsd_analyzer.py`, `post_docking_analysis/affinity_analyzer.py`, `post_docking_analysis/plugin_manager.py`, and `post_docking_analysis/plugins/` (after migrating any still-used references).
- **FR-012**: System MUST make a binary decision on PandaMap (remove fully or repair fully) and reflect it in menu, CLI flags, imports, and docs.

**Configuration (P3):**
- **FR-013**: System MUST consume one canonical config loader (`config_manager` against `config/schema.yaml`) from every CLI entrypoint, including `simplified_cli`.
- **FR-014**: System MUST apply uniform precedence: CLI args > config file > defaults.

**Documentation (P3):**
- **FR-015**: System MUST update `README.md`, `post_docking_analysis/README.md`, `SIMPLIFIED_IMPLEMENTATION_STATUS.md`, `REDESIGN_PROPOSAL.md`, and `INTERACTIVE_CLI_DECISION_TREE.md` to match the post-remediation behavior.

**Tests (P3):**
- **FR-016**: System MUST include unit tests for every hard-fail branch added during the recent hardening.
- **FR-017**: System MUST include per-engine parser and contract tests for {GNINA, Vina, Smina, AutoDock4}.
- **FR-018**: System MUST include an interactive smoke test that drives `workflow/interactive.py` against a non-GNINA fixture and asserts the contract is enforced.

### Key Entities

- **Run Manifest** (`run_tracking/run_manifest.json`): canonical record of one pipeline invocation; fields include `stage_contract` (requested vs enforced settings), `step_status[]`, `outputs_index[]`, engine identity, timestamps, and per-optional-feature classification.
- **Stage Contract**: the set of mandatory and optional stages with their requested vs effective flags; produced by the simplified pipeline at `_legacy_simplified_pipeline_impl.py:156,471`.
- **Engine Adapter**: per-engine module (currently `gnina_hpc_adapter.py`) responsible for layout detection, log path conventions, and HPC submission helpers.
- **Analysis Target**: a string token like `analyze.interactions.clean` routed through `workflow/execution.py:run_analysis_target`; the menu in `workflow/interactive.py:64-80` lists all of them.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: 100% of interactive analysis runs across {GNINA, Vina, Smina, AutoDock4} produce a `run_manifest.json` whose `stage_contract` shows `enforced == requested` for RMSD, visualizations, ProLIF, and LigPlot (verified by smoke harness).
- **SC-002**: 0 imports of `rmsd_analyzer`, `affinity_analyzer`, `plugin_manager`, or `_legacy_*_impl` remain in the repo after the cleanup story closes (verified by `grep`).
- **SC-003**: Every supported engine has a passing parser test, a passing hard-fail contract test, and a passing HPC-adapter smoke test (verified by `pytest`).
- **SC-004**: Every interactive analysis completion message names the manifest path and contract diff (verified by capturing stdout in the interactive smoke test).
- **SC-005**: One config file produces identical `run_manifest.json` (modulo timestamps) when applied to `simplified_cli`, `unified_pipeline`, `analyze clean`, and the interactive workflow (verified by diff harness).
- **SC-006**: Quickstart commands in every updated doc execute successfully against a fresh checkout (verified by a documentation smoke test).
- **SC-007**: PandaMap reaches a terminal state (fully removed or fully working) with no mention of "broken", "doesn't work", or "marked for removal" remaining in the repo.

## Critical Files To Modify

**Contract enforcement & interactive alignment (P1):**
- `workflow/execution.py` (lines 2249-2250, 2354-2363, 2366-2376, 2482-2483)
- `workflow/interactive.py` (lines 64-80, 2507-2553)
- `post_docking_analysis/_legacy_simplified_pipeline_impl.py` (any non-GNINA branches)

**Engine parity (P1):**
- `post_docking_analysis/engine_detector.py`
- `post_docking_analysis/docking_parser.py`
- `post_docking_analysis/gnina_hpc_adapter.py` → split into per-engine adapters under `post_docking_analysis/hpc_adapters/`
- `post_docking_analysis/_legacy_multi_engine_pipeline_impl.py`

**Observability (P2):**
- `workflow/interactive.py:2507-2553` (completion message)
- `post_docking_analysis/_legacy_simplified_pipeline_impl.py:193,324,403` (manifest schema)

**Cleanup (P2):**
- Delete: `post_docking_analysis/rmsd_analyzer.py`, `post_docking_analysis/affinity_analyzer.py`, `post_docking_analysis/plugin_manager.py`, `post_docking_analysis/plugins/`
- Rename: `_legacy_simplified_pipeline_impl.py` → `simplified_pipeline_impl.py` (collapse wrapper)
- Rename: `_legacy_multi_engine_pipeline_impl.py` → `multi_engine_pipeline_impl.py` (collapse wrapper)
- `post_docking_analysis/pipeline.py:1366` (extend disable to class import if class targets are removed)
- PandaMap decision: `pandamap_integration.py`, `publication_pandamap.py`, `pandamap_runner.py` (remove) OR repair with a passing test

**Configuration (P3):**
- `post_docking_analysis/config_manager.py`, `post_docking_analysis/config/schema.yaml`
- `post_docking_analysis/simplified_cli.py` (consume the loader)
- `workflow/interactive.py` (consume the loader)

**Docs (P3):**
- `README.md`, `post_docking_analysis/README.md`, `post_docking_analysis/SIMPLIFIED_IMPLEMENTATION_STATUS.md`, `post_docking_analysis/REDESIGN_PROPOSAL.md`, `post_docking_analysis/INTERACTIVE_CLI_DECISION_TREE.md`

**Tests (P3):**
- `test/test_dockforge_smoke.py` (extend; new shortcut test already at line 1625 to mirror)
- New: `test/test_post_docking_contract.py` (per hard-fail branch)
- New: `test/test_engine_parity.py` (per-engine parsers + manifests)
- New: `test/test_interactive_smoke.py` (drive workflow/interactive.py with each engine)

## Existing Utilities To Reuse

- `_build_simplified_pipeline` (`workflow/execution.py:2243`) — the single chokepoint for `SimplifiedPostDockingPipeline` construction; FR-001 is implemented by making `_prepare_non_gnina_context` route through it instead of constructing directly.
- `_run_clean_interaction_pipeline_target` (`workflow/execution.py:1751`) — the clean pipeline entrypoint; FR-002 and FR-009 reuse this rather than introducing a new path.
- `run_tracking` writers (`_legacy_simplified_pipeline_impl.py:193,324,403`) — manifest/step_status/outputs_index already structured; FR-007/FR-008 only extend the printer in interactive.py and the classification taxonomy in the writers.
- `engine_detector.detect_engines` (`engine_detector.py`) — FR-004 extends rather than replaces.
- Smoke runner hook (`test/test_dockforge_smoke.py:4900`) — FR-016/FR-017/FR-018 register new tests against the same hook.
- Spec-kit templates (`.specify/templates/tasks-template.md`, `.specify/templates/plan-template.md`) — use to break this spec into task and plan files when implementation starts.

## Verification

End-to-end verification once all stories close:

1. **Per-engine fixture run**:
   - For each engine in {GNINA, Vina, Smina, AutoDock4}, run `python -m post_docking_analysis.simplified_cli --project-dir test/fixtures/<engine>_project --output /tmp/<engine>_out`.
   - Assert `/tmp/<engine>_out/run_tracking/run_manifest.json` exists; `jq '.stage_contract'` shows all mandatory stages `enforced == requested`.

2. **Interactive smoke**:
   - Run `pytest test/test_interactive_smoke.py -v`.
   - For each engine, pipe scripted answers to `python main.py`; assert exit code 0, manifest exists, stdout contains manifest path.

3. **Contract regression**:
   - Run `pytest test/test_post_docking_contract.py -v`.
   - Each hard-fail branch (no scores, RDKit unavailable, ProLIF zero maps, LigPlus root missing, LigPlot zero diagrams) has a dedicated test.

4. **Engine parity**:
   - Run `pytest test/test_engine_parity.py -v`.
   - Per-engine parser tests + manifest-schema-equality tests.

5. **Cleanup verification** (P2):
   - `grep -rn "_legacy_" post_docking_analysis/` → empty.
   - `grep -rn "from post_docking_analysis.rmsd_analyzer\|from post_docking_analysis.affinity_analyzer\|plugin_manager" .` → empty (or only in CHANGELOG/docs).
   - `grep -rn -i "pandamap" .` → either zero hits (removed) or only matches in working integration paths with a passing test (repaired).

6. **Config uniformity** (P3):
   - Run a diff harness: same YAML applied to all four entrypoints; manifests match modulo timestamps.

7. **Docs walkthrough** (P3):
   - Manually execute every code block in the updated docs against a fresh checkout; all commands succeed.

8. **`analyze clean` shortcut**:
   - `python main.py analyze clean --project-dir <fixture>` exits 0 and produces the documented command artifacts (`command_manifest.json/csv`, `command_log.txt`, `run_clean_interaction_pipeline.sh`).
   - Same scenario reached via interactive menu (`analyze.interactions.clean`) produces an equivalent run.
