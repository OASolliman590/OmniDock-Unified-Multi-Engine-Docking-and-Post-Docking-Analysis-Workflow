# Spec 038 T002: entry-point inventory

Date: 2026-10-19. Base: HEAD 8070791 (branch chore/install-agentic-bioinformatics-stack).
Callers were found by grep over code, tests, docs and setup/MANIFEST/check_wheel. `docs/` release checklists and `CHANGELOG.md` are historical and were not edited.

Action vocabulary: **keep** (supported), **shim** (deprecated, prints notice and forwards to `main.py` or refuses), **move-code-then-shim** (scientific helpers moved to `docking/preparation/`, script entry point shimmed), **deferred** (not converted; needs a decision, see STOP list).

## 1. Console scripts (`setup.py`)

| Console script | Target | Callers | Scientific code? | Action (this change) |
|---|---|---|---|---|
| `pdb-prepare-wizard` | `workflow.cli:main` | setup.py; test_spec031_operational (asserts this line) | no (CLI) | keep |
| `pdb-wizard-workflow` | `workflow.cli:main` | setup.py | no | keep |
| `pdb-wizard-interactive` | `workflow.interactive:run_interactive_workflow` | setup.py; workflow shell | no (shell) | keep (shell owned by Agent C) |
| `pdb-wizard-legacy-interactive` | was `interactive_pipeline:run_interactive_pipeline` | setup.py; `workflow/execution.py` imports `run_interactive_pipeline` (library) | yes (legacy interactive PDB flow) | **shim**: now `interactive_pipeline:main` -> `pdb legacy-interactive` |
| `pdb-wizard-cli` | `cli_pipeline:main` | setup.py; INSTALLATION_GUIDE, README (old) | yes (legacy batch-of-PDB run) | **shim**: `cli_pipeline:main` -> `pdb run`; `--create-config` refused |
| `pdb-wizard-batch` | `batch_pdb_preparation:main` | setup.py; INSTALLATION_GUIDE, README (old) | yes (batch preparation) | **shim**: `batch_pdb_preparation:main` -> `pdb batch` |
| `pdb-wizard-prepare-docking` | was `docking.cli:prepare_docking_main` | setup.py; `workflow/execution.py` calls `prepare_docking_main` for `prep project` | no (docking project build; library stays in `docking/cli.py`) | **shim**: `workflow.legacy_entrypoints:prepare_docking` -> `prep project` (flags verified identical) |
| `pdb-wizard-dock` | was `docking.cli:dock_main` | setup.py; `workflow/execution.py` calls `dock_main` for `dock run`; test_docking_correctness | no (library stays in `docking/cli.py`) | **shim**: `workflow.legacy_entrypoints:dock` -> `dock run`; 14 legacy-only flags refused |
| `post-docking-analysis` | `post_docking_analysis.cli:main` | post_docking_analysis/USAGE.md; README "Post-Docking Analysis Usage" (`python -m post_docking_analysis`); skills/post-docking-contract-review; test_dockforge_smoke imports `_prompt_cli_engine_fallback` | **yes** (second post-docking path with its own splitting, apo, visualization and analysis steps) | **deferred** (STOP). Its behaviour is not covered by `analyze comparative`; routing is decided by Spec 038 R3 / T004 (Agent B). Entry point left unchanged. |

## 2. Top-level scripts and modules

| File | Holds | Callers (code, tests) | Scientific? | Imported by workflow (names) | Action |
|---|---|---|---|---|---|
| `main.py` | argv normaliser for legacy aliases (`interactive`, `cli`, `batch`, `prepare-docking`, `dock`) | every entry point, tests | no | — | keep (only supported entry) |
| `cli_pipeline.py` | `run_single_pdb_cli`, `run_batch_analysis`, `parse_pdb_list`, `load_config`, `_select_ligand_instance` (ligand selection), `create_sample_config`; `main` (argparse) | `workflow/execution.py` (`run_pdb_cli`, `run_pdb_collect`); test_preparation_correctness (`run_single_pdb_cli`) | **yes** (ligand selection, cleaning and analysis orchestration) | `parse_pdb_list`, `load_config`, `run_batch_analysis`, `run_single_pdb_cli`, `EXCEL_AVAILABLE` | **move-code-then-shim**: Excel writer -> `docking/preparation/excel_export.py`; `main` -> shim. Orchestration and ligand selection **not moved** (STOP, see section 4). |
| `core_pipeline.py` | `MolecularDockingPipeline` (fetch, clean, hetatm, pocket analysis, PLIP text), `extract_residue_level_coordinates`, `parse_plip_text_report`, `centre_method_from_provenance`, `excel_cell_value`, `EXCEL_AVAILABLE` | workflow/execution (`MolecularDockingPipeline`); cli/interactive/batch; tests (test_preparation_correctness, test_spec031_remediation, test_spec036_docking_params) | **yes** | `MolecularDockingPipeline` | **move-code** (R2): binding-site extraction -> `docking/preparation/binding_site_extraction.py`; Excel helpers -> `docking/preparation/excel_export.py`. `MolecularDockingPipeline` class stays (STOP). Old names re-exported. |
| `interactive_pipeline.py` | `run_interactive_pipeline` (legacy interactive flow), `run_single_pdb_analysis`, prompts; `add_to_excel_workbook` (duplicate); `main` (added) | workflow/execution (`run_legacy_pdb_interactive`, `run_pdb_collect`) | **yes** (interactive selection flow) | `run_interactive_pipeline`, `run_single_pdb_analysis` | **move-code-then-shim**: Excel writer (was duplicated) -> `excel_export`; `main` -> shim. Flow stays (library). |
| `batch_pdb_preparation.py` | `BatchPDBPreparationPipeline`; `main` (argparse) | workflow/execution (`run_pdb_batch`); INSTALLATION_GUIDE, README | **yes** (batch preparation) | `BatchPDBPreparationPipeline` | **shim** (`main` -> `pdb batch`). Class stays (STOP). |
| `autodock_preparation.py` | `PreparationConfig`, `profile_requires_autodocktools` (-> moved); `AutoDockPreparationPipeline` (class, stays); `main` (argparse) | workflow/execution (`run_autodock_prepare`); batch_pdb_preparation; test_spec033_receptor_prep; test_dockforge_smoke (`import autodock_preparation`) | **yes** (ligand/receptor preparation orchestration) | `AutoDockPreparationPipeline`, `PreparationConfig` | **move-code-then-shim**: config -> `docking/preparation/preparation_config.py` (re-exported); `main` -> shim. Pipeline class stays (STOP). |
| `test_pipeline.py` | `run_smoke_test`, unified all-engines smoke runner (runs `main.py` subprocesses) | test/test_dockforge_smoke.py (`from test_pipeline import run_smoke_test`); docs checklists | no (test harness) | — | **move** to `test/test_pipeline.py` (repo root path fixed; smoke import still resolves via `test/` on sys.path) |
| `prep_autodock_enhanced.sh` | shell backend: receptor/ligand loop, PDBQT/PDB2PQR calls | `autodock_preparation.py` (`enhanced_script`); test_spec033_receptor_prep; test_spec037_r4a_prepare_once; setup.py build step; MANIFEST.in; check_wheel | **yes** (receptor/ligand loop) | — | **keep as internal backend**; header now says it is not a user entry point. Moving the loop to Python is **not done** (STOP, see section 4). |
| `webui_cli.py` | `main` -> `webui.app` launcher | `workflow/cli.py` (`main.py webui` imports it); test/test_webui.py; docs/WEB_UI_GUIDE.md; start_webui.bat (indirect) | no | — | **not in this assignment** (webui is Agent C area). Recommend: `main.py webui` is already the supported launcher; `webui_cli.py` becomes a shim once `workflow/cli.py` stops importing it. |
| `install_vina_dock_dependencies.sh` | one-off conda installer for a named `vina_dock` env | none in code; README/INSTALLATION not linking it | no | — | **not in this assignment**; environment helper, not an entry point. Recommend the lead decide keep/remove. |
| `scripts/check_wheel.py` | wheel acceptance check | Spec 038 gate; test_spec038 wheel test | no | — | keep |

## 3. What workflow already imports from top-level modules (must keep working)

From `workflow/execution.py` (all import paths verified by the test suite):
- `cli_pipeline`: `parse_pdb_list`, `load_config`, `run_batch_analysis`, `run_single_pdb_cli`, `EXCEL_AVAILABLE`
- `core_pipeline`: `MolecularDockingPipeline`
- `interactive_pipeline`: `run_interactive_pipeline`, `run_single_pdb_analysis`
- `batch_pdb_preparation`: `BatchPDBPreparationPipeline`
- `autodock_preparation`: `AutoDockPreparationPipeline`, `PreparationConfig`
- `docking.cli`: `prepare_docking_main`, `dock_main`

Also imported by tests or other modules: `core_pipeline.excel_cell_value`, `core_pipeline.extract_residue_level_coordinates`, `core_pipeline.MolecularDockingPipeline`, `cli_pipeline.run_single_pdb_cli`, `autodock_preparation.PreparationConfig`, `test_pipeline.run_smoke_test`.

All of these still resolve. `test/test_spec038_entry_points.py` asserts identity for the moved names.

## 4. Not moved (STOP: needs Scientific Lead or architecture decision)

These contain scientific computation or orchestration that R2 names only as candidates, or that is too large to move without changing structure. They are left in place with no change to behaviour:
1. `MolecularDockingPipeline` class (`core_pipeline.py`): PDB fetch, cleaning, HETATM handling, pocket analysis. The binding-site and Excel parts moved; the rest is a large class used by `workflow/execution.py`.
2. `run_single_pdb_cli`, `run_batch_analysis`, `_select_ligand_instance` and the non-anchor constants (`cli_pipeline.py`). Ligand selection is a scientific choice.
3. `BatchPDBPreparationPipeline` (`batch_pdb_preparation.py`) and `AutoDockPreparationPipeline` (`autodock_preparation.py`).
4. The receptor/ligand loop in `prep_autodock_enhanced.sh` (R2 candidate). Moving it changes execution order or tool calls, so it is not done here.
5. `post-docking-analysis` console script and the second post-docking path (R3 / T004, Agent B).

## 5. Summary

- Console scripts: 9. Keep: 3 (`pdb-prepare-wizard`, `pdb-wizard-workflow`, `pdb-wizard-interactive`). Shim: 5 (`pdb-wizard-legacy-interactive`, `pdb-wizard-cli`, `pdb-wizard-batch`, `pdb-wizard-prepare-docking`, `pdb-wizard-dock`). Deferred: 1 (`post-docking-analysis`).
- Top-level scripts in scope: 7. Shim `main` in 4 modules. Move `test_pipeline.py` to `test/`. Keep `prep_autodock_enhanced.sh` as internal backend. Scientific helpers moved: binding-site extraction, Excel writer (deduplicated), preparation config.
- Known gap found while mapping flags: `pdb-wizard-dock` (`dock_main`) accepts 14 flags that `main.py dock run` does not (`--ligand-qc-*`, `--no-ligand-admet-filters`, `--autodocktools-*`). The shim refuses them. Users who relied on them must move to a `main.py` path that accepts them, which does not exist yet. Decision needed.
