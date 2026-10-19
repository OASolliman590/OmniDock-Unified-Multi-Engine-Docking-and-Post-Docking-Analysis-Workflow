# Spec 038 — Entry-Point Consolidation (three modules behind one CLI)

**Status:** directed by the Scientific Lead on 2026-10-19 ("the pipeline must be seamless and coherent"; consolidation after the merge).

**Goal:** one documented entry point, `main.py` / `pdb-prepare-wizard`, which drives three modules, with no parallel paths that can compute the same thing differently:
1. **Input / preparation:** `docking/preparation/`
2. **Docking:** `docking/`
3. **Post-docking:** `post_docking_analysis/`

`workflow/` (CLI, interactive shell, state, DAG) and `webui/` orchestrate those modules and compute nothing scientific themselves.

## Current redundancy (2026-10-19)
- Seven console scripts in `setup.py`, including the legacy interactive pipeline, `cli_pipeline`, `batch_pdb_preparation` and others. Top-level scripts include `cli_pipeline.py`, `core_pipeline.py` (1,269 lines), `batch_pdb_preparation.py`, `autodock_preparation.py`, `interactive_pipeline.py`, `test_pipeline.py` and `prep_autodock_enhanced.sh`.
- Two post-docking engines: `simplified_pipeline_impl.py` (5,183 lines) and `multi_engine_pipeline_impl.py` (6,538 lines), plus `pipeline.py`, `unified_pipeline.py`, `cli.py` and `simplified_cli.py`. The workflow already routes through a "simplified bridge" into the multi-engine implementation.
- Protonation overrides (Spec 037) are CLI-only.

## Requirements
- **R1 — One CLI.** `main.py` (`workflow`, `pdb`, `prep`, `dock`, `analyze`, `webui`) is the only supported entry point. Every other console script and top-level script becomes either:
  - (a) a thin deprecated shim that prints the equivalent `main.py` command and forwards to it, or
  - (b) removed when nothing uses it.
  An inventory decides which, per script, with callers found by grep and tests. Shims stay one release and are listed in `USAGE.md`.
- **R2 — Scientific code lives only in the three modules.** Any scientific computation still in top-level scripts moves into the matching module behind a function the CLI calls. Candidates: binding-site extraction and Excel export in `core_pipeline.py`, preparation configuration in `autodock_preparation.py`, and the shell receptor/ligand loop.
- **R3 — One post-docking path.** One analysis implementation is chosen and the other becomes a compatibility layer or is retired. The choice rests on which path carries the accepted Spec 031–037 methods, which is expected to be the multi-engine implementation used by the DAG. Its single-engine handling stays `single_engine_native_v1`. Duplicate functions (score parsing, best-pose writers, normalisers, summary writers) are removed in favour of the shared ones (`pose_selection.py`, `consensus.py`, `replicates.py`, `docking_parser.py`). Behaviour must not change: before/after golden outputs on the smoke fixtures must match, apart from documented intended differences.
- **R4 — Overrides everywhere.** Per-receptor and per-ligand protonation overrides (Spec 037 R1) can be set and reviewed in the interactive shell and the web UI, both calling the same resolver.
- **R5 — One documented path.** `USAGE.md` and `README.md` show the single sequence for one target and one ligand, and for many proteins and many ligands. The `--help` text matches.

## Non-goals
No change to any scientific method, threshold or default. This spec is structural. Any behaviour difference found during consolidation is reported to the Scientific Lead, not silently resolved.

## Acceptance
- An inventory table (script or module, callers, action) is committed.
- Every removed or shimmed path has a test showing the shim forwards, or that nothing imports the removed path.
- Golden-output comparison on the smoke fixtures before and after shows no unexplained difference.
- Full suite and both smoke tests pass, and the wheel check passes.
- Scientific Lead review.
