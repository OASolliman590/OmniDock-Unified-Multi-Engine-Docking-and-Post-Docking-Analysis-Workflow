# OmniDock technical and developmental audit

Audited on 2026-09-12 at commit `cfe33056c77d817d19439a9e9affe225d00ca488`. This is a whole-repository assessment, not a diff review. No tracked implementation files were changed. Findings marked **reproduced** use synthetic fixtures, without executing molecular docking engines. P1 means a high-priority defect that can corrupt results, misreport execution, or prevent a supported workflow; P2 means a significant reliability or distribution problem.

## Assessment

The repository is an ambitious research prototype with valuable infrastructure. It has four engine adapters, explicit pairing and preparation contracts, manifests, preflight checks, a scope-driven artifact graph, checkpoint metadata, and a SQLite export option. These are useful foundations. The current implementation does not yet justify an unattended, scientifically validated production claim.

The maintainability problem is the spread of scientific decisions across orchestration, selection, reporting, and compatibility paths. For example, score ordering is independently implemented in the consensus module, the top-pose atlas, and the multi-engine implementation. An interface that declares a metric's direction once would prevent several scientific defects found in this audit.

The repository contains 105 Python files and 60,677 Python lines, including a 6,446-line smoke script. The multi-engine and simplified implementation files contain 5,913 and 5,196 lines respectively. Size alone is not a defect; here those files combine parsing, scientific selection, caching, file export, visualization, metadata and execution, which makes a scientific correction difficult to apply consistently. No tracked `.github` CI workflow was found at the audited commit.

## T1 — P1: normal reruns can use stale scientific data (reproduced through actual pipeline)

`raw_scores` declares only pairlist and engine-scope metadata as inputs in [multi_engine_pipeline_impl.py:1401](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/post_docking_analysis/multi_engine_pipeline_impl.py#L1401). Actual engine output files are not dependencies. The artifact executor hashes path, modification time and size for declared inputs, rather than content, in [artifact_graph.py:171](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/post_docking_analysis/artifact_graph.py#L171).

The integrated probe uses the repository's actual Vina fixture, parser and graph registration. First parse: −7.2. Change the PDBQT score to −11.2, leaving pairlist/scope unchanged. Normal graph request: `cache_hit`, still −7.2. `force=True`: −11.2. This is a correctness failure in the normal incremental path, independent of the engine runner's separate stale-file problem.

Other dependency omissions have the same shape: validation consumes reference and pose files but declares only the normalized CSV; biology file content is not directly declared; several nodes accept directories whose own stat does not summarize nested file content. A separate manifest writer already computes SHA-256 checksums, but those checksums are not connected to this cache contract.

**Remedy:** inventory content-addressed inputs; include raw poses/logs, chemistry/reference/receptor assets, effective parameters, executable/model identity, code/schema version and all downstream dependencies in a run fingerprint. Validate outputs before marking a node reusable. Mutation tests must change each scientific input independently, including reference coordinates and biology table contents. Cache reuse must produce the same result as recomputation.

## T2 — P1: graph completion does not guarantee valid dependency/output contracts (reproduced)

In [artifact_graph.py:294](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/post_docking_analysis/artifact_graph.py#L294), failures propagate only from nodes with status `failed`, not from `blocked_by_failure`. In an A → B → C chain, the probe gives A=`failed`, B=`blocked_by_failure`, C=`completed`. Descendants can therefore execute with absent or stale prerequisites. The unified runner does inspect the original failed node and returns failure, so this counterexample does **not** establish that the whole run is reported successful; it establishes incorrect descendant execution and artifacts.

Separately, a callback returning `False` without creating its declared output is marked `completed` at [artifact_graph.py:269](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/post_docking_analysis/artifact_graph.py#L269). The probe's requested output does not exist. `None`/`False` return values, output existence and output validity have no enforced completion contract. Optional exceptions are represented as `skipped_optional`; missing output and legitimate absence need distinct states.

**Remedy:** return a typed stage result; propagate blocked states transitively; check required external inputs before execution and validate declared output schemas before commit. Optional nodes must expose explicit unavailable results, and consumers must declare whether they tolerate those results. Preserve overall execution failure separately from scientific validation failure.

## T3 — P1/P2: concurrent workflow updates lose data (reproduced)

[workflow/state.py:232](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/workflow/state.py#L232) and [284](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/workflow/state.py#L284) read state, update a private snapshot and save it. An `RLock` protects the read and write individually, not the transaction. A deterministic two-thread probe synchronizes both reads, then writes two independent context keys. Only one survives.

This can affect the actual background task manager, which permits concurrent tasks. An atomic file replacement prevents partial JSON; it cannot prevent lost updates. An in-process lock also cannot coordinate separate CLI processes.

**Remedy:** lock the complete read–modify–write operation immediately; use a process-safe transaction/lock for multiple CLI processes. For the durable design, use a transactional run/task store with immutable event records and derive the current state from it. Keep context convenience separate from experimental provenance.

## T4 — P2: installed console commands omit their Python modules (static, high confidence)

[setup.py:23](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/setup.py#L23) explicitly sets `packages=find_packages()` but supplies no `py_modules`. Several console entries point to root-level modules, including `main`, `cli_pipeline`, `interactive_pipeline` and `batch_pdb_preparation`; root preparation modules are also used by packaged workflow code. Those files are not packages and are not explicitly included as Python modules.

Setuptools disables automatic module discovery when `packages` is explicitly supplied. Consequently, this configuration does not include those standalone modules in a normal built distribution; editable/source-tree execution can conceal the problem. This conclusion is from configuration and [official package-discovery behavior](https://setuptools.pypa.io/en/stable/userguide/package_discovery.html), not an installed-wheel experiment in this environment.

**Remedy:** move runnable modules into one installable `omnidock` package, use a `pyproject.toml`, and test a built wheel from a clean directory. Declare preparation, analysis and visualization extras. Resolve the current contradiction between mandatory `plip` in `requirements.txt` and optional PLIP claims in documentation. Pin and test scientific tool versions/container digests; broad minimum Python dependency versions are not a reproducible engine environment.

## T5 — P2: the test entrypoints do not match the documented pytest expectation (reproduced)

The smoke module contains 100 `_smoke_*` functions and a manual `main()`. There are no `test_*` functions or pytest configuration enabling `_smoke_*` discovery. `python -m pytest --collect-only -q`, with unrelated globally installed pytest plugins disabled, reported **no tests collected**. This matches [pytest's documented discovery rules](https://docs.pytest.org/en/stable/explanation/goodpractices.html). It does not mean no tests exist: the repository explicitly supports executing its custom smoke script.

The documented custom command `python test/test_dockforge_smoke.py --skip-prep-matrix --skip-all-engines` passed the first three groups, then stopped at `checkpoint listing should be newest-first`. The likely explanation is timestamp granularity: `_now_iso()` removes microseconds, while checkpoints created within the same second are sorted only by `created_at`, preserving old-first insertion order. This is a lower-priority issue, but it prevents the remainder of the monolithic run from executing.

To inspect important areas despite that early exit, eight existing smoke functions were invoked independently: seven passed; the graph-registration check failed because Windows symbolic-link privileges were unavailable. The passing set includes consensus, geometric alignment and multi-pose parsing contracts. Some passing tests explicitly encode incorrect scientific assumptions, such as accepting truncated atom mappings and ascending weighted-hybrid selection. Test success therefore cannot currently establish scientific validity.

**Remedy:** convert the harness into independently discoverable tests and add CI for a clean Linux scientific environment plus supported Windows behavior. Maintain a tiny suite of scientific invariants, parser golden files from real versioned tools, a release-level multi-engine integration suite and an external scientific benchmark. Assert outcomes from an independently specified oracle, not the implementation's current choices.

## T6 — P2: fresh project initialization assumes symbolic-link privileges (reproduced)

[docking/project_layout.py:282](https://github.com/OASolliman590/OmniDock-Unified-Multi-Engine-Docking-and-Post-Docking-Analysis-Workflow/blob/cfe33056c77d817d19439a9e9affe225d00ca488/docking/project_layout.py#L282) unconditionally creates a compatibility symlink when the legacy directory is absent. A fresh project failed with Windows `WinError 1314` even outside the filesystem sandbox. The integrated cache probe therefore used the already-supported existing-directory branch by precreating the legacy directory; no implementation was monkeypatched for that result.

**Remedy:** make the compatibility alias optional, with a clear redirect file or documented platform fallback. A convenience alias should not prevent creating a canonical project. Do not claim OS-independent operation solely from Python portability.

## Design direction

Preserve the real engine-adapter seam and reduce the number of places that make scientific decisions. Introduce small public interfaces for:

1. `prepare_structure` — chemical identity, receptor/ligand state, exclusions, atom mapping, coordinates and a preparation audit trail.
2. `compile_run` / `execute_run` — one immutable resolved run specification shared by local and HPC adapters, plus a manifest of successfully validated outputs.
3. `read_poses` — structured coordinates plus preserved molecular graph, source pose identity and score records.
4. `validate_poses` — explicit chemical, physical and reference-pose tests, each with passed/failed/not-evaluable and a reason.
5. `rank_candidates` — typed metric direction, units and scope; coverage-aware aggregation; output labels backed by validation evidence.
6. `publish_report` — renders validated records and provenance, without deciding scientific rankings again.

These are modules with deep implementations behind narrow interfaces. They need not become separate services. A single installable Python package plus adapters is enough initially. Avoid a wholesale rewrite before the counterexamples are locked into regression tests. Each scientific correction should update one implementation used by every command surface.

## Verification record and limits

- `technical_probes.py`: four behavioral counterexamples plus repository inventory; all reproduced. Temporary fixtures only.
- `integrated_technical_probes.py`: actual pipeline cache defect; eight selected existing smoke checks, 7 pass / 1 Windows privilege failure. See `integrated-technical-results.json`.
- Full existing smoke command: stopped at checkpoint ordering; not a successful complete suite. See `existing-smoke-run.log`.
- Default pytest discovery: no collected tests after disabling unrelated globally installed plugins; the first attempt was interrupted by a host plugin/DLL policy issue.
- Chemistry and engine-specific execution checks are documented in the companion scientific audits. No heavy dependency installation, real docking benchmark, HPC submission, wheel installation, timing benchmark or experimental affinity validation was performed. The host lacks RDKit, Open Babel, Meeko/PLIP and Biopython in the inspected Python environment; their absence here is an audit limitation, not by itself a repository defect.

Companion files: `preparation-audit.md`, `docking-audit.md`, `postdocking-audit.md`, and `REVOLUTION_PLAN.md`.
