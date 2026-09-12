# Docking execution and engine contract audit

Repository commit: `cfe33056c77d817d19439a9e9affe225d00ca488`. Review date: 2026-09-12. Scope: active `workflow.execution -> docking.cli -> docking.runners` paths, local/Slurm/HTCondor deployment, pairlist/grid specification, score parsing. No pipeline implementation was modified. No AGENTS.md was found in the repository by the initial file search.

**Assessment:** the shared runner and pairlist architecture is a useful base, but its current run/result contract permits stale results, reports failed campaigns as successful, and changes scientific settings between local and deployed execution. Repair these before comparing engines, scaling campaigns, or interpreting improved ranks as improved science.

## Reproduced findings

The executable audit probe is [docking_probe.py](docking_probe.py); its captured results are [docking-probe-results.json](docking-probe-results.json). Run from the repository root with `python audit/docking_probe.py`. It creates only audit fixtures under `audit/probe-data/docking`; external engines are mocked for failure/output tests, and deployment scripts are generated without submitting jobs. Thus these are demonstrations of orchestration behavior, not docking accuracy benchmarks.

Severity: P1 = repair before scientific reuse/scaling; P2 = substantive correctness or robustness repair. Confidence is high unless qualified.

### D1 — P1: failed or changed runs can reuse and re-analyze old poses

- `docking/models.py:50-52` identifies a job solely as `receptor_site_ligand`. `docking/project_layout.py:199-241` does not incorporate round, seed, prepared-state hash, box, or engine settings into output paths.
- `docking/runners/base.py:130-142` skips any existing pose file, even an empty one. Planning the same pair with different parameters still identifies the same file as completed.
- `base.py:95-125` runs into the final output path, preserves an existing nonempty pair log, and collects scores after failures. `vina.py:63-69` and `autodock4.py:350-356` scan the pose directory without consulting successful jobs from the current invocation.
- **Probe:** a zero-byte pose was `skipped`. A mocked failed Vina rerun left an old `-9.0` result in `normalized_scores.csv`, with `OLD LOG` retained as its referenced pair log. The job manifest simultaneously said `failed`.
- **Impact:** an interrupted rerun or modified preparation can appear to have the previous favorable result; exhaustive and screening runs overwrite/reuse the same evidence. This can invalidate scientific attribution even when the numeric score is a real older score.
- **Fix:** immutable attempt directories keyed by experiment/job IDs, input hashes and resolved protocol; atomic publish after output validation; explicit links from score records to successful attempts. Resume only when that exact attempt has a validated completion record. Collect from the accepted-job manifest, not directory globs.
- **Qualification:** the comparative importer filters unknown pair tags (`post_docking_analysis/multi_engine_pipeline_impl.py:2439,2480,2520`). That protects against some unrelated files, but cannot distinguish old/new results for the same tag.

### D2 — P1: local execution can report a failed campaign as completed

- `docking/cli.py:582-618` prints failed counts but unconditionally returns `0`.
- `workflow/execution.py:1324-1333` maps this return code directly to workflow status `completed`.
- `docking/runners/base.py:110` defines success using only return code; it does not require a parseable, nonempty pose. `subprocess.run` startup exceptions also escape before the aggregate manifest is written.
- **Probe:** a Vina process returning `1` caused `dock_main` to return `0`. A GNINA process returning `0` with a numeric stdout table and no pose file was marked `completed` and produced normalized scores referencing absent SDF files.
- **Fix:** per-engine output validators and a typed result state (`succeeded`, `failed`, `invalid_output`, `interrupted`); checkpoint each attempt; catch startup errors; return nonzero for failed requested jobs, with a distinct partial-success report.
- **Regression criterion:** failed, zero-byte, truncated, missing-pose, and startup-error fixtures must never become completed scientific observations. Fewer poses than `num_modes` is not itself failure: that engine option is a maximum, also constrained by the energy window. [Official Vina FAQ](https://autodock-vina.readthedocs.io/en/latest/faq.html).

### D3 — P1: local and HPC paths silently use different scientific protocols

- Local execution resolves/applies the parameter schema and transforms box rows at `docking/cli.py:516-535`; it saves that schema at `:596`.
- Deployment loads the original pairlist at `:711-715` and builds fresh runtime defaults at `:759-762`. It does not load/apply the saved schema or transformed rows. Its CLI options at `:622-700` have no shared seed/box-scale/box-padding/schema controls.
- **Probe:** a local dry run with advanced parameters, `--seed 42 --box-scale 2`, generated a 40 Å box and `--seed 42`. Immediately deploying that project generated a 20 Å box with no seed flag; both entry points returned success.
- **Impact:** cluster results are not replications of locally reviewed jobs; comparison of search outcomes or engine performance is confounded by different boxes and stochastic protocols.
- **Fix:** resolve a single immutable `ExperimentSpec` before choosing an executor; local, Slurm and HTCondor consume exactly the same resolved jobs. Resource/scheduler settings must be separate from scientific settings. Save requested and effective box dimensions, explicit seeds, and all engine options.
- **Qualification:** HPC profile runtime dictionaries can supply a seed. The defect is failure to preserve/share the local protocol, not an assertion that no HPC seed can ever be configured. Exact reproducibility requires the same seed, inputs and parameters, as explained in the [Vina FAQ](https://autodock-vina.readthedocs.io/en/latest/faq.html).

### D4 — P1: portable AD4 deployment leaves absolute local paths inside generated GPF/DPF files

- The fallback path writes GPF/DPF during local planning (`autodock4.py:325-347`), embedding absolute receptor, ligand, map and parameter paths (`:149-165,196-209`).
- `docking/deployment.py:699-700,750-752` maps the command and result filenames, but does not rewrite the contents of those parameter files. HTCondor uses the same planning/mapping pattern at `:484-485,537-539`.
- **Probe:** an AD4 bundle with a different execution root mapped the command to that root, while its generated `grid.gpf` still contained the original local project path. This is portable-file generation failure independent of engine availability.
- **Impact:** after transfer to a host where the local root does not exist, AutoGrid/AutoDock cannot locate inputs/maps; if that local path does exist, the job may accidentally read an unintended copy.
- **Fix:** stage relative, whitespace-safe filenames inside a job working directory and execute there, or generate the files on the execution host from the resolved spec. Test all transitive file references after relocation.
- **Qualification:** the ADT-generator branch (`autodock4.py:256-324`) generates GPF/DPF on the execution host and has a different path lifecycle; this finding specifically targets the built-in fallback, used when ADT generators are not configured.

### D5 — P1: requested QC gates fail open on implementation/runtime errors

- `docking/preflight.py:123-183` wraps both ligand and receptor QC in `try/except` and downgrades audit exceptions to warnings. Return `ok` is based only on the errors list (`:219-228`).
- **Probe:** a requested receptor gate throwing `RuntimeError` returned `ok=True`, no errors, and a warning.
- **Impact:** broken dependencies, parsing failures, or internal QC defects can admit structures whose validity was never evaluated.
- **Fix:** represent `pass`, `fail`, and `not_evaluated` separately. A required gate must block on `not_evaluated`; an explicit recorded override can permit exploratory execution. Do not conflate optional drug-likeness alerts with structural validity.

### D6 — P2: GNINA CPU/GPU options do not implement the selected environment

- `docking/execution_environment.py:116-120` sets `use_gpu=False` for local CPU mode, but `docking/runners/gnina.py:71-78` never passes `--no_gpu`.
- The same block uses `if device ... elif cpu`, so choosing a GPU device discards the configured CPU count.
- **Probe:** `{use_gpu: False, cpu: 4}` omitted `--no_gpu`; `{use_gpu: True, device: 0, cpu: 4}` omitted `--cpu`.
- **Impact:** native GNINA can still use an available GPU during a requested CPU run; GPU jobs can exceed the intended CPU allocation. The metadata no longer describes the executed resources.
- **Fix:** emit `--no_gpu` when requested; pass `--cpu` independently of GPU device; derive its maximum from the scheduler allocation. Test each environment/flag combination against a pinned GNINA version. GNINA documents separate `--cpu`, `--device`, and `--no_gpu` controls. [Official GNINA usage](https://github.com/gnina/gnina#usage).

### D7 — P2: custom AutoDock parameters affect docking but not grid generation

- `autodock4.py:124-129,186-196` accepts a custom parameter file and writes it into DPF.
- `_write_gpf` at `:149-168` never writes a `parameter_file` entry; the ADT GPF command at `:266-282` also omits it, while the DPF command at `:293-294` includes it.
- **Probe:** `custom-atom-parameters.dat` appeared in DPF and not GPF.
- **Impact:** a custom atom type may fail grid generation; changed atom/force-field parameters can make intermolecular grids and docking calculations inconsistent. This does not assert that the default internal and bound parameter sets are always mismatched.
- **Fix:** use one parameter-set object/hash for both stages; check that every receptor/ligand atom type is supported; test standard and custom types. AutoGrid uses internal parameters unless its GPF supplies a parameter file. [AutoDock 4.2.6 guide, AutoGrid keywords](https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf).

### D8 — P2: AD4 silently shrinks the requested grid

- `autodock4.py:71-74` rounds `size / spacing` to an integer but does not require an even count. A 20 Å box at 0.375 Å spacing produces `npts 53 53 53`.
- AutoGrid truncates odd counts to the next lower even integer, so this becomes 52 intervals, or **19.5 Å**, not the requested 20 Å. This is not necessarily a fatal error. The behavior is explicit in [AutoGrid `check_size.cpp`, lines 72-75](https://github.com/ccsb-scripps/AutoGrid/blob/master/check_size.cpp#L72-L75).
- **Fix:** choose the next valid even grid interval count that covers the requested box, enforce version-specific limits, and record the effective dimensions. Test boundaries and anisotropic boxes. Some discretization difference is unavoidable; hidden contraction is avoidable.

### D9 — P2: GNINA's local normalized table includes duplicate runner-log pseudo-pairs

- `base.py:92,106-109` writes `tag.runner.log` with captured stdout in the same directory as `tag.log`.
- `post_docking_analysis/generate_scores_csv.py:179-181,260-270` discovers every `*.log`; unmatched pairlist names are warned about but their scores are still retained. `gnina.py:91-114` normalizes them without a known-pair filter.
- **Probe:** one GNINA job created both `R.pdbqt_site_1_L.pdbqt` and `R.pdbqt_site_1_L.pdbqt.runner` score rows. The latter had blank receptor/ligand metadata and referenced a nonexistent `.runner.sdf`.
- **Fix:** parse only manifest-designated engine outputs; parse SDF property records for pose-associated scores, and retain logs as diagnostic/provenance data. Require a resolvable pose ID for every accepted score.
- **Qualification:** the comparative raw importer filters these unknown tags when a pair index is available. The malformed local `normalized_scores.csv` contract is still real; do not overstate this as inevitable double-counting in every final report.

### D10 — P2: geometry and missing-metric validation are insufficient

- `docking/preflight.py` validates asset presence and optionally QC, but never checks finite box coordinates, positive finite dimensions, pocket coverage, or ligand fit. `parameter_schema.py:102-105` uses comparisons that admit NaN, and no `isfinite` check.
- **Probe:** a pair with NaN center, negative x size, infinite y size, and zero z size passed preflight. A schema with NaN scale and infinite padding had no validation errors.
- `post_docking_analysis/docking_parser.py:84-96` changes an absent Smina `minimizedRMSD` into two measured-looking zero values. **Probe:** an affinity-only Smina pose parsed to `rmsd_lb=0.0,rmsd_ub=0.0`.
- **Fix:** finite/unit-aware geometry validation plus ligand/reference-box containment checks; preserve missing metrics as null. Separate RMSD-to-input, RMSD-to-top-pose bounds, symmetry-corrected reference RMSD, and aligned shape RMSD. These describe different quantities and should not share a generic field simply because they have Å units.
- **Scientific reason:** a valid pocket-centered box must allow the relevant ligand coordinates; too-large boxes require greater search effort. [Official Vina search-space guidance](https://autodock-vina.readthedocs.io/en/latest/faq.html). For Smina local minimization, the engine author's explanation confirms that minimizedRMSD measures movement caused by refinement, not universally the Vina lower/upper-bound metric. [David Koes, Smina support discussion](https://sourceforge.net/p/smina/discussion/help/thread/158a12d859/).

## Additional technical observations

- **Windows portability:** initial fixture creation actually failed with `WinError 1314` at `docking/project_layout.py:299`, which unconditionally creates a directory symlink. The probe precreated a legacy compatibility directory to continue testing without modifying production code. README line 144 includes Windows installation guidance. Graceful compatibility handling and an explicit supported-platform matrix are needed.
- **Remote path semantics:** `deployment.py:663-666` uses host `Path` for a remote Linux root; under Windows, `/remote/project` is not a fully absolute Windows path and is joined to the local project. `_map_project_path` also preserves backslash tails. Use `PurePosixPath` for Linux execution targets and host paths only for local storage.
- **Pair identity:** `pairlist_builder.py:78-124` resolves normalized aliases with `setdefault`, silently choosing the first match if multiple files collapse onto the same normalized name/PDB ID. Exact filenames work; ambiguous aliases should error with candidates. Generated site metadata is effectively one site per PDB ID, with fixed cubic default size (`:198-208,253-260`); it cannot represent multiple site instances or ligand-dependent box extent without custom pairlists.
- **Parameter scope:** the shared presets change only Vina-family `exhaustiveness/num_modes`; AutoDock4 reads separate GA controls. Do not interpret a common "exhaustive" label as equal search effort across different algorithms. Derive engine-specific stopping budgets from validation curves.
- **Failure-state coverage:** `test/test_dockforge_smoke.py` has many workflow and synthetic-output checks, but the reproduced failures above show that contract coverage is incomplete. Large smoke files alone do not establish scientific correctness.

## Strengths worth preserving

1. A shared explicit pairlist and small runner interface make cross-executor parity achievable without rewriting every analysis module.
2. Slurm/HTCondor generation already reuses `plan_jobs`, and scripts track shell failures. Preserve this common planning seam while strengthening the experiment/result types.
3. AD4 has engine-specific prepared-asset resolution, torsional-degree extraction, explicit GA settings, and both built-in and ADT-assisted parameter generation.
4. Native Vina/Smina command construction uses argument lists; Vina output parsing preserves its separate lower/upper RMSD values. Smina supports custom scoring selection.
5. Runtime manifests, pair-intent metadata, explicit cocrystal-benchmark flags, and optional seeds are useful foundations for provenance.
6. Ensemble receptor support is honestly documented as a scaffold (`models.py:135-145`, `base.py:40-58`); it is not falsely implemented as a working multiconformer engine.

## Scientific redesign priorities

**First establish one reproducible experiment contract.** Include prepared receptor state, ligand protomer/tautomer/stereoisomer/conformer, coordinate-frame ID, site/reference identity, requested/effective grid, engine build/container digest, scoring function/CNN model identity, explicit seed, search budget, and output atom mapping. Executors change resource allocation and transport only. Scores carry a direction, unit, semantic definition, model identity and exact pose/attempt link.

**Separate pose generation, pose selection, and affinity ranking.** GNINA, Smina and Vina share ancestry; agreement between them is not three independent experimental confirmations. GNINA's default pose ordering uses CNNscore, while CNNaffinity and empirical energy answer different questions. Preserve all scores and validate a selection policy per target rather than flattening them into one generic binding-energy field. [GNINA repository and usage](https://github.com/gnina/gnina).

**Calibrate each target before screening.** Use held-out reference complexes, symmetry-aware heavy-atom redocking RMSD in the receptor frame, cross-docking when receptor-state transfer matters, multiple explicit seeds, and actives/decoys or measured ligands for enrichment/ranking evaluation. Report uncertainty and failure rates, not just the most favorable seed. Choose exhaustiveness/grid padding from convergence and pose-recovery data. The distinction between sampling failure and scoring failure is explained in the [Vina FAQ](https://autodock-vina.readthedocs.io/en/latest/faq.html).

**Make grids auditable scientific assets.** Store a typed box with source reference/pocket, coordinate frame and visualization; verify reference/ligand coverage, finite values and engine-effective dimensions. Cache AD4 maps by receptor-state hash, site, spacing, atom-type union and parameter hash instead of recomputing identical receptor/site maps for each ligand.

**Use independent replicates and receptor states deliberately.** Add explicit replicate and receptor-state axes to the task matrix. Do not bolt ensemble execution onto the current filename tag or interpret arbitrary docking energies as calibrated Boltzmann weights; aggregation needs validated score meaning and an explicit scientific policy.

## Acceptance tests for the next implementation

1. The same resolved job produces identical scientific arguments/config content across local, Slurm and HTCondor, including after relocation.
2. Changing receptor/ligand bytes, seed, box, scoring model, or engine version changes run identity and cannot hit an old completion cache.
3. Failed/interrupted processes and missing, truncated or invalid pose files never enter accepted result tables; completed rows always link to a present validated pose.
4. Required QC exceptions block launch; optional filters retain explicit outcomes and overrides.
5. AD4 default/custom parameter files are consistent in GPF/DPF, grid intervals are valid, and emitted dimensions cover the requested box.
6. Parsers run against real pinned-engine output fixtures, including no-CNN GNINA, Smina missing/-1 RMSD, AD4 partial DLG, and multi-model poses. Missing values remain missing and one accepted pose has one score record per metric.
7. A small curated end-to-end benchmark suite measures redocking success, score-direction correctness, output-coordinate preservation, replicate stability and enrichment. No real-engine or receptor/ligand chemistry benchmark was run in this audit.
