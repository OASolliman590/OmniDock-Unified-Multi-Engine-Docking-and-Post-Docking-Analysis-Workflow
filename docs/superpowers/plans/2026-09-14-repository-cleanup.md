# OmniDock Repository Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Clean OmniDock's repository, simplify non-scientific orchestration modules, and replace contradictory documentation with one verified active documentation set.

**Architecture:** Preserve public CLI and deployment interfaces while putting repository checks, parser construction, and common deployment planning behind small module interfaces. Historical evidence remains accessible under explicit archive/history labels; active guides route from one documentation index.

**Tech Stack:** Python 3.10+, argparse, setuptools, pytest, GitHub Actions, Markdown, PowerShell/Git Bash for local verification.

**Spec:** `docs/superpowers/specs/2026-09-14-repository-cleanup-design.md`

## Global Constraints

- Fixed review base: `8d4434d3a33eb83b1e12cad82944b02c83270e47`.
- Grok implements one bounded batch at a time without staging or committing; Sol reviews, verifies, and commits.
- Preserve every scientific invariant listed in the design specification.
- Do not edit scientific preparation/post-docking execution files reserved by the parallel scientific-integrity task.
- Preserve public CLI flags, defaults, console scripts, and deployment function signatures unless the specification names an additive compatibility-safe change.
- Do not delete research/manuscript data, benchmark outputs, curated audits, structured probe results, credentials, untracked work, or files outside this isolated checkout.
- Never claim historical validation as current; never invent benchmark or scientific guarantees.

---

### Task 1: Repository contract checker

**Files:**
- Create: `scripts/check_repository.py`
- Create: `test/test_repository_contracts.py`

**Interfaces:**
- Consumes: a repository root and Git's tracked-file list.
- Produces: `check_repository(root: Path) -> list[Finding]` and CLI exit code `0` only when no finding exists.

- [ ] **Step 1: Write the failing checker-interface test**

```python
def test_repository_checker_reports_broken_links_generated_logs_and_dangling_links(tmp_path):
    from scripts.check_repository import check_repository
    # Create a tiny Git fixture with README.md -> missing.md, audit/run.log,
    # and a dangling symlink when the platform supports symlinks.
    findings = check_repository(tmp_path)
    assert {finding.code for finding in findings} >= {"broken-markdown-link", "tracked-generated-artifact"}
```

- [ ] **Step 2: Run the focused test and verify RED**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test/test_repository_contracts.py -q`
Expected: collection fails because `scripts.check_repository` does not exist.

- [ ] **Step 3: Implement the minimal checker**

```python
@dataclass(frozen=True)
class Finding:
    code: str
    path: str
    detail: str

def check_repository(root: Path) -> list[Finding]:
    tracked = tracked_paths(root)
    return [*broken_markdown_links(root, tracked), *tracked_generated_artifacts(tracked), *dangling_symlinks(root, tracked)]
```

The Markdown resolver must ignore external schemes and in-page anchors, decode `%20`, strip fragments, and resolve paths relative to the source document. Generated-artifact detection is limited to the named audit verification outputs in the specification; it must not flag structured probe JSON or user data.

- [ ] **Step 4: Run focused tests and verify GREEN**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test/test_repository_contracts.py -q`

- [ ] **Step 5: Verify its current findings**

Run it locally and record the expected base findings before Task 2; do not weaken the checker to make them disappear. CI integration waits until Task 2 makes the repository contract green.

- [ ] **Step 6: Review and commit**

Review the diff and focused test output. Commit message: `test: add repository integrity checks`.

### Task 2: Evidence-preserving repository hygiene

**Files:**
- Delete: `4-K_pneumoniae_/5-Post_Docking_Analysis`
- Delete: `audit/build-correctness.log`
- Delete: `audit/correctness-tests.log`
- Delete: `audit/correctness-tests.xml`
- Delete: `audit/docking-dry-run.log`
- Delete: `audit/wheel-verification.log`
- Delete: `autodock/prep_autodock.sh`
- Delete: `autodock/prep_ligands_custom.sh`
- Create: `audit/verification-history.md`
- Create: `specs/README.md`
- Modify: `audit/postdocking-audit.md`
- Modify: `docs/correctness-update.md`
- Modify: `docs/dockforge_migration_notes.md`
- Modify: `README.md`
- Modify: `.gitignore`
- Modify: `.github/workflows/correctness.yml`
- Move: root Sertraline `.docx` files to `closing_thesis/manuscripts/`
- Move: root `multiscale_overlapped_experimental_vs_triple_engine.{png,svg}` to `closing_thesis/figures/`
- Modify: `closing_thesis/README_closing_thesis_option2_serthcl.md`
- Modify: `closing_thesis/figures/FIGURE_INDEX.md`

**Interfaces:**
- Consumes: historical verification artifacts at the fixed base.
- Produces: a curated verification-history document with original blob hashes/recovery commands and a repository that passes the checker.

- [ ] **Step 1: Capture deletion evidence before editing**

Run `git ls-files -s` and `git hash-object` for every removal/move. Record original paths, blob hashes, commands, dates/outcomes, and the fixed-base recovery command in `audit/verification-history.md`.

- [ ] **Step 2: Apply only evidence-backed removals and moves**

Remove the unreferenced dangling private symlink and generated logs/XML. Retire the two broken legacy shell helpers: repository/package/reference searches must confirm that no active caller remains, and migration notes must route users to `prep_autodock_enhanced.sh` or the Python workflow. Keep audit reports, probes, JSON results, manuscripts, and figures. Move research assets without modifying their contents.

- [ ] **Step 3: Repair historical links and labels**

Replace links to deleted generated artifacts with the curated history. Change the missing probe-artifact link in `audit/postdocking-audit.md` to plain historical text that states the artifact was not tracked. Add an explicit historical-snapshot index for every numbered spec directory.

- [ ] **Step 4: Tighten ignore rules**

Ignore canonical runtime roots (`0-Input/` through `7-Reports/`, `.meta/`, `sessions/`), generated analysis outputs, local HPC profiles, and verification logs while keeping public templates and test fixtures trackable.

Add `python scripts/check_repository.py` to `.github/workflows/correctness.yml` only after the repository passes it.

- [ ] **Step 5: Verify**

Run:

```text
python scripts/check_repository.py
bash -n prep_autodock_enhanced.sh autodock/install_autodock.sh
git diff --check
git status --short
```

Confirm the checker passes and every deletion has its recorded evidence/recovery path.

- [ ] **Step 6: Review and commit**

Commit message: `chore: remove generated clutter and organize research assets`.

### Task 3: Deepen docking deployment planning

**Files:**
- Modify: `test/test_docking_correctness.py`
- Modify: `docking/deployment.py`

**Interfaces:**
- Consumes: existing `generate_condor_deployment(...)` and `generate_slurm_deployment(...)` signatures.
- Produces: the same generated files/paths plus `deployment_manifest["scheduler"]` for both adapters; shared private planning helpers remove duplicated implementation.

- [ ] **Step 1: Add a failing scheduler/common-manifest contract test**

```python
@pytest.mark.parametrize(
    ("factory", "scheduler", "options"),
    [(generate_slurm_deployment, "slurm", {"slurm_options": {"partition": "test"}}),
     (generate_condor_deployment, "condor", {})],
)
def test_deployment_manifests_identify_scheduler_and_common_context(project, factory, scheduler, options):
    root, pair = project
    result = factory(project_root=root, engines=["vina"], pairlist_rows=[pair],
                     runtime_by_engine={"vina": {"seed": 42}}, round_id="r1",
                     docking_mode="screen", **options)
    assert result["scheduler"] == scheduler
    assert result["generation_project_root"] == str(root.resolve())
    assert result["pair_count"] == 1
```

- [ ] **Step 2: Run the focused test and verify RED**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test/test_docking_correctness.py::test_deployment_manifests_identify_scheduler_and_common_context -q`
Expected: Slurm assertion fails because the field is absent.

- [ ] **Step 3: Extract shared private planning helpers**

Use a private context record and helpers for common project setup, runtime mapping, job planning, payload rows, and final manifest writing. Scheduler adapters keep their validation and script rendering. Do not introduce a scheduler class hierarchy: two adapters make the seam real, but their differences do not require a public abstraction.

- [ ] **Step 4: Run focused and deployment regression tests**

Run:

```text
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test/test_docking_correctness.py -q
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test/test_dockforge_smoke.py -q -k "deployment or hpc or condor or slurm"
```

- [ ] **Step 5: Compare generated artifacts**

Generate one local and one remote-root manifest for each scheduler in temporary directories. Diff key sets, job commands, mapped paths, submit modes, and executable permissions against the pre-refactor characterization; only the additive Slurm scheduler field may differ.

- [ ] **Step 6: Review and commit**

Commit message: `refactor: share docking deployment planning`.

### Task 4: Isolate CLI parser construction

**Files:**
- Create: `workflow/cli_parser.py`
- Modify: `workflow/cli.py`
- Modify: `test/test_repository_contracts.py`

**Interfaces:**
- Consumes: constants/choices passed from `workflow.cli` where circular imports would otherwise occur.
- Produces: `workflow.cli_parser.build_parser()` and compatibility export `workflow.cli.build_parser`.

- [ ] **Step 1: Write the failing parser-module contract**

```python
def test_cli_parser_module_accepts_representative_public_commands():
    from workflow.cli_parser import build_parser
    parser = build_parser()
    assert parser.parse_args(["workflow", "status", "--project-dir", "p"]).workflow_command == "status"
    assert parser.parse_args(["dock", "run", "--project-dir", "p"]).dock_command == "run"
    assert parser.parse_args(["analyze", "comparative", "--project-dir", "p"]).analyze_command == "comparative"
```

- [ ] **Step 2: Run the focused test and verify RED**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test/test_repository_contracts.py::test_cli_parser_module_accepts_representative_public_commands -q`
Expected: import fails because `workflow.cli_parser` does not exist.

- [ ] **Step 3: Move parser declaration behind one interface**

Split the 465-line builder into private group functions such as `_add_workflow_group`, `_add_pdb_group`, `_add_prep_group`, `_add_dock_group`, and `_add_analysis_group`. Reuse common argument builders; do not create duplicate parsers or change flag/default/help text. Keep a compatibility import in `workflow.cli`.

- [ ] **Step 4: Verify parser compatibility**

Compare parser actions recursively before/after using `(command path, option strings, required, nargs, default, choices, type name)`. Run representative `python main.py ... --help` commands and existing CLI smoke contracts.

- [ ] **Step 5: Run regression tests**

Run:

```text
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test/test_repository_contracts.py -q
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test/test_workflow_correctness.py -q
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test/test_dockforge_smoke.py -q -k "cli or workflow"
```

- [ ] **Step 6: Review and commit**

Commit message: `refactor: separate CLI declaration from dispatch`.

### Task 5: Consolidate and verify all documentation

**Files:**
- Create: `docs/README.md`
- Create: `docs/architecture.md`
- Create: `docs/testing.md`
- Create: `CONTRIBUTING.md`
- Modify: `README.md`
- Modify: `INSTALLATION_GUIDE.md`
- Modify: `DEPENDENCIES.md`
- Modify: `USAGE.md`
- Modify: `PDB_PREPARATION_USAGE.md`
- Modify: `PDBQT_PREPARATION_EXPLAINED.md`
- Modify: `HPC_DEPLOYMENT_GUIDE.md`
- Modify: `POST_DOCKING_ANALYSIS_GUIDE.md`
- Modify: `PLIP_INTERACTION_TYPES_GUIDE.md`
- Modify: `post_docking_analysis/README.md`
- Modify: `post_docking_analysis/VISUALIZATION_GUIDE.md`
- Modify: `post_docking_analysis/examples/README.md`
- Modify: `docs/dockforge_migration_notes.md`
- Modify: `docs/dockforge_release_candidate_summary_20260331.md`
- Modify: `docs/dockforge_release_checklist.md`
- Modify: `docs/dockforge_release_checklist_run_20260331.md`
- Modify: `CHANGELOG.md`
- Move: obsolete package-local design/status/usage documents to `post_docking_analysis/docs/_archive/`

**Interfaces:**
- Consumes: actual CLI help, package metadata, current source layout, historical correctness notes, and public-safe HPC templates.
- Produces: one indexed active documentation set with valid relative links and explicitly historical archives.

- [ ] **Step 1: Inventory every tracked Markdown file**

For each document, record `active`, `historical`, or `specialized research bundle`. Check its commands, links, version claims, output paths, dependency claims, and whether another active guide owns the same instruction.

- [ ] **Step 2: Rewrite the active landing/index documents**

Keep `README.md` compact and route detailed topics through `docs/README.md`. Remove unsupported production-readiness percentages and stale module names. State Python 3.10+, the legacy distribution/console naming, required external tools by workflow, and the current scientific limitations.

- [ ] **Step 3: Rewrite focused guides from live interfaces**

Use `python main.py --help` and nested help output as the source for examples. Cover preparation, local docking, Slurm, NMRBox/HTCondor template use, analysis modes, score directions, native-frame RMSD, `not_evaluable`, provenance/coverage, optional tools, migration, and clean test/build gates. Do not include private connection details.

- [ ] **Step 4: Archive obsolete documents without changing their historical claims**

Move the decision tree, pairlist integration summary, redesign proposal, simplified implementation status, and obsolete package-local usage guide under `post_docking_analysis/docs/_archive/`. Add an archive README stating that current usage lives in `POST_DOCKING_ANALYSIS_GUIDE.md` and that archived files describe earlier states.

- [ ] **Step 5: Verify all documentation and examples**

Run:

```text
python scripts/check_repository.py
python main.py --help
python main.py workflow --help
python main.py prep --help
python main.py dock --help
python main.py analyze --help
python -m post_docking_analysis --help
git diff --check
```

Confirm every tracked Markdown relative link resolves. Confirm every numbered spec and dated release record is labeled historical through its index/header without falsifying the original record.

- [ ] **Step 6: Review and commit**

Commit message: `docs: consolidate OmniDock guidance and history`.

### Task 6: Final distribution, regression, and independent review

**Files:**
- Modify only files required by concrete final-gate or review findings.

**Interfaces:**
- Consumes: the complete branch diff from the fixed base.
- Produces: verified commits and a reviewable draft pull request; no merge.

- [ ] **Step 1: Run repository/static gates**

```text
python scripts/check_repository.py
python -m compileall -q .
git diff --check 8d4434d...HEAD
```

Run `bash -n` for every tracked `.sh` file on a host with Bash.

- [ ] **Step 2: Run full regression**

Run `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest test -q` in a clean environment with RDKit available. On a host where Application Control blocks RDKit, record the exact collection error and run the non-RDKit suite without changing tests or converting failures to skips.

- [ ] **Step 3: Build and inspect the distribution**

```text
python -m build
python scripts/check_wheel.py
```

Install the wheel into a clean virtual environment and invoke/import all console entrypoints and package resources outside the source checkout.

- [ ] **Step 4: Measure before/after**

Record tracked file count/bytes, generated-artifact count/bytes, broken-link count, active-versus-historical document counts, and largest affected function sizes. Explain organization-only moves separately from actual footprint reduction.

- [ ] **Step 5: Run independent Standards and Spec reviews**

Use `git diff 8d4434d...HEAD` and `git log 8d4434d..HEAD --oneline`. Standards sources are `CONTRIBUTING.md`, the repository checker, CI, and the Fowler smell baseline. Spec source is `docs/superpowers/specs/2026-09-14-repository-cleanup-design.md`. Fix every Critical/Important or hard-spec finding and rerun gates.

- [ ] **Step 6: Push and open a draft pull request**

Push `codex/cleanup-simplify-docs`, open a draft PR, and let Linux/Windows Python 3.10/3.12 CI run. The PR description must list actual removals/moves/refactors, before/after evidence, local limitations, CI outcomes, and sequencing with the separate scientific-integrity PR. Do not merge.
