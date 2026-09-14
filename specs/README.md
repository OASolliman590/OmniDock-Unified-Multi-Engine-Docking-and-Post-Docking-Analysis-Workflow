# Specs directory

This directory holds **historical implementation snapshots**. It is not the
current requirements set, not the current architecture description, and not the
active user documentation.

Treat every numbered subdirectory `001` through `029` and everything under
`_archive/` as a completed-cycle record of an earlier design. Claims inside
those files describe the state at the time they were written. They must not be
rewritten as current product requirements or as current scientific validation.

## Current work versus these snapshots

- Current repository-cleanup work lives under `docs/superpowers/` (plan, design
  specification, and progress notes).
- Active user documentation currently lives in the repository-root guides
  (`README.md`, `INSTALLATION_GUIDE.md`, `USAGE.md`, `PDB_PREPARATION_USAGE.md`,
  `PDBQT_PREPARATION_EXPLAINED.md`, `HPC_DEPLOYMENT_GUIDE.md`,
  `POST_DOCKING_ANALYSIS_GUIDE.md`, and related files) plus `docs/` historical
  notes. A later cleanup batch will add a documentation index; this file does
  not link to a path that does not exist yet.

## Numbered historical snapshots (`001`–`029`)

| Directory | Snapshot title (from the historical spec heading / directory name) |
| --- | --- |
| `001-unify-and-optimize` | Unified visualization overhaul |
| `002-pipeline-unification` | Pipeline unification and deduplication |
| `003-simplified-pipeline-dedup` | Simplified pipeline deduplication |
| `004-visualization-rmsd-overhaul` | Visualization + RMSD overhaul |
| `005-interactive-cli-prompts` | Interactive simplified CLI prompts |
| `006-run-log-error-remediation` | Simplified pipeline run-log error remediation |
| `007-multi-engine-docking-orchestration` | Multi-engine docking orchestration |
| `008-engine-aware-post-docking-analysis` | Engine-aware post-docking analysis |
| `009-docking-preparation-automation` | Docking preparation automation |
| `010-unified-workflow-orchestrator` | Unified workflow orchestrator |
| `011-docking-first-project-workflow` | Docking-first project workflow |
| `012-backlog-reconciliation` | Backlog reconciliation after workflow foundation |
| `013-post-docking-closeout` | Post-docking closeout |
| `014-multi-engine-hpc-docking-orchestration` | Multi-engine HPC docking orchestration |
| `015-ligand-normalization-and-raw-ligand-unification` | Ligand normalization and raw-ligand unification |
| `016-post-docking-ux-output-restructure` | Post-docking UX output restructure |
| `017-dockbox-consensus-maturation` | DockBox-style consensus and post-docking maturation |
| `018-dockforge-omnidock-platform` | Omni-DockForge platform |
| `019-top-pose-ligand-performance` | Top-pose ligand performance atlas |
| `020-post-docking-scientific-foundation` | Post-docking scientific foundation |
| `021-unified-output-topology` | Unified output topology |
| `022-engine-detection-and-solo-routing` | Engine detection and solo routing |
| `023-engine-scope-filter` | Engine scope filter |
| `024-artifact-dag-pipeline` | Artifact DAG pipeline |
| `025-visualization-report-suite` | Visualization report suite |
| `026-pipeline-correctness-fixes` | Pipeline correctness fixes |
| `027-consensus-validation-correctness` | Consensus validation correctness |
| `028-visualization-corrections` | Visualization corrections and scientific guardrails |
| `029-post-docking-remediation` | Post-docking analysis remediation and contract uniformity |

## `_archive/` historical snapshots

`_archive/` is also a historical snapshot area, not a current requirements tree.

| Directory | Snapshot title |
| --- | --- |
| `_archive/014-autodock4-engine-parity` | AutoDock4 engine parity and unified ligand prep |

No numbered or archived spec directory is omitted from this index.
