# OmniDock documentation

This page is the canonical navigation index. Active guides describe the current
checkout. Historical records preserve what was proposed or observed at a dated
point; they are not current guarantees.

## Active guides

- [Project landing page](../README.md)
- [Installation](../INSTALLATION_GUIDE.md)
- [Dependencies](../DEPENDENCIES.md)
- [CLI recipes](../USAGE.md)
- [Architecture](architecture.md)
- [Testing and current evidence](testing.md)
- [PDB and receptor/ligand preparation](../PDB_PREPARATION_USAGE.md)
- [PDBQT preparation and mapping](../PDBQT_PREPARATION_EXPLAINED.md)
- [HPC deployment](../HPC_DEPLOYMENT_GUIDE.md)
- [Post-docking analysis](../POST_DOCKING_ANALYSIS_GUIDE.md)
- [Interaction types](../PLIP_INTERACTION_TYPES_GUIDE.md)
- [Contributing](../CONTRIBUTING.md)

## Current correctness and compatibility records

- [Scientific correctness update](correctness-update.md) describes implemented
  behavior, known limitations, and the dated local verification record.
- [DockForge migration notes](dockforge_migration_notes.md) record compatibility
  seams and route users to the current guides.
- [Changelog](../CHANGELOG.md) records repository changes; old roadmap statements
  are historical, not promises.

## Historical material

- [Audit index](../audit/AUDIT.md) and
  [verification history](../audit/verification-history.md) preserve findings and
  dated evidence.
- [Numbered specifications](../specs/README.md) are design snapshots. Their
  plans/tasks are not an alternative current manual.
- Dated release candidate/checklist records are explicitly snapshot-labeled:
  [summary](dockforge_release_candidate_summary_20260331.md),
  [checklist](dockforge_release_checklist.md), and
  [checklist run](dockforge_release_checklist_run_20260331.md).
- [Archived post-docking documents](../post_docking_analysis/docs/_archive/README.md)
  describe superseded package states.

## Specialized research bundles

- [Closing-thesis research bundle](../closing_thesis/README_closing_thesis_option2_serthcl.md)
  and its [figure index](../closing_thesis/figures/FIGURE_INDEX.md) are
  project-specific research material, not general OmniDock documentation.
- [Post-docking package notes](../post_docking_analysis/README.md),
  [visualization notes](../post_docking_analysis/VISUALIZATION_GUIDE.md), and
  [examples](../post_docking_analysis/examples/README.md) are specialized entry
  points subordinate to the top-level analysis guide.

When a historical file conflicts with live CLI help, package metadata, or an
active guide, use the live surface and active guide and open an issue for the
contradiction.
