# DockForge migration and compatibility notes

This is a compatibility record, not the current usage manual. Follow the
[documentation index](README.md), [CLI recipes](../USAGE.md), and live
`python main.py ... --help` output for current operation.

## Names and entrypoints

The project is documented as OmniDock/DockForge while the Python distribution
and console entrypoints retain legacy `pdb-prepare-wizard` and `pdb-wizard-*`
names. In particular, `pdb-wizard-workflow` is declared as
`workflow.cli:main`. Existing installed commands remain compatibility surfaces;
new repository examples use `python main.py`.

## Layout seams

The canonical project uses `project_manifest.json`, `pairlist.csv`, per-engine
directories, numbered working/analysis/visualization/report roots, and workflow
state. `docking_legacy`, raw GNINA adapters, and
`5-Post_Docking_Analysis` are compatibility paths. New analysis writes to the
canonical `5-Analysis` root; on systems without symlink support the old path may
be a directory containing a routing notice.

Do not move or rename generated data solely to imitate an older tree. Initialize
or stage a canonical project, retain source/provenance records, and let the
layout/adapters resolve supported inputs.

## Behavioral compatibility changes

Current correctness checks intentionally reject or mark unevaluable cases that
older flows accepted permissively. Resume binds input contents, executable
identity, protocol, and validated outputs. Preparation requires authoritative
chemistry/mappings. Cross-engine comparison preserves score directions and
normalizes like evidence. RMSD uses the native receptor frame. Missing mappings,
poses, provenance, or comparisons remain `not_evaluable` and reduce coverage.

These changes can make an old project incomplete or unevaluable rather than
reproducing an unsupported result. Preserve the original project, restage from
authoritative inputs, rerun affected stages, and compare provenance. See the
[correctness update](correctness-update.md) and
[post-docking guide](../POST_DOCKING_ANALYSIS_GUIDE.md).

## Historical documents

Numbered specs, audit records, dated release checklists, and
[archived post-docking documents](../post_docking_analysis/docs/_archive/README.md)
describe earlier states. They remain useful migration evidence but do not
override active guides or current code.
