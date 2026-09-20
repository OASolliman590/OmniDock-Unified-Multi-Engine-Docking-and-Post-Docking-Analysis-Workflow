# Architecture

## Command declaration and dispatch

`main.py` is the public compatibility shim: it normalizes legacy arguments and
delegates to `workflow.cli.main`. `workflow/cli_parser.py` declares the unified
`workflow`, `pdb`, `prep`, `dock`, and `analyze` parser tree;
`workflow/cli.py` dispatches parsed commands. The workflow layer owns state
updates and adapts older implementations in the root modules, `docking/`, and
`post_docking_analysis/`. Package metadata retains legacy console names;
`pdb-wizard-workflow` points to `workflow.cli:main` while `pdb-prepare-wizard`
points to `main:main`. Treat parser help as the public command contract and the
dispatched modules as implementation details.

## State and artifact flow

```text
raw structures
  -> preparation records + prepared structures
  -> project_manifest.json + pairlist.csv
  -> per-engine configs, logs, poses, scores, completion evidence
  -> normalized score/geometry artifacts
  -> analysis sessions, reports, visualizations, interaction evidence
```

`workflow/state.py` records workflow context; `docking/project_layout.py` owns
canonical and compatibility paths; `docking/` owns command generation,
execution, deployment bundles, and output validation. Resume decisions bind
input content, executable identity, protocol, and validated output—not merely an
existing filename. Failed attempts do not inherit stale outputs as success.

The canonical project places engine artifacts below `engines/<engine>/` and
human-facing analysis under `5-Analysis/`, with numbered working, visualization,
and report directories created by the layout module. The `docking_legacy`
profile and `5-Post_Docking_Analysis` compatibility shim are migration seams,
not new output authorities.

## Boundaries

- Preparation owns chemical graph acquisition, atom mapping, protonation/tool
  selection, format export, and fail-closed QC records.
- Docking owns boxes, engine-specific parameters, execution identity, scheduler
  manifests, attempts, and output validation.
- Deployment generates portable artifacts locally; Slurm or HTCondor executes
  them remotely through public-safe profiles.
- Analysis owns engine detection, direction-aware score parsing/normalization,
  geometry evaluation, coverage, interactions, visualizations, and reports.

Adapters at these boundaries pass files plus explicit provenance. A consumer
must not infer chemical identity from coordinates, success from a filename, or
scientific validation from a generated chart.

## Public compatibility seams

Legacy console entrypoints, older raw GNINA layouts, engine-specific adapters,
the `docking_legacy` profile, and the post-docking compatibility path remain for
migration. New integrations should use `main.py`, the project manifest,
`pairlist.csv`, engine output contracts, and the numbered analysis topology.
See [migration notes](dockforge_migration_notes.md).

## Scientific ownership and failure semantics

The preparation record owns authoritative graph/mapping provenance. Docking
owns whether an engine attempt actually completed. Analysis owns whether an
available result is comparable. These states are distinct:

- a failed or unverifiable preparation/docking step fails closed;
- an unavailable scientific comparison is `not_evaluable` with a reason;
- omitted or unevaluable cases reduce reported coverage;
- a validated software gate only covers its named checks.

Energy scores sort ascending. GNINA CNNscore/CNNaffinity and consensus scores
sort descending. Cross-engine consensus normalizes by engine/scoring function
and target/site; raw unlike scores are not interchangeable. Redocking RMSD uses
the common receptor frame without ligand-on-reference fitting.

The scientific behavior is summarized in
[correctness-update.md](correctness-update.md). That file is not a substitute
for target-specific scientific review or a versioned benchmark.
