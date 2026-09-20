# Contributing

Use Python 3.10+ and an isolated environment. Install development and chemistry
dependencies with:

```bash
python -m pip install -e '.[dev,chemistry]'
```

Read the [architecture](docs/architecture.md) before changing module boundaries
and follow the complete [testing/evidence guide](docs/testing.md). Set the POSIX
or PowerShell CI-parity environment shown there before running pytest. Then run
the focused tests for the change plus the documented full gate, including:

```bash
python -m pytest test -q
python scripts/check_repository.py
python -m build
python scripts/check_wheel.py
git diff --check
```

Run `bash -n` on changed shell scripts in Bash/WSL/Linux. Do not claim current
CI, engine, GPU, scheduler, or scientific validation from a historical log or a
dry run.

## Compatibility boundaries

Preserve documented console entrypoints, project manifests, `pairlist.csv`,
engine output contracts, workflow state, and migration seams unless the change
explicitly includes a compatibility plan. New behavior should enter through the
unified `main.py` command tree and canonical project topology. Update active
guides and live help together; label dated design/release records as historical
rather than rewriting their claims.

## Scientific boundaries

Scientific changes require explicit provenance and fail-closed behavior. Do not
infer graph identity from coordinates, compare unlike raw scoring functions,
hide ligand displacement through superposition, convert missing evidence into a
pass, or omit `not_evaluable`/coverage accounting. Document engine score
direction, reference frame, mappings, exclusions, and failure semantics.

Regression tests are not prospective validation. Do not add claims about
affinity accuracy, enrichment, efficacy, production readiness, or clinical use
without versioned evidence and a separately reviewed scientific scope. Preserve
the limitations in [docs/correctness-update.md](docs/correctness-update.md).

## Repository hygiene and security

- keep generated data, caches, logs, environments, and build artifacts out of
  version control;
- use public-safe profile templates and never commit credentials, private hosts,
  usernames, account IDs, keys, or personal absolute paths;
- keep links relative and run the repository checker;
- include exact commands, environment, commit/tree, skips, and coverage when
  recording evidence;
- keep pull requests focused and explain behavior, compatibility, scientific
  implications, tests, and remaining limitations.
