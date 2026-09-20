# Post-docking examples

These scripts illustrate lower-level Python/package integration and may require
optional dependencies or local input data. For supported CLI recipes, use the
top-level [CLI guide](../../USAGE.md) and
[post-docking guide](../../POST_DOCKING_ANALYSIS_GUIDE.md).

Before running an example:

1. inspect its input/output paths and configuration;
2. use a disposable project or explicit output directory;
3. confirm required engine, chemistry, interaction, and visualization tools;
4. retain provenance, coverage, and `not_evaluable` reasons;
5. do not treat demonstration output as a scientific benchmark.

Start with `python main.py analyze --help` for canonical projects or
`python -m post_docking_analysis --help` for the direct package interface.
