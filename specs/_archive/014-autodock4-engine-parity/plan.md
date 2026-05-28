# Implementation Plan

1. Extend engine typing/registry and workflow target enums with `autodock4`.
2. Add AutoDock4 parser support (`.dlg`) and new runner implementation.
3. Thread AD4 args through docking CLI and workflow CLI argument forwarders.
4. Update interactive shell engine choices/defaults and workflow-init defaults.
5. Add hybrid ligand materialization for `ligands_engine_specific/autodock4`.
6. Extend built-in HPC profiles with AD4 runtime/slurm defaults.
7. Validate with `py_compile`, CLI help checks, and an AD4 dry-run smoke project.
