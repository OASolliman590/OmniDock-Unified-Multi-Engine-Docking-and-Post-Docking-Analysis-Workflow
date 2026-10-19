"""Deprecated entry points, kept for one release (Spec 038 R1).

Every function here does one of two things:
- prints a one-line deprecation notice naming the equivalent `main.py` command and
  forwards the request to the workflow CLI, or
- refuses with a clear message (exit code 2) when the legacy request has no
  `main.py` equivalent.

No scientific code lives here. Forwarded commands run through `workflow.cli`, the
same path `python main.py ...` uses.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Sequence

REFUSED = 2

# Flags accepted by the legacy `pdb-wizard-dock` (docking.cli.dock_main) that the
# `main.py dock run` parser does not accept. They are refused up front with a message.
LEGACY_ONLY_DOCK_FLAGS = (
    "--autodocktools-python",
    "--autodocktools-prepare-gpf4",
    "--autodocktools-prepare-dpf4",
    "--ligand-qc-allow-brenk",
    "--ligand-qc-allow-pains",
    "--ligand-qc-allow-reactive",
    "--ligand-qc-max-formal-charge-abs",
    "--ligand-qc-max-lipinski-violations",
    "--ligand-qc-max-logp",
    "--ligand-qc-max-molecular-weight",
    "--ligand-qc-max-rotatable-bonds",
    "--ligand-qc-max-tpsa",
    "--ligand-qc-min-heavy-atom-count",
    "--no-ligand-admet-filters",
)


def _argv(argv: Optional[Sequence[str]]) -> List[str]:
    return list(sys.argv[1:] if argv is None else argv)


def _notice(legacy: str, replacement: str) -> None:
    print(
        f"DEPRECATED: {legacy} will be removed in the next release. "
        f"Use: python main.py {replacement}",
        file=sys.stderr,
    )


def _refuse(legacy: str, reason: str, replacement: str) -> int:
    print(
        f"ERROR: {legacy} cannot forward this request: {reason}. "
        f"Use: python main.py {replacement}",
        file=sys.stderr,
    )
    return REFUSED


def _forward(main_argv: List[str]) -> int:
    from workflow.cli import main as workflow_main

    return workflow_main(main_argv)


def _parser(prog: str) -> argparse.ArgumentParser:
    return argparse.ArgumentParser(prog=prog, add_help=True)


# ---------------------------------------------------------------- pdb-wizard-legacy-interactive
def legacy_interactive(argv: Optional[Sequence[str]] = None) -> int:
    label = "pdb-wizard-legacy-interactive (interactive_pipeline.py)"
    args = _argv(argv)
    replacement = "pdb legacy-interactive [-o OUTPUT]"
    if args:
        return _refuse(label, f"takes no arguments (got: {' '.join(args)})", replacement)
    _notice(label, replacement)
    return _forward(["pdb", "legacy-interactive"])


# ---------------------------------------------------------------- pdb-wizard-cli (cli_pipeline.py)
def legacy_cli(argv: Optional[Sequence[str]] = None) -> int:
    label = "pdb-wizard-cli (cli_pipeline.py)"
    replacement = "pdb run -p <pdbs> -o <output> [-c <config>]"
    parser = _parser("pdb-wizard-cli")
    parser.add_argument("-p", "--pdbs", type=str)
    parser.add_argument("-o", "--output", type=str, default="pipeline_output")
    parser.add_argument("-c", "--config", type=str)
    parser.add_argument("--create-config", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(_argv(argv))
    if args.create_config:
        return _refuse(label, "--create-config has no main.py equivalent", replacement)
    if not args.pdbs:
        return _refuse(label, "-p/--pdbs is required", replacement)
    _notice(label, replacement + ("  (-v has no effect and is not forwarded)" if args.verbose else ""))
    forward = ["pdb", "run", "-p", args.pdbs, "-o", args.output]
    if args.config:
        forward += ["-c", args.config]
    return _forward(forward)


# ---------------------------------------------------------------- pdb-wizard-batch (batch_pdb_preparation.py)
def legacy_batch(argv: Optional[Sequence[str]] = None) -> int:
    label = "pdb-wizard-batch (batch_pdb_preparation.py)"
    replacement = "pdb batch -c <config> [-o <output>]"
    parser = _parser("pdb-wizard-batch")
    parser.add_argument("-c", "--config", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, default="batch_docking_preparation")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(_argv(argv))
    _notice(label, replacement + ("  (-v has no effect and is not forwarded)" if args.verbose else ""))
    return _forward(["pdb", "batch", "-c", args.config, "-o", args.output])


# ---------------------------------------------------------------- autodock_preparation.py CLI
def legacy_autodock_prepare(argv: Optional[Sequence[str]] = None) -> int:
    label = "autodock_preparation.py"
    replacement = "pdb prepare-both --receptors-input <dir> --ligands-input <dir> --force-field <FF> --ph <pH> ..."
    parser = _parser("autodock_preparation.py")
    parser.add_argument("--config")
    parser.add_argument("--ligands-input")
    parser.add_argument("--receptors-input")
    parser.add_argument("--ligands-output")
    parser.add_argument("--receptors-output")
    parser.add_argument("--force-field", choices=["AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON"], default=None)
    parser.add_argument("--ph", type=float, default=None)
    parser.add_argument("--protonation-policy", choices=["ph_model", "explicit_state", "as_input"], default="ph_model")
    parser.add_argument("--protonation-state-map", default="")
    parser.add_argument("--ligand-profile", default="engine_aware_full")
    parser.add_argument("--selected-engines", default="")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--create-config", action="store_true")
    args = parser.parse_args(_argv(argv))

    if args.create_config:
        return _refuse(label, "--create-config has no main.py equivalent", replacement)
    if args.config:
        return _refuse(label, "--config (JSON preparation config) has no main.py equivalent; pass the same values as flags", replacement)
    for input_path in (args.ligands_input, args.receptors_input):
        if input_path and Path(input_path).is_file():
            return _refuse(label, f"file input {input_path!r} is not supported; put it in a directory", replacement)

    _notice(label, replacement)
    forward = [
        "pdb", "prepare-both",
        "--receptors-input", args.receptors_input or "./receptors_raw",
        "--ligands-input", args.ligands_input or "./ligands_raw",
        "--receptors-output", args.receptors_output or "./receptors_prep",
        "--ligands-output", args.ligands_output or "./ligands_prep",
        "--protonation-policy", args.protonation_policy,
        "--ligand-backend", args.ligand_profile,
    ]
    if args.force_field is not None:
        forward += ["--force-field", args.force_field]
    if args.ph is not None:
        forward += ["--ph", str(args.ph)]
    if args.protonation_state_map:
        forward += ["--protonation-state-map", str(Path(args.protonation_state_map).expanduser().resolve())]
    if args.selected_engines:
        forward += ["--selected-engines", args.selected_engines]
    if args.force:
        forward.append("--force")
    return _forward(forward)


# ---------------------------------------------------------------- pdb-wizard-prepare-docking
def prepare_docking(argv: Optional[Sequence[str]] = None) -> int:
    label = "pdb-wizard-prepare-docking"
    _notice(label, "prep project ...  (same flags)")
    return _forward(["prep", "project", *_argv(argv)])


# ---------------------------------------------------------------- pdb-wizard-dock
def dock(argv: Optional[Sequence[str]] = None) -> int:
    label = "pdb-wizard-dock"
    replacement = "dock run ..."
    args = _argv(argv)
    for token in args:
        for flag in LEGACY_ONLY_DOCK_FLAGS:
            if token == flag or token.startswith(flag + "="):
                return _refuse(label, f"{flag} is not accepted by `dock run`", replacement)
    _notice(label, replacement)
    return _forward(["dock", "run", *args])
