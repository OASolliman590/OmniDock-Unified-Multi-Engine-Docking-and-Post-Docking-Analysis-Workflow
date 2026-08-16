"""Catalog of runnable targets and their web form schemas.

The vocabulary is imported from ``workflow.interactive`` at runtime rather
than restated here, so adding a target to the terminal UI surfaces it in
the web UI with no change to this module (FR-007, SC-002).

Every generated argv is checked against the real ``workflow.cli`` parser
before a job is created, so a template that drifts from the CLI fails
loudly here instead of at run time.
"""

from __future__ import annotations

import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

# Imported as a module, not as `from ... import ANALYSIS_LABELS`, so the
# constants are read at call time. A `from` import binds the value once at
# import and the web UI would then keep serving a stale vocabulary after
# the terminal UI's lists change (SC-002).
from workflow import interactive

# Interaction aliases that the pipeline reroutes to the clean contract
# (spec 029 FR-002). The UI states this rather than silently substituting.
CLEAN_ROUTED_TARGETS = {
    "analyze.interactions.pandamap",
    "analyze.interactions.prolif",
    "analyze.interactions.ligplot",
}
CLEAN_TARGET = "analyze.interactions.clean"


class TargetError(Exception):
    """Raised for an unknown target or invalid form values."""


@dataclass
class FormField:
    name: str
    label: str
    type: str  # text | path | select | multiselect | bool | int | float
    required: bool = False
    default: object = None
    choices: list[dict] | None = None
    help: str = ""
    flag: str = ""  # CLI flag; defaults to --<name with underscores as hyphens>

    def cli_flag(self) -> str:
        return self.flag or "--" + self.name.replace("_", "-")


@dataclass
class LaunchTarget:
    key: str
    label: str
    group: str
    cli_path: list[str]
    fields: list[FormField] = field(default_factory=list)
    routes_to: str | None = None

    def to_dict(self) -> dict:
        return {
            "key": self.key,
            "label": self.label,
            "group": self.group,
            "routes_to": self.routes_to,
            "fields": [asdict(f) for f in self.fields],
        }


def _engine_choices() -> list[dict]:
    """Engine options, straight from the terminal UI's list."""
    return [{"value": value, "label": label} for value, label in interactive.ENGINE_CHOICES if value != "all"]


def _cli_path_for(key: str) -> list[str]:
    """Translate an analysis-target token into CLI subcommand tokens.

    ``analyze.stage.structure_quality`` -> ``["analyze", "stage", "structure-quality"]``.
    Underscores become hyphens because argparse subcommands are declared
    that way (``favorite-engine``, ``structure-quality``).
    """
    return [token.replace("_", "-") for token in key.split(".")]


def _common_analysis_fields() -> list[FormField]:
    """Fields shared by targets built on _add_analysis_common_args."""
    return [
        FormField("output", "Output directory", "path", help="Defaults to the project's analysis root."),
        FormField("config_file", "Config file", "path", help="Optional post-docking analysis YAML/JSON."),
        FormField("engine", "Engine", "select", choices=_engine_choices(),
                  help="Engine for stage-level analysis in canonical projects."),
        FormField("favorite_engine", "Favorite engine", "select", choices=_engine_choices()),
        FormField("complex_query", "Complex filter", "text",
                  help="e.g. protein=2FVD;ligand=Sorafenib. Use 0/skip/all to disable."),
        FormField("pairlist", "Pairlist CSV", "path"),
        FormField("ligplus_root", "LigPlus root", "path"),
        FormField("enable_poseview", "Enable PoseView stage", "bool", default=False),
    ]


def _build_analysis_targets() -> dict[str, LaunchTarget]:
    """One target per ANALYSIS_LABELS entry -- never a hardcoded list."""
    targets: dict[str, LaunchTarget] = {}

    for key, label in interactive.ANALYSIS_LABELS.items():
        cli_path = _cli_path_for(key)
        fields = _common_analysis_fields()

        if key == "analyze.favorite_engine":
            fields = [
                FormField("favorite_engine", "Favorite engine", "select", required=True,
                          choices=_engine_choices()),
                FormField("output", "Output directory", "path"),
                FormField("config_file", "Config file", "path"),
                FormField("complex_query", "Complex filter", "text"),
            ]
        elif key == "analyze.comparative":
            fields = [
                FormField("output", "Output directory", "path"),
                FormField("config_file", "Config file", "path"),
                FormField("consensus_mode", "Consensus mode", "select",
                          default="dockbox_geometric",
                          choices=[{"value": v, "label": v} for v in
                                   ["dockbox_geometric", "weighted_hybrid",
                                    "strict_consensus", "favorite_guardrails"]]),
                FormField("rescoring_scope", "Rescoring scope", "select",
                          default="top_n_per_protein",
                          choices=[{"value": v, "label": v} for v in
                                   ["top_n_per_protein", "top_n_global"]]),
                FormField("rescoring_top_n", "Rescoring top N", "int", default=3),
                FormField("complex_query", "Complex filter", "text"),
            ]
        elif key == CLEAN_TARGET:
            fields = [
                FormField("output", "Output directory", "path"),
                FormField("config_file", "Config file", "path"),
                FormField("dataset_root", "Dataset root", "path",
                          help="Contains pose_ensemble_appraisal and chemistry-fix scripts."),
            ]

        targets[key] = LaunchTarget(
            key=key,
            label=label,
            group="analysis",
            cli_path=cli_path,
            fields=fields,
            routes_to=CLEAN_TARGET if key in CLEAN_ROUTED_TARGETS else None,
        )

    return targets


def _build_pipeline_targets() -> dict[str, LaunchTarget]:
    """Preparation and docking commands, which are not in ANALYSIS_LABELS."""
    engines = _engine_choices()
    return {
        # Note: `prep pairlist` has no --output; it writes into the project.
        "prep.pairlist": LaunchTarget(
            key="prep.pairlist", label="Build pairlist", group="preparation",
            cli_path=["prep", "pairlist"],
            fields=[
                FormField("mode", "Pairlist mode", "select",
                          choices=[{"value": v, "label": v} for v in
                                   ["cocrystal_only", "cocrystal_plus_all",
                                    "cocrystal_plus_nonreference",
                                    "curated_cartesian", "curated_per_protein"]]),
                FormField("prepared_proteins", "Prepared proteins dir", "path"),
                FormField("prepared_ligands", "Prepared ligands dir", "path"),
                FormField("default_site_id", "Default site id", "text", default="site_1"),
                FormField("default_box_size", "Default box size", "float", default=20.0),
                FormField("freeze", "Freeze into pairlist.csv", "bool", default=False),
            ],
        ),
        "prep.project": LaunchTarget(
            key="prep.project", label="Prepare docking folders", group="preparation",
            cli_path=["prep", "project"],
            fields=[
                FormField("engines", "Engines", "multiselect", choices=engines),
                FormField("output", "Output directory", "path"),
            ],
        ),
        "dock.run": LaunchTarget(
            key="dock.run", label="Run docking", group="docking",
            cli_path=["dock", "run"],
            fields=[
                FormField("engines", "Engines", "multiselect", choices=engines),
                FormField("mode", "Mode", "select", default="standard",
                          choices=[{"value": "standard", "label": "standard"},
                                   {"value": "score-only", "label": "score-only"}]),
                FormField("favorite_engine", "Favorite engine", "select", choices=engines),
                FormField("skip_completed", "Skip completed rows", "bool", default=False),
            ],
        ),
        "dock.dry-run": LaunchTarget(
            key="dock.dry-run", label="Docking dry run (generate commands only)", group="docking",
            cli_path=["dock", "dry-run"],
            fields=[
                FormField("engines", "Engines", "multiselect", choices=engines),
                FormField("mode", "Mode", "select", default="standard",
                          choices=[{"value": "standard", "label": "standard"},
                                   {"value": "score-only", "label": "score-only"}]),
            ],
        ),
        "dock.deploy": LaunchTarget(
            key="dock.deploy", label="HPC: generate deployment assets", group="hpc",
            cli_path=["dock", "deploy"],
            fields=[
                FormField("hpc_profile", "HPC profile", "text", required=True),
                FormField("engines", "Engines", "multiselect", choices=engines),
                FormField("remote_project_dir", "Remote project dir", "text"),
            ],
        ),
        "dock.sync": LaunchTarget(
            key="dock.sync", label="HPC: sync project to cluster", group="hpc",
            cli_path=["dock", "sync"],
            fields=[
                FormField("hpc_profile", "HPC profile", "text", required=True),
                FormField("remote_project_dir", "Remote project dir", "text"),
                FormField("dry_run", "Dry run", "bool", default=False),
            ],
        ),
        "dock.submit": LaunchTarget(
            key="dock.submit", label="HPC: submit synced deployment", group="hpc",
            cli_path=["dock", "submit"],
            fields=[
                FormField("hpc_profile", "HPC profile", "text", required=True),
                FormField("engines", "Engines", "multiselect", choices=engines),
                FormField("mode", "Mode", "select",
                          choices=[{"value": "screen", "label": "screen"},
                                   {"value": "exhaustive", "label": "exhaustive"}]),
                FormField("round", "Deployment round", "text"),
            ],
        ),
        "workflow.init": LaunchTarget(
            key="workflow.init", label="Initialize workflow state", group="workflow",
            cli_path=["workflow", "init"],
            fields=[
                FormField("engines", "Engines", "multiselect", choices=engines),
                FormField("project_name", "Project name", "text"),
                FormField("favorite_engine", "Favorite engine", "select", choices=engines),
            ],
        ),
    }


def all_targets() -> dict[str, LaunchTarget]:
    catalog = _build_analysis_targets()
    catalog.update(_build_pipeline_targets())
    return catalog


def get_target(key: str) -> LaunchTarget:
    try:
        return all_targets()[key]
    except KeyError:
        raise TargetError(f"Unknown target: {key}") from None


def form_schema(key: str, context: dict | None = None) -> dict:
    """Form schema with defaults pre-filled from the project's context."""
    target = get_target(key)
    schema = target.to_dict()
    ctx = context or {}

    for entry in schema["fields"]:
        name = entry["name"]
        if entry.get("default") not in (None, ""):
            continue
        if name == "favorite_engine" and ctx.get("favorite_engine"):
            entry["default"] = ctx["favorite_engine"]
        elif name == "engines" and ctx.get("selected_engines"):
            entry["default"] = list(ctx["selected_engines"])
        elif name == "engine" and ctx.get("favorite_engine"):
            entry["default"] = ctx["favorite_engine"]
        elif name == "output" and ctx.get("analysis_output_dir"):
            entry["default"] = ctx["analysis_output_dir"]

    return schema


# ------------------------------------------------------------------ argv


def validate(key: str, values: dict) -> dict[str, str]:
    """Return {field_name: message} for anything wrong. Empty means valid."""
    target = get_target(key)
    errors: dict[str, str] = {}
    valid_engines = {c["value"] for c in _engine_choices()}

    for f in target.fields:
        raw = values.get(f.name)
        provided = raw not in (None, "", [], False)

        if f.required and not provided:
            errors[f.name] = f"{f.label} is required."
            continue
        if not provided:
            continue

        if f.type == "path":
            if not Path(str(raw)).expanduser().exists():
                errors[f.name] = f"Path does not exist: {raw}"
        elif f.type in ("int", "float"):
            try:
                float(raw)
            except (TypeError, ValueError):
                errors[f.name] = f"{f.label} must be a number."
        elif f.type == "select" and f.choices:
            allowed = {c["value"] for c in f.choices}
            if str(raw) not in allowed:
                errors[f.name] = f"Unknown value '{raw}' for {f.label}."
        elif f.type == "multiselect":
            items = raw if isinstance(raw, list) else [raw]
            if f.name == "engines":
                unknown = [i for i in items if i not in valid_engines]
                if unknown:
                    errors[f.name] = f"Unknown engine(s): {', '.join(map(str, unknown))}"

    return errors


def build_argv(key: str, values: dict, project_path: str) -> list[str]:
    """Build the argv for a job. Always python -c main.py ... (FR-010)."""
    target = get_target(key)
    argv = [sys.executable, "main.py", *target.cli_path, "--project-dir", str(project_path)]

    for f in target.fields:
        raw = values.get(f.name)
        if raw in (None, "", [], False):
            continue

        if f.type == "bool":
            argv.append(f.cli_flag())
        elif f.type == "multiselect":
            items = raw if isinstance(raw, list) else [raw]
            argv += [f.cli_flag(), ",".join(str(i) for i in items)]
        elif f.name == "output":
            # argparse declares this as -o/--output
            argv += ["--output", str(raw)]
        else:
            argv += [f.cli_flag(), str(raw)]

    return argv


def assert_argv_parses(argv: list[str]) -> None:
    """Check the generated argv against the real CLI parser.

    A template that drifts from ``workflow/cli.py`` fails here rather than
    producing a job that dies on startup.
    """
    from workflow.cli import build_parser

    parser = build_parser()
    try:
        parser.parse_args(argv[2:])  # drop [python, main.py]
    except SystemExit as exc:  # argparse exits on error
        raise TargetError(
            f"Generated command is not accepted by the CLI parser: {' '.join(argv[2:])}"
        ) from exc
