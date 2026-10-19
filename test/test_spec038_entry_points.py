"""Spec 038 R1, R2 and R5: one CLI, moved scientific code, one documented path.

- R1: each deprecated entry point forwards to its `main.py` equivalent or refuses clearly.
- R2: old import paths resolve to the moved functions (identity); nothing imports removed paths.
- R5: USAGE.md, README.md and the `main.py --help` examples use commands and flags the parser accepts.
- Packaging: `python -m build` plus scripts/check_wheel.py pass (build outputs are removed afterwards).
"""

from __future__ import annotations

import ast
import argparse
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
SEARCH_ROOTS = (REPO, REPO / "test")


# ---------------------------------------------------------------- helpers
def _parser() -> argparse.ArgumentParser:
    from workflow.cli import build_parser

    return build_parser()


def _subparser_action(parser):
    return next((a for a in parser._actions if isinstance(a, argparse._SubParsersAction)), None)


def _resolve(parser, tokens):
    """Walk `group command ...` through the subparsers; return (leaf parser, index after the path)."""
    node, index = parser, 0
    while index < len(tokens) and not tokens[index].startswith("-"):
        action = _subparser_action(node)
        if action is None or tokens[index] not in action.choices:
            break
        node = action.choices[tokens[index]]
        index += 1
    return node, index


def _check_command(tokens) -> list[str]:
    """Return a list of problems for one `main.py` argv (path, flags, required options)."""
    problems = []
    leaf, consumed = _resolve(_parser(), tokens)
    if consumed == 0:
        if tokens and all(t in ("-h", "--help") for t in tokens):
            return []  # `python main.py --help`
        return [f"no known command in {tokens!r}"]
    known = {opt for action in leaf._actions for opt in action.option_strings}
    for token in tokens[consumed:]:
        if token.startswith("--") and token.split("=", 1)[0] not in known:
            problems.append(f"unknown flag {token} for {' '.join(tokens[:consumed])}")
    for action in leaf._actions:
        if action.required and action.option_strings:
            if not any(opt in tokens for opt in action.option_strings):
                problems.append(f"missing required {action.option_strings[0]} for {' '.join(tokens[:consumed])}")
    return problems


def _documented_commands(path: Path) -> list[list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    commands = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("python main.py"):
            buf = line
            while buf.endswith("\\") and i + 1 < len(lines):
                i += 1
                buf = buf[:-1].rstrip() + " " + lines[i].strip()
            buf = re.sub(r"\s+#.*$", "", buf)
            tokens = [t.strip("[]") for t in shlex.split(buf)][2:]
            tokens = [t for t in tokens if t]
            commands.append(tokens)
        i += 1
    return commands


def _module_file(name: str) -> Path | None:
    rel = Path(*name.split("."))
    for root in SEARCH_ROOTS:
        for candidate in (root / f"{rel}.py", root / rel / "__init__.py"):
            if candidate.exists():
                return candidate
    return None


def _top_level_names(path: Path) -> tuple[set[str], bool]:
    """Names a module binds at top level, and whether it has a star import (names then unknown)."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names: set[str] = set()
    has_star = False
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, (ast.FunctionDef,)) and node.name == "__getattr__":
            has_star = True  # PEP 562 lazy names
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
        elif isinstance(node, ast.ImportFrom) and any(a.name == "*" for a in node.names):
            has_star = True
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                for sub in ast.walk(target):
                    if isinstance(sub, ast.Name):
                        names.add(sub.id)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                names.add((alias.asname or alias.name).split(".")[0])
        elif isinstance(node, ast.Try):
            for sub in ast.walk(node):
                if isinstance(sub, (ast.Import, ast.ImportFrom)):
                    for alias in sub.names:
                        names.add((alias.asname or alias.name).split(".")[0])
                elif isinstance(sub, ast.Name) and isinstance(sub.ctx, ast.Store):
                    names.add(sub.id)
    return names, has_star


class _Recorder:
    def __init__(self, code=0):
        self.calls = []
        self.code = code

    def __call__(self, argv):
        self.calls.append(list(argv))
        return self.code


@pytest.fixture
def forwarded(monkeypatch):
    recorder = _Recorder()
    import workflow.cli

    monkeypatch.setattr(workflow.cli, "main", recorder)
    return recorder


# ---------------------------------------------------------------- R1: shims forward or refuse
def test_legacy_cli_forwards_to_pdb_run(forwarded, capsys):
    from workflow.legacy_entrypoints import legacy_cli

    assert legacy_cli(["-p", "7CMD,6WX4", "-o", "out", "-c", "cfg.json", "-v"]) == 0
    assert forwarded.calls == [["pdb", "run", "-p", "7CMD,6WX4", "-o", "out", "-c", "cfg.json"]]
    err = capsys.readouterr().err
    assert "DEPRECATED" in err and "python main.py pdb run" in err and "has no effect" in err


def test_legacy_cli_refuses_create_config_and_missing_pdbs(forwarded, capsys):
    from workflow.legacy_entrypoints import REFUSED, legacy_cli

    assert legacy_cli(["--create-config"]) == REFUSED
    assert legacy_cli([]) == REFUSED
    assert forwarded.calls == []
    assert "ERROR" in capsys.readouterr().err


def test_legacy_batch_forwards_to_pdb_batch(forwarded):
    from workflow.legacy_entrypoints import legacy_batch

    assert legacy_batch(["-c", "pdb_batch_config.yaml"]) == 0
    assert forwarded.calls == [["pdb", "batch", "-c", "pdb_batch_config.yaml", "-o", "batch_docking_preparation"]]


def test_legacy_autodock_forwards_with_flag_mapping(forwarded, tmp_path):
    from workflow.legacy_entrypoints import legacy_autodock_prepare

    rc = legacy_autodock_prepare([
        "--receptors-input", "rin", "--ligands-input", "lin",
        "--receptors-output", "rout", "--ligands-output", "lout",
        "--force-field", "AMBER", "--ph", "7.4",
        "--ligand-profile", "openbabel_only", "--selected-engines", "vina", "--force",
    ])
    assert rc == 0
    assert forwarded.calls == [[
        "pdb", "prepare-both",
        "--receptors-input", "rin", "--ligands-input", "lin",
        "--receptors-output", "rout", "--ligands-output", "lout",
        "--protonation-policy", "ph_model", "--ligand-backend", "openbabel_only",
        "--force-field", "AMBER", "--ph", "7.4", "--selected-engines", "vina", "--force",
    ]]


def test_legacy_autodock_refuses_config_create_config_and_file_inputs(forwarded, tmp_path):
    from workflow.legacy_entrypoints import REFUSED, legacy_autodock_prepare

    raw_file = tmp_path / "ligand.sdf"
    raw_file.write_text("x", encoding="utf-8")
    assert legacy_autodock_prepare(["--config", "c.json"]) == REFUSED
    assert legacy_autodock_prepare(["--create-config"]) == REFUSED
    assert legacy_autodock_prepare(["--ligands-input", str(raw_file)]) == REFUSED
    assert forwarded.calls == []


def test_legacy_interactive_forwards_and_refuses_arguments(forwarded):
    from workflow.legacy_entrypoints import REFUSED, legacy_interactive

    assert legacy_interactive([]) == 0
    assert forwarded.calls == [["pdb", "legacy-interactive"]]
    assert legacy_interactive(["-o", "x"]) == REFUSED
    assert len(forwarded.calls) == 1


def test_prepare_docking_forwards_to_prep_project(forwarded):
    from workflow.legacy_entrypoints import prepare_docking

    argv = ["--prepared-proteins", "p", "--prepared-ligands", "l", "--excel", "e.xlsx", "--output", "o"]
    assert prepare_docking(argv) == 0
    assert forwarded.calls == [["prep", "project", *argv]]


def test_dock_forwards_to_dock_run_and_refuses_legacy_only_flags(forwarded, capsys):
    from workflow.legacy_entrypoints import REFUSED, dock

    assert dock(["--project-dir", "p", "--engines", "vina"]) == 0
    assert forwarded.calls == [["dock", "run", "--project-dir", "p", "--engines", "vina"]]
    assert dock(["--project-dir", "p", "--ligand-qc-max-tpsa=5"]) == REFUSED
    assert len(forwarded.calls) == 1
    assert "--ligand-qc-max-tpsa" in capsys.readouterr().err


@pytest.mark.parametrize(
    "legacy_flags",
    [
        ["--project-dir", "p", "--autodocktools-prepare-gpf4", "x"],
        ["--project-dir", "p", "--no-ligand-admet-filters"],
    ],
)
def test_dock_refuses_each_legacy_only_group(forwarded, legacy_flags):
    from workflow.legacy_entrypoints import REFUSED, dock

    assert dock(legacy_flags) == REFUSED
    assert forwarded.calls == []


def test_forwarded_argv_parses_with_the_real_parser(forwarded):
    """Every shim forward must be accepted by the workflow parser (no execution)."""
    from workflow.legacy_entrypoints import dock, legacy_batch, legacy_cli, prepare_docking

    legacy_cli(["-p", "7CMD", "-o", "out"])
    legacy_batch(["-c", "cfg.yaml"])
    prepare_docking(["--prepared-proteins", "p", "--prepared-ligands", "l", "--excel", "e", "--output", "o"])
    dock(["--project-dir", "p", "--engines", "vina", "--dry-run"])
    parser = _parser()
    for argv in forwarded.calls:
        parser.parse_args(argv)  # raises SystemExit on a bad forward


def test_legacy_script_entry_points_return_exit_codes(capsys):
    import cli_pipeline

    assert cli_pipeline.main(["--create-config"]) == 2
    assert "ERROR" in capsys.readouterr().err


def test_console_scripts_point_at_supported_or_shim_callables():
    setup_tree = ast.parse((REPO / "setup.py").read_text(encoding="utf-8"))
    scripts = None
    for node in ast.walk(setup_tree):
        if isinstance(node, ast.keyword) and node.arg == "entry_points":
            for key, value in zip(node.value.keys, node.value.values):
                if isinstance(key, ast.Constant) and key.value == "console_scripts":
                    scripts = [elt.value for elt in value.elts]
    assert scripts is not None
    table = dict(item.split("=", 1) for item in scripts)

    expected = {
        "pdb-prepare-wizard": "workflow.cli:main",
        "pdb-wizard-workflow": "workflow.cli:main",
        "pdb-wizard-interactive": "workflow.interactive:run_interactive_workflow",
        "pdb-wizard-legacy-interactive": "interactive_pipeline:main",
        "pdb-wizard-cli": "cli_pipeline:main",
        "pdb-wizard-batch": "batch_pdb_preparation:main",
        "pdb-wizard-prepare-docking": "workflow.legacy_entrypoints:prepare_docking",
        "pdb-wizard-dock": "workflow.legacy_entrypoints:dock",
        # Documented exception: owned by Spec 038 R3 / T004, not converted in this change.
        "post-docking-analysis": "post_docking_analysis.cli:main",
    }
    assert table == expected
    import importlib

    for name, target in table.items():
        module_name, attribute = target.split(":")
        assert callable(getattr(importlib.import_module(module_name), attribute)), name


def test_shim_module_is_deprecation_only_and_holds_no_scientific_code():
    source = (REPO / "workflow" / "legacy_entrypoints.py").read_text(encoding="utf-8")
    for forbidden in ("core_pipeline", "docking.preparation", "numpy", "Bio.PDB", "openpyxl"):
        assert forbidden not in source, forbidden


# ---------------------------------------------------------------- R2: moved code and import identity
def test_old_import_paths_resolve_to_moved_functions():
    import autodock_preparation
    import cli_pipeline
    import core_pipeline
    import interactive_pipeline
    from docking.preparation import binding_site_extraction as bse
    from docking.preparation import excel_export as xl
    from docking.preparation import preparation_config as pc

    assert core_pipeline.extract_residue_level_coordinates is bse.extract_residue_level_coordinates
    assert core_pipeline.parse_plip_text_report is bse.parse_plip_text_report
    assert core_pipeline.centre_method_from_provenance is bse.centre_method_from_provenance
    assert core_pipeline.load_protein_structure is bse.load_protein_structure
    assert core_pipeline.excel_cell_value is xl.excel_cell_value
    assert core_pipeline.EXCEL_AVAILABLE is xl.EXCEL_AVAILABLE
    assert cli_pipeline.extract_residue_level_coordinates is bse.extract_residue_level_coordinates
    assert cli_pipeline.excel_cell_value is xl.excel_cell_value
    assert cli_pipeline.centre_method_from_provenance is bse.centre_method_from_provenance
    assert cli_pipeline.add_to_excel_workbook is xl.add_to_excel_workbook
    assert interactive_pipeline.add_to_excel_workbook is xl.add_to_excel_workbook
    assert autodock_preparation.PreparationConfig is pc.PreparationConfig
    assert autodock_preparation.profile_requires_autodocktools is pc.profile_requires_autodocktools


def test_moved_definitions_are_removed_from_the_legacy_modules():
    core_source = (REPO / "core_pipeline.py").read_text(encoding="utf-8")
    for name in ("def extract_residue_level_coordinates", "def parse_plip_text_report",
                 "def centre_method_from_provenance", "def excel_cell_value"):
        assert name not in core_source, name
    cli_source = (REPO / "cli_pipeline.py").read_text(encoding="utf-8")
    interactive_source = (REPO / "interactive_pipeline.py").read_text(encoding="utf-8")
    assert "def add_to_excel_workbook" not in cli_source
    assert "def add_to_excel_workbook" not in interactive_source
    ad_source = (REPO / "autodock_preparation.py").read_text(encoding="utf-8")
    assert "class PreparationConfig" not in ad_source
    assert "def profile_requires_autodocktools" not in ad_source


def test_class_methods_keep_their_signatures_and_seams(monkeypatch):
    from core_pipeline import MolecularDockingPipeline
    import docking.preparation.binding_site_extraction as bse

    pipeline = object.__new__(MolecularDockingPipeline)
    pipeline.plip_available = False

    # extract_active_site_coords: the instance's distance detector and the resolver are both honoured.
    pipeline.distance_based_interaction_detection = lambda *a, **k: [np.array([0.0, 0.0, 0.0])]
    monkeypatch.setattr(
        bse, "resolve_binding_site_center",
        lambda **kwargs: {"status": "completed", "binding_site_center": [1.0, 2.0, 3.0]},
    )
    coords, count = pipeline.extract_active_site_coords("x.pdb", "LIG", "A", 1, method="distance")
    assert count == 1
    np.testing.assert_allclose(coords, [1.0, 2.0, 3.0])

    # distance_based_interaction_detection: the instance's structure loader is honoured.
    pipeline2 = object.__new__(MolecularDockingPipeline)
    pipeline2._load_structure = lambda *a, **k: []
    with pytest.raises(ValueError, match="Ligand LIG not found"):
        pipeline2.distance_based_interaction_detection("x.pdb", "LIG", "A", 1)


# Modules whose public names this change moved, shimmed or re-exported.
TOUCHED_MODULES = {
    "core_pipeline", "cli_pipeline", "interactive_pipeline", "batch_pdb_preparation",
    "autodock_preparation", "test_pipeline", "main", "docking.cli",
    "docking.preparation.binding_site_extraction", "docking.preparation.excel_export",
    "docking.preparation.preparation_config", "workflow.legacy_entrypoints",
}


def test_nothing_imports_a_removed_path():
    """Every `from X import names` of a touched module must resolve to a file that defines the name.

    Scope: targets in TOUCHED_MODULES. Other packages (post_docking_analysis, webui) are owned by other
    specs; their pre-existing import issues are not this change's removals.
    """
    skip_dirs = {".git", "specs", "build", "dist", "__pycache__", ".pytest_cache", ".ruff_cache", "node_modules"}
    problems = []
    for path in REPO.rglob("*.py"):
        if any(part in skip_dirs or part.endswith(".egg-info") for part in path.relative_to(REPO).parts):
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                if node.module not in TOUCHED_MODULES:
                    continue
                target = _module_file(node.module)
                if target is None:
                    problems.append(f"{path.relative_to(REPO)}: {node.module} no longer exists")
                    continue
                defined, has_star = _top_level_names(target)
                if has_star:
                    continue  # names come from a star import; checked at runtime by the import tests
                for alias in node.names:
                    if alias.name == "*":
                        continue
                    if alias.name in defined:
                        continue
                    if _module_file(f"{node.module}.{alias.name}") is not None:
                        continue
                    problems.append(f"{path.relative_to(REPO)}: {node.module} has no {alias.name}")
    assert problems == []


def test_root_test_pipeline_is_gone_and_the_smoke_imports_the_moved_copy():
    assert not (REPO / "test_pipeline.py").exists()
    assert (REPO / "test" / "test_pipeline.py").exists()
    smoke = (REPO / "test" / "test_dockforge_smoke.py").read_text(encoding="utf-8")
    assert "from test_pipeline import run_smoke_test" in smoke


def test_prep_autodock_enhanced_is_documented_as_internal():
    head = (REPO / "prep_autodock_enhanced.sh").read_text(encoding="utf-8")[:1200]
    assert "NOT A USER ENTRY POINT" in head


# ---------------------------------------------------------------- R5: docs and help match the parser
def test_usage_shows_the_single_sequence_in_order():
    text = (REPO / "USAGE.md").read_text(encoding="utf-8")
    steps = [
        "workflow init", "workflow protonation-policy", "pdb prepare-protein", "pdb prepare-ligand",
        "prep pairlist", "prep project", "dock run", "analyze comparative", "analyze md-inputs",
    ]
    positions = []
    for step in steps:
        index = text.find(f"python main.py {step}")
        assert index >= 0, step
        positions.append(index)
    assert positions == sorted(positions)


def test_usage_and_readme_document_the_many_entity_path():
    usage = (REPO / "USAGE.md").read_text(encoding="utf-8")
    readme = (REPO / "README.md").read_text(encoding="utf-8")
    assert "## Many proteins and many ligands" in usage
    assert "--receptor <1IEP>" in usage and "--ligand <STI>" in usage
    assert "workflow protonation-policy --receptor NAME" in readme


@pytest.mark.parametrize("doc", ["USAGE.md", "README.md"])
def test_every_documented_command_and_flag_exists(doc):
    commands = _documented_commands(REPO / doc)
    assert commands, f"no documented commands found in {doc}"
    problems = []
    for tokens in commands:
        problems.extend(_check_command(tokens))
    assert problems == []


def test_main_help_examples_match_the_parser():
    epilog = _parser().epilog or ""
    commands = []
    for line in epilog.splitlines():
        stripped = line.strip()
        if stripped.startswith("python main.py"):
            commands.append([t for t in shlex.split(stripped)][2:])
    assert len(commands) >= 9
    problems = []
    for tokens in commands:
        problems.extend(_check_command(tokens))
    assert problems == []


def test_no_known_wrong_flags_remain_in_the_docs():
    usage = (REPO / "USAGE.md").read_text(encoding="utf-8")
    assert "--energy_range" not in usage
    readme = (REPO / "README.md").read_text(encoding="utf-8")
    assert "python test_pipeline.py" not in readme
    assert "\npython main.py batch -c" not in readme
    assert "\npython main.py cli " not in readme


# ---------------------------------------------------------------- packaging
@pytest.mark.integration
def test_wheel_builds_and_check_wheel_passes():
    dist = REPO / "dist"
    existed = {p: p.exists() for p in (REPO / "dist", REPO / "build")}
    egg_before = set(REPO.glob("*.egg-info"))
    try:
        if dist.exists():
            shutil.rmtree(dist)
        build = subprocess.run(
            [sys.executable, "-m", "build", "--wheel", "--outdir", str(dist)],
            cwd=REPO, capture_output=True, text=True, timeout=900,
        )
        assert build.returncode == 0, build.stderr[-2000:]
        check = subprocess.run(
            [sys.executable, str(REPO / "scripts" / "check_wheel.py")],
            cwd=REPO, capture_output=True, text=True, timeout=900,
        )
        assert check.returncode == 0, (check.stdout + check.stderr)[-2000:]
    finally:
        if not existed[REPO / "dist"] and dist.exists():
            shutil.rmtree(dist)
        if not existed[REPO / "build"] and (REPO / "build").exists():
            shutil.rmtree(REPO / "build")
        for egg in set(REPO.glob("*.egg-info")) - egg_before:
            shutil.rmtree(egg, ignore_errors=True)
