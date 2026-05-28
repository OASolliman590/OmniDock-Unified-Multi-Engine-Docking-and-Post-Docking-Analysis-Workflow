from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import statistics
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import yaml

from docking.models import (
    MAX_PREPARATION_PH,
    MIN_PREPARATION_PH,
    PairlistRow,
    validate_ligand_preparation_profile,
    validate_preparation_ph,
)
from docking.preparation.ligand_quality import audit_project_ligands
from docking.preparation.receptor_quality import audit_project_receptors
from docking.preparation.pairlist_builder import list_prepared_asset_names
from docking.project_layout import (
    deployment_root,
    detect_layout_profile,
    ensure_project_layout,
    load_manifest,
    load_pairlist,
    pair_intent_path,
    pairlist_path,
    post_docking_root,
    shared_receptors_dir,
)
from docking.hpc_profiles import load_hpc_profile, resolve_remote_project_dir, resolve_ssh_target
from docking.remote_ops import infer_round_mode_from_deployment_root
from .ids import format_run_timestamp
from .execution import (
    BackgroundTaskManager,
    run_analysis_comparative,
    run_analysis_favorite,
    run_analysis_target,
    run_autodock_prepare,
    run_docking,
    run_docking_deploy,
    run_docking_submit,
    run_docking_sync,
    run_pdb_collect,
    run_prepare_pairlist,
    run_prepare_project,
    run_workflow_clone_checkpoint,
    run_workflow_init,
)
from .state import ensure_state, list_background_tasks, list_steps, summarize_state


ENGINE_CHOICES = [
    ("all", "All Panels"),
    ("gnina", "GNINA"),
    ("vina", "Vina"),
    ("smina", "Smina"),
    ("autodock4", "AutoDock4"),
]

ANALYSIS_LABELS = {
    "analyze.favorite_engine": "Favorite-engine continuation",
    "analyze.comparative": "Comparative analysis",
    "analyze.stage.hierarchical": "Binding affinity summary",
    "analyze.stage.polypharmacology": "Polypharmacology",
    "analyze.stage.rmsd": "RMSD analysis only",
    "analyze.stage.reports": "Reports only",
    "analyze.stage.visualizations": "Visualizations only",
    "analyze.stage.structure_quality": "Structure quality only",
    "analyze.interactions.pandamap": "PandaMap target (routes to clean interactions pipeline)",
    "analyze.interactions.prolif": "ProLIF target (routes to clean interactions pipeline)",
    "analyze.interactions.ligplot": "LigPlot target (routes to clean interactions pipeline)",
    "analyze.interactions.poseview": "PoseView only",
    "analyze.interactions.clean": "Clean interactions pipeline (chem-fix + PLIP + ProLIF)",
    "analyze.visuals.py3dmol": "py3Dmol only",
    "analyze.visuals.pymol": "PyMOL only",
}

TIMELINE_STEPS = [
    ("pdb.collect", "Fetch proteins from PDB"),
    ("pdb.prepare_both", "Prepare proteins/ligands"),
    ("prep.pairlist", "Build pairlist"),
    ("prep.project", "Prepare docking folders"),
    ("workflow.clone_checkpoint", "Checkpoint & Revise"),
    ("dock.run", "Run docking"),
    ("dock.deploy", "HPC deploy"),
    ("dock.sync", "HPC sync"),
    ("dock.submit", "HPC submit"),
    ("analyze.comparative", "Post-docking analysis"),
]

_TASK_MANAGER = BackgroundTaskManager(max_workers=2)

_BACK_TOKENS = {"back", "menu"}
_EXIT_TOKENS = {"exit", "quit"}
_INTERACTION_TARGETS_ROUTED_TO_CLEAN = {
    "analyze.interactions.pandamap",
    "analyze.interactions.prolif",
    "analyze.interactions.ligplot",
}


class _PromptNavigation(Exception):
    def __init__(self, action: str):
        super().__init__(action)
        self.action = action


def _questionary():
    try:
        import questionary
    except ImportError as exc:
        raise RuntimeError(
            "Interactive workflow requires questionary. Install it with `pip install questionary` or update the project environment."
        ) from exc
    return questionary


def _load_optional_ssh_systems(project_root: Path) -> Tuple[List[Dict[str, str]], str]:
    """
    Load optional SSH-system presets from .workflow/ssh_systems.{yaml,yml,json}.
    Returns (systems, source_path).
    """
    root = Path(project_root).expanduser().resolve()
    candidates = [
        root / ".workflow" / "ssh_systems.yaml",
        root / ".workflow" / "ssh_systems.yml",
        root / ".workflow" / "ssh_systems.json",
    ]
    source_path = ""
    payload: object = None
    for candidate in candidates:
        if not candidate.exists():
            continue
        source_path = str(candidate)
        with open(candidate, "r", encoding="utf-8") as handle:
            if candidate.suffix.lower() in {".yaml", ".yml"}:
                payload = yaml.safe_load(handle)
            else:
                payload = json.load(handle)
        break

    if payload is None:
        return [], source_path
    if not isinstance(payload, dict):
        return [], source_path

    systems_payload = payload.get("systems")
    normalized: List[Dict[str, str]] = []

    def _append_system(name: str, cfg: object) -> None:
        if not str(name or "").strip():
            return
        mapping = dict(cfg) if isinstance(cfg, dict) else {}
        normalized.append(
            {
                "name": str(name).strip(),
                "label": str(mapping.get("label", "") or "").strip(),
                "kind": str(mapping.get("kind", "") or "").strip().lower(),
                "hpc_profile": str(mapping.get("hpc_profile", "") or "").strip(),
                "hpc_profile_file": str(mapping.get("hpc_profile_file", "") or "").strip(),
                "ssh_target": str(mapping.get("ssh_target", "") or "").strip(),
                "remote_project_dir": str(mapping.get("remote_project_dir", "") or "").strip(),
            }
        )

    if isinstance(systems_payload, list):
        for item in systems_payload:
            if not isinstance(item, dict):
                continue
            _append_system(str(item.get("name", "") or "").strip(), item)
    elif isinstance(systems_payload, dict):
        for name, cfg in systems_payload.items():
            _append_system(str(name), cfg)

    return normalized, source_path


def _probe_ssh_connectivity(ssh_target: str, timeout_seconds: int = 6) -> Tuple[bool, str]:
    target = str(ssh_target or "").strip()
    if not target:
        return False, "no SSH target provided"
    try:
        probe = subprocess.run(
            [
                "ssh",
                "-o",
                "BatchMode=yes",
                "-o",
                f"ConnectTimeout={int(timeout_seconds)}",
                target,
                "exit",
            ],
            capture_output=True,
            text=True,
            timeout=max(2, int(timeout_seconds) + 2),
            check=False,
        )
    except FileNotFoundError:
        return False, "ssh client is not installed on this machine"
    except subprocess.TimeoutExpired:
        return False, f"probe timed out after {timeout_seconds}s"

    stderr_text = str(probe.stderr or "").strip()
    if probe.returncode == 0:
        return True, "SSH login probe succeeded"
    if probe.returncode == 255 and stderr_text:
        return False, stderr_text
    return False, f"ssh probe failed with return code {probe.returncode}"


def _as_float(value: str) -> float:
    token = str(value or "").strip()
    if token in {"", "undefined", "UNDEFINED", "nan", "NaN"}:
        return 0.0
    try:
        return float(token)
    except Exception:
        return 0.0


def _parse_condor_slot_snapshot(raw_output: str) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    for line in str(raw_output or "").splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 9:
            continue
        rows.append(
            {
                "state": str(parts[1]).strip().lower(),
                "activity": str(parts[2]).strip().lower(),
                "cpus": _as_float(parts[3]),
                "memory_mb": _as_float(parts[4]),
                "disk_kb": _as_float(parts[5]),
                "gpus": _as_float(parts[6]),
                "gpu_capability": _as_float(parts[7]),
                "total_gpus": _as_float(parts[8]),
            }
        )
    return rows


def _recommended_gpu_capability(capabilities: List[float]) -> float:
    valid = sorted(cap for cap in capabilities if cap > 0)
    if not valid:
        return 0.0
    if any(cap >= 8.9 for cap in valid):
        return 8.9
    if any(cap >= 7.5 for cap in valid):
        return 7.5
    return round(valid[0], 1)


def _probe_nmrbox_condor_fit(ssh_target: str, timeout_seconds: int = 10) -> Tuple[bool, Dict[str, object], str]:
    target = str(ssh_target or "").strip()
    if not target:
        return False, {}, "no SSH target provided"
    remote_cmd = (
        "condor_status -af Name State Activity TotalSlotCpus TotalSlotMemory Disk GPUs GPUs_Capability TotalGpus "
        "2>/dev/null || "
        "condor_status -af Name State Activity Cpus Memory Disk Gpus GPUs_Capability TotalGpus 2>/dev/null"
    )
    try:
        probe = subprocess.run(
            [
                "ssh",
                "-o",
                "BatchMode=yes",
                "-o",
                f"ConnectTimeout={int(timeout_seconds)}",
                target,
                remote_cmd,
            ],
            capture_output=True,
            text=True,
            timeout=max(4, int(timeout_seconds) + 4),
            check=False,
        )
    except FileNotFoundError:
        return False, {}, "ssh client is not installed on this machine"
    except subprocess.TimeoutExpired:
        return False, {}, f"condor_status probe timed out after {timeout_seconds}s"

    if probe.returncode != 0:
        details = str(probe.stderr or "").strip() or f"ssh return code {probe.returncode}"
        return False, {}, details

    slots = _parse_condor_slot_snapshot(probe.stdout)
    if not slots:
        return False, {}, "condor_status returned no parseable slots"

    idle_slots = [
        row
        for row in slots
        if row.get("state") == "unclaimed" and row.get("activity") == "idle"
    ]
    gpu_idle_slots = [
        row
        for row in idle_slots
        if row.get("gpus", 0.0) > 0 or row.get("total_gpus", 0.0) > 0
    ]

    cpu_values = [row.get("cpus", 0.0) for row in idle_slots if row.get("cpus", 0.0) > 0]
    memory_values = [row.get("memory_mb", 0.0) for row in idle_slots if row.get("memory_mb", 0.0) > 0]
    disk_values = [row.get("disk_kb", 0.0) for row in idle_slots if row.get("disk_kb", 0.0) > 0]
    capability_values = [row.get("gpu_capability", 0.0) for row in gpu_idle_slots if row.get("gpu_capability", 0.0) > 0]

    cpu_rec = int(max(4, min(16, round(statistics.median(cpu_values) if cpu_values else 8))))
    mem_mb_rec = statistics.median(memory_values) if memory_values else 32768.0
    mem_gb_rec = int(max(16, min(64, round(mem_mb_rec / 1024.0))))
    disk_kb_rec = statistics.median(disk_values) if disk_values else float(20 * 1024 * 1024)
    disk_gb_rec = int(max(20, min(80, round(disk_kb_rec / (1024.0 * 1024.0)))))

    gpu_capability_rec = _recommended_gpu_capability(capability_values)
    requirement = ""
    if gpu_capability_rec > 0:
        requirement = f"(GPUs >= 1) && (GPUs_Capability >= {gpu_capability_rec:.1f})"

    payload: Dict[str, object] = {
        "slots_total": len(slots),
        "idle_slots_total": len(idle_slots),
        "idle_gpu_slots_total": len(gpu_idle_slots),
        "recommended_cpus": cpu_rec,
        "recommended_memory": f"{mem_gb_rec}GB",
        "recommended_disk": f"{disk_gb_rec}GB",
        "recommended_gpus": 1 if gpu_idle_slots else 0,
        "recommended_requirements": requirement,
        "recommended_gpu_capability": gpu_capability_rec,
    }
    return True, payload, "condor pool probe succeeded"


def _collapse_interaction_alias_targets(selected_targets: List[str]) -> tuple[List[str], List[str]]:
    """
    Enforce clean-interaction contract routing for legacy standalone interaction targets.
    Returns (normalized_targets, collapsed_alias_targets).
    """
    normalized: List[str] = []
    collapsed: List[str] = []
    seen: set[str] = set()
    for target in selected_targets:
        mapped_target = target
        if target in _INTERACTION_TARGETS_ROUTED_TO_CLEAN:
            collapsed.append(target)
            mapped_target = "analyze.interactions.clean"
        if mapped_target in seen:
            continue
        normalized.append(mapped_target)
        seen.add(mapped_target)
    return normalized, collapsed


def _resolve_run_tracking_manifest_path(result) -> Optional[Path]:
    outputs = getattr(result, "outputs", {})
    if not isinstance(outputs, dict):
        outputs = {}

    explicit_keys = ("run_manifest_file", "run_tracking_manifest_file")
    for key in explicit_keys:
        token = str(outputs.get(key, "") or "").strip()
        if not token:
            continue
        candidate = Path(token).expanduser()
        try:
            resolved = candidate.resolve()
        except Exception:
            continue
        if resolved.exists() and resolved.is_file():
            return resolved

    output_dir_token = str(outputs.get("output_dir", "") or "").strip()
    if not output_dir_token:
        return None

    output_dir = Path(output_dir_token).expanduser()
    try:
        output_dir = output_dir.resolve()
    except Exception:
        return None

    candidates = [
        output_dir / "run_tracking" / "run_manifest.json",
        output_dir / "favorite_engine" / "run_tracking" / "run_manifest.json",
    ]
    for candidate in candidates:
        if candidate.exists() and candidate.is_file():
            return candidate.resolve()
    return None


def _format_stage_contract_diff(stage_contract: Dict[str, object]) -> str:
    if not isinstance(stage_contract, dict) or not stage_contract:
        return "Stage contract: unavailable"

    tracked = ("run_rmsd", "run_visualizations", "run_prolif", "run_ligplot")
    parts: List[str] = []
    for stage in tracked:
        requested_key = f"requested_{stage}"
        effective_key = f"effective_{stage}"
        if requested_key not in stage_contract and effective_key not in stage_contract:
            continue
        requested_value = bool(stage_contract.get(requested_key, False))
        effective_value = bool(stage_contract.get(effective_key, False))
        parts.append(f"{stage} requested={requested_value} enforced={effective_value}")

    if not parts:
        return "Stage contract: unavailable"
    return "Stage contract: " + " | ".join(parts)


def _print_run_tracking_summary(result) -> None:
    manifest_path = _resolve_run_tracking_manifest_path(result)
    if manifest_path is None:
        print("Run tracking manifest: unavailable for this target")
        print("Stage contract: unavailable")
        return

    print(f"Run tracking manifest: {manifest_path}")
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"Stage contract: unavailable (failed to read manifest: {exc})")
        return

    stage_contract = payload.get("stage_contract")
    if not isinstance(stage_contract, dict):
        print("Stage contract: unavailable")
        return
    print(_format_stage_contract_diff(stage_contract))


def _header(title: str) -> None:
    print(f"\n{title}")
    print("-" * len(title))
    print("Type 'back' or select Back at any prompt to return. Type 'exit' to quit.\n")


def _normalize_prompt_token(value: object) -> str:
    return _strip_wrapping_quotes(str(value or "").strip()).lower()


def _maybe_navigation(value: object) -> None:
    token = _normalize_prompt_token(value)
    if token in _EXIT_TOKENS:
        raise _PromptNavigation("exit")
    if token in _BACK_TOKENS:
        raise _PromptNavigation("back")


def _select(
    questionary,
    message: str,
    choices: List[tuple[str, str]],
    default: Optional[str] = None,
    include_back: bool = True,
) -> str:
    select_choices = [
        questionary.Choice(title=label, value=value, checked=value == default)
        for value, label in choices
    ]
    if include_back:
        select_choices.append(questionary.Choice(title="Back", value="back"))
    result = questionary.select(
        message,
        choices=select_choices,
        default=default,
    ).ask()
    if result is None:
        raise _PromptNavigation("back")
    _maybe_navigation(result)
    return str(result)


def _strip_wrapping_quotes(value: str) -> str:
    text = value.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {"'", '"'}:
        return text[1:-1].strip()
    return text


def _split_csv(value: str) -> List[str]:
    return [item.strip().upper() for item in value.split(",") if item.strip()]


def _text(questionary, message: str, default: Optional[str] = None, required: bool = False) -> str:
    while True:
        value = questionary.text(message, default=default or "").ask()
        if value is None:
            raise _PromptNavigation("back")
        value = _strip_wrapping_quotes((value or "").strip())
        _maybe_navigation(value)
        if value:
            return value
        if default is not None:
            return default
        if not required:
            return ""
        print("This value is required. Type 'back' to return or 'exit' to quit.")


def _optional_text(questionary, message: str, default: Optional[str] = None) -> str:
    value = questionary.text(message, default=default or "").ask()
    if value is None:
        raise _PromptNavigation("back")
    _maybe_navigation(value)
    return _strip_wrapping_quotes((value or "").strip())


def _confirm(questionary, message: str, default: bool = True) -> bool:
    result = questionary.confirm(message, default=default).ask()
    if result is None:
        raise _PromptNavigation("back")
    return bool(result)


def _float_input(questionary, message: str, default: float) -> float:
    while True:
        raw = _text(questionary, message, default=str(default), required=True)
        try:
            return float(raw)
        except ValueError:
            print("Please enter a valid numeric value (example: 0.0, -7.5).")


def _int_input(questionary, message: str, default: int, minimum: Optional[int] = None) -> int:
    while True:
        raw = _text(questionary, message, default=str(default), required=True)
        try:
            value = int(raw)
        except ValueError:
            print("Please enter a whole number.")
            continue
        if minimum is not None and value < minimum:
            print(f"Please enter an integer >= {minimum}.")
            continue
        return value


def _parse_rmsd_workers(raw: str) -> int:
    token = str(raw or "").strip().lower()
    if not token:
        return 0
    if token in {"auto", "default"}:
        return 0
    if token in {"single", "single-thread", "single_thread"}:
        return 1
    if token in {"parallel", "multi", "multi-thread", "multi_thread", ">1", ">=2"}:
        return 2
    if token.startswith(">"):
        try:
            parsed = int(token[1:])
            return 2 if parsed <= 1 else parsed
        except ValueError as exc:
            raise ValueError("Invalid RMSD worker count") from exc
    value = int(token)
    if value < 0:
        raise ValueError("RMSD worker count must be >= 0")
    return value


def _rmsd_workers_input(questionary, default: int = 0) -> int:
    while True:
        raw = _text(
            questionary,
            "RMSD workers (0=auto, 1=single-thread, >1 parallel workers; also accepts auto/single/parallel)",
            default=str(default),
            required=True,
        )
        try:
            return _parse_rmsd_workers(raw)
        except ValueError:
            print("Please enter 0, 1, an integer >1, or one of: auto, single, parallel.")


def _checkbox(
    questionary,
    message: str,
    values: List[object],
    preselected: Optional[List[str]] = None,
) -> List[str]:
    preselected_set = {str(value) for value in (preselected or [])}
    choices = []
    for item in values:
        if isinstance(item, tuple) and len(item) >= 2:
            value, label = str(item[0]), str(item[1])
        else:
            value, label = str(item), str(item)
        choices.append(questionary.Choice(title=label, value=value, checked=value in preselected_set))
    selected = questionary.checkbox(
        message,
        choices=choices,
    ).ask()
    if selected is None:
        raise _PromptNavigation("back")
    return list(selected or [])


def _analysis_engine_scope_picker(
    questionary,
    manifest: Dict[str, object],
    detection: Dict[str, object],
    *,
    preselected: Optional[List[str]] = None,
    preset_name: str = "",
) -> Tuple[List[str], str]:
    valid_engines = [
        str(engine).strip().lower()
        for engine in (detection.get("valid_engines") or [])
        if str(engine).strip()
    ]
    if len(valid_engines) <= 1:
        return valid_engines, preset_name

    presets = manifest.get("engine_presets") if isinstance(manifest.get("engine_presets"), dict) else {}
    selected = list(preselected or valid_engines)
    selected = [engine for engine in selected if engine in valid_engines] or list(valid_engines)
    active_preset = str(preset_name or "")

    while True:
        if presets:
            mode = _select(
                questionary,
                "Engine scope source",
                [
                    ("manual", "Choose engines manually"),
                    ("preset", "Load a saved preset"),
                ],
                default="preset" if active_preset else "manual",
            )
            if mode == "preset":
                preset_options = [(name, f"{name}: {', '.join(map(str, values))}") for name, values in presets.items()]
                active_preset = _select(
                    questionary,
                    "Engine preset",
                    preset_options,
                    default=active_preset or next(iter(presets.keys())),
                )
                raw_preset = presets.get(active_preset) or []
                selected = [
                    str(engine).strip().lower()
                    for engine in raw_preset
                    if str(engine).strip().lower() in valid_engines
                ]
                if not selected:
                    print("Preset does not contain any currently valid engines. Please choose again.")
                    continue
                return selected, active_preset

        engine_rows = []
        engines_payload = detection.get("engines") if isinstance(detection.get("engines"), dict) else {}
        for engine in valid_engines:
            payload = engines_payload.get(engine) or {}
            coverage = payload.get("coverage_pct")
            coverage_label = f"{float(coverage):.1f}% coverage" if coverage not in (None, "") else "coverage unknown"
            engine_rows.append((engine, f"{engine.upper()} ({coverage_label})"))

        selected = _checkbox(
            questionary,
            "Select analysis engines (space to toggle, enter to confirm)",
            engine_rows,
            preselected=selected,
        )
        if not selected:
            print("Please keep at least one engine selected.")
            continue
        return selected, ""


def _resolve_engine_selection(selection: str) -> List[str]:
    if selection == "all":
        return ["gnina", "vina", "smina", "autodock4"]
    return [selection]


def _project_layout(project_root: Path) -> Dict[str, Path]:
    return ensure_project_layout(project_root, "docking_legacy")


def _shared_naming_dir(project_root: Path) -> Path:
    root = Path(project_root).expanduser().resolve()
    analysis_root = post_docking_root(root)
    naming_dir = analysis_root / "_naming"
    naming_dir.mkdir(parents=True, exist_ok=True)
    return naming_dir


def _resolve_existing_path(project_root: Path, candidate: str) -> str:
    token = _strip_wrapping_quotes(str(candidate or "").strip())
    if not token:
        return ""
    path = Path(token).expanduser()
    if not path.is_absolute():
        path = (project_root / path).resolve()
    else:
        path = path.resolve()
    return str(path) if path.exists() else ""


def _default_rerun_manifest(project_root: Path, manifest: Dict[str, object]) -> str:
    candidate = _resolve_existing_path(project_root, str(manifest.get("latest_rerun_manifest_file", "") or ""))
    return candidate


def _load_local_deployment_manifest(project_root: Path, round_id: str, mode: str) -> Dict[str, object]:
    manifest_path = deployment_root(project_root) / round_id / mode / "deployment_manifest.json"
    if not manifest_path.exists():
        return {}
    with open(manifest_path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return payload if isinstance(payload, dict) else {}


def _suggest_hpc_round_label(project_root: Path, mode: str, inferred_round: str = "", inferred_mode: str = "") -> str:
    if inferred_round and inferred_mode == mode:
        return inferred_round
    prefix = f"hpc_{mode}_"
    next_index = 1
    root = deployment_root(project_root)
    pattern = re.compile(rf"^{re.escape(prefix)}(\d+)$")
    if root.exists():
        for child in root.iterdir():
            match = pattern.match(child.name)
            if match:
                next_index = max(next_index, int(match.group(1)) + 1)
    return f"{prefix}{next_index:03d}"


def _create_project(questionary) -> Path:
    _header("Create Docking Project")
    working_dir = Path(_text(questionary, "Working directory", default=str(Path.cwd()), required=True)).expanduser().resolve()
    project_name = _text(questionary, "Project name", required=True)
    engine_selection = _select(questionary, "Prepare docking folders for which panels?", ENGINE_CHOICES, default="all")
    project_root = working_dir / project_name
    result = run_workflow_init(
        str(project_root),
        engines=_resolve_engine_selection(engine_selection),
        project_name=project_name,
        layout_profile="docking_legacy",
    )
    print(f"Project initialized at {result.root}")
    return Path(result.root)


def _resume_project(questionary, initial_root: Optional[Path] = None) -> Path:
    suggested_root = Path(initial_root).expanduser().resolve() if initial_root else Path.cwd().resolve()
    choice = _select(
        questionary,
        "Project root selection",
        [
            ("suggested", "Use detected project root (recommended)"),
            ("paste", "Paste another project path"),
        ],
        default="suggested",
    )
    if choice == "suggested":
        root = suggested_root
    else:
        root = Path(
            _text(
                questionary,
                "Project root",
                default="",
                required=True,
            )
        ).expanduser().resolve()
    ensure_state(root)
    if not (root / "4-Docking" / "project_manifest.json").exists() and not (root / "project_manifest.json").exists():
        run_workflow_init(str(root), layout_profile="docking_legacy")
    return root


def _suggest_startup_project_root(initial_root: Optional[Path] = None) -> Path:
    """
    Suggest a local project root when interactive workflow starts.

    Priority:
    1) explicit initial_root argument
    2) repository root inferred from this file location
    3) current working directory
    """
    if initial_root is not None:
        return Path(initial_root).expanduser().resolve()

    repo_root = Path(__file__).resolve().parents[1]
    if repo_root.exists():
        return repo_root
    return Path.cwd().resolve()


def _suggest_prepare_ligand4_path() -> str:
    for command_name in ("prepare_ligand4.py", "prepare_ligand4"):
        resolved = shutil.which(command_name)
        if resolved:
            return str(Path(resolved).expanduser().resolve())

    for env_key in ("AUTODOCKTOOLS_PREPARE_LIGAND4", "ADT_PREPARE_LIGAND4"):
        env_value = str(os.environ.get(env_key, "") or "").strip()
        if not env_value:
            continue
        candidate = Path(env_value).expanduser()
        if candidate.exists() and candidate.is_file():
            return str(candidate.resolve())

    mgltools_root = str(os.environ.get("MGLTOOLS_PATH", "") or "").strip()
    if mgltools_root:
        candidate = Path(mgltools_root).expanduser() / "MGLToolsPckgs" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"
        if candidate.exists() and candidate.is_file():
            return str(candidate.resolve())

    repo_root = Path(__file__).resolve().parents[1]
    for candidate in [
        repo_root / ".workflow" / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py",
        repo_root / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py",
    ]:
        if candidate.exists() and candidate.is_file():
            return str(candidate.resolve())

    return ""


def _suggest_prepare_receptor4_path(ligand_script: str = "") -> str:
    ligand_hint = str(ligand_script or "").strip()
    if ligand_hint:
        ligand_path = Path(ligand_hint).expanduser()
        derived = ligand_path.with_name("prepare_receptor4.py")
        if derived.exists() and derived.is_file():
            return str(derived.resolve())

    for command_name in ("prepare_receptor4.py", "prepare_receptor4"):
        resolved = shutil.which(command_name)
        if resolved:
            return str(Path(resolved).expanduser().resolve())

    for env_key in ("AUTODOCKTOOLS_PREPARE_RECEPTOR4", "ADT_PREPARE_RECEPTOR4"):
        env_value = str(os.environ.get(env_key, "") or "").strip()
        if not env_value:
            continue
        candidate = Path(env_value).expanduser()
        if candidate.exists() and candidate.is_file():
            return str(candidate.resolve())

    mgltools_root = str(os.environ.get("MGLTOOLS_PATH", "") or "").strip()
    if mgltools_root:
        candidate = Path(mgltools_root).expanduser() / "MGLToolsPckgs" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py"
        if candidate.exists() and candidate.is_file():
            return str(candidate.resolve())

    repo_root = Path(__file__).resolve().parents[1]
    for candidate in [
        repo_root / ".workflow" / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py",
        repo_root / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py",
    ]:
        if candidate.exists() and candidate.is_file():
            return str(candidate.resolve())

    return ""


def _show_status(project_root: Path) -> None:
    _header("Workflow Status")
    _show_workflow_timeline(project_root)
    print("")
    for line in summarize_state(project_root):
        print(line)


def _colorize(text: str, color_code: str) -> str:
    return f"\033[{color_code}m{text}\033[0m"


def _step_symbol(status: str) -> str:
    normalized = str(status or "").strip()
    if normalized == "validated":
        return _colorize("✔", "32")
    if normalized == "completed":
        return _colorize("✔", "32")
    if normalized in {"in_progress", "running"}:
        return _colorize("▶", "33")
    if normalized in {"failed", "blocked"}:
        return _colorize("✖", "31")
    if normalized == "needs_review":
        return _colorize("!", "35")
    if normalized == "skipped":
        return _colorize("⏭", "90")
    return "○"


def _derive_timeline_status(steps: Dict[str, Dict[str, object]], target: str) -> str:
    if target in steps:
        return str(steps[target].get("status", "not_started")).strip() or "not_started"
    if target == "pdb.prepare_both":
        statuses = []
        for fallback_target in ("pdb.prepare_both", "pdb.prepare_protein", "pdb.prepare_ligand"):
            if fallback_target in steps:
                statuses.append(str(steps[fallback_target].get("status", "")))
        if statuses:
            if any(state in {"failed", "blocked"} for state in statuses):
                return "needs_review"
            if any(state in {"in_progress", "running"} for state in statuses):
                return "in_progress"
            if all(state in {"completed", "validated"} for state in statuses):
                return "validated"
            if any(state == "partial" for state in statuses):
                return "needs_review"
    if target == "analyze.comparative":
        analyze_targets = [name for name in steps if str(name).startswith("analyze.")]
        if analyze_targets:
            analyze_statuses = [str(steps[name].get("status", "")) for name in analyze_targets]
            if any(state in {"in_progress", "running"} for state in analyze_statuses):
                return "in_progress"
            if any(state in {"failed", "blocked"} for state in analyze_statuses):
                return "needs_review"
            if any(state == "partial" for state in analyze_statuses):
                return "needs_review"
            if all(state in {"completed", "validated"} for state in analyze_statuses):
                return "validated"
    return "not_started"


def _show_workflow_timeline(project_root: Path) -> None:
    steps = list_steps(project_root)
    print("\nWorkflow Timeline")
    print("-----------------")
    for index, (target, label) in enumerate(TIMELINE_STEPS, start=1):
        status = _derive_timeline_status(steps, target)
        symbol = _step_symbol(status)
        status_text = status
        if status in {"validated", "completed"}:
            status_text = _colorize(status, "32")
        elif status in {"in_progress", "running"}:
            status_text = _colorize(status, "33")
        elif status in {"failed", "blocked", "needs_review"}:
            status_text = _colorize(status, "31")
        print(f"{index:>2}. {symbol} {label} [{status_text}]")
    active_tasks = list_background_tasks(project_root, include_finished=False)
    if active_tasks:
        print("\nActive tasks:")
        for task in active_tasks[:5]:
            label = str(task.get("label", "Task"))
            progress = int(task.get("progress", 0) or 0)
            status = str(task.get("status", "unknown"))
            print(f" - {label}: {status} ({progress}%)")


def _show_background_tasks(project_root: Path) -> None:
    _header("Background Tasks")
    rows = list_background_tasks(project_root, include_finished=True)
    if not rows:
        print("No background tasks recorded yet.")
        return
    for row in rows:
        task_id = str(row.get("task_id", ""))
        label = str(row.get("label", ""))
        status = str(row.get("status", "unknown"))
        progress = int(row.get("progress", 0) or 0)
        note = str(row.get("note", ""))
        updated_at = str(row.get("updated_at", ""))
        print(f"- {task_id} | {label} | status={status} | progress={progress}% | updated={updated_at}")
        if note:
            print(f"  note: {note}")
        error = str(row.get("error", "")).strip()
        if error:
            preview = error.strip().splitlines()
            if preview:
                print(f"  error: {preview[-1]}")


def _create_checkpoint_workspace(questionary, project_root: Path) -> None:
    _header("Checkpoint & Revise")
    source_root = _text(
        questionary,
        "Source project root",
        default=str(Path(project_root).expanduser().resolve()),
        required=True,
    )
    source_path = Path(source_root).expanduser().resolve()
    target_default = source_path.parent / f"{source_path.name}_checkpoint"
    target_root = _text(
        questionary,
        "New checkpoint project root",
        default=str(target_default),
        required=True,
    )

    print("Creating checkpoint workspace clone...")
    try:
        result = run_workflow_clone_checkpoint(
            source_project_dir=source_root,
            target_project_dir=target_root,
        )
    except Exception as exc:
        print(f"❌ Could not create checkpoint workspace: {exc}")
        return

    print(f"✅ Checkpoint workspace created: {result.outputs.get('target_project_dir', '')}")
    print("   Original project remains untouched.")


def _collect_pdbs(questionary, project_root: Path) -> None:
    _header("Fetch Proteins From PDB")
    selection_mode = _select(
        questionary,
        "Ligand extraction mode",
        [
            ("heuristic", "Smart case-by-case detection"),
            ("auto", "Auto parsing (legacy)"),
            ("interactive", "One-by-one parsing (manual)"),
        ],
        default="heuristic",
    )
    preserve_metals = _confirm(questionary, "Preserve common metals during cleaning (ZN, MN, MG, ...)?", default=True)
    preserve_cofactors = _confirm(questionary, "Preserve common cofactors during cleaning (NAD/NDP/NAP/FAD...)?", default=True)
    preserve_full_receptor = _confirm(
        questionary,
        "Keep full protein unchanged (skip cleaning entirely)?",
        default=False,
    )
    preserve_residues = _split_csv(
        _text(
            questionary,
            "Extra residues to preserve (comma-separated, optional)",
            default="",
            required=False,
        )
    )
    preferred_ligands: List[str] = []
    if selection_mode == "heuristic":
        preferred_ligands = _split_csv(
            _text(
                questionary,
                "Preferred anchor ligands (comma-separated, optional)",
                default="",
                required=False,
            )
        )

    while True:
        pdbs = _text(questionary, "PDB ID(s) (single, comma-separated, or file path)", required=True)
        result = run_pdb_collect(
            str(project_root),
            pdbs,
            selection_mode=selection_mode,
            preserve_metals=preserve_metals,
            preserve_cofactors=preserve_cofactors,
            preserve_full_receptor=preserve_full_receptor,
            preserve_residues=preserve_residues,
            preferred_ligands=preferred_ligands,
        )
        print(f"Collected into {result.outputs.get('raw_proteins_dir', '')}")
        if result.notes:
            print("Collection notes:")
            for note in result.notes:
                print(f" - {note}")
        if not _confirm(questionary, "Fetch another PDB?", default=False):
            return


def _prepare_assets(questionary, project_root: Path) -> None:
    _header("Prepare Proteins And Ligands")
    layout = _project_layout(project_root)
    manifest = load_manifest(project_root)
    selected_engines = [str(engine).strip().lower() for engine in (manifest.get("engines", []) or []) if str(engine).strip()]
    action = _select(
        questionary,
        "What do you want to prepare?",
        [
            ("pdb.prepare_both", "Prepare proteins and ligands"),
            ("pdb.prepare_protein", "Prepare proteins only"),
            ("pdb.prepare_ligand", "Prepare ligands only"),
        ],
        default="pdb.prepare_both",
    )
    receptors_input = _text(questionary, "Raw proteins input", default=str(layout["raw_proteins"])) if action != "pdb.prepare_ligand" else ""
    ligands_input = _text(questionary, "Raw ligands input", default=str(layout["raw_ligands"])) if action != "pdb.prepare_protein" else ""
    receptors_output = _text(questionary, "Prepared proteins output", default=str(layout["prepared_proteins"])) if action != "pdb.prepare_ligand" else ""
    ligands_output = _text(questionary, "Prepared ligands output", default=str(layout["prepared_ligands"])) if action != "pdb.prepare_protein" else ""
    force_field = _text(questionary, "Force field", default="AMBER", required=True)
    while True:
        raw_ph = _text(questionary, "pH", default="7.4", required=True)
        try:
            parsed_ph = float(raw_ph)
        except ValueError:
            print("pH must be a numeric value.")
            continue
        ph_ok, normalized_ph, ph_error = validate_preparation_ph(parsed_ph)
        if ph_ok:
            ph = normalized_ph
            break
        print(f"❌ {ph_error} (allowed range: {MIN_PREPARATION_PH:.1f}-{MAX_PREPARATION_PH:.1f})")
    while True:
        ligand_backend = _select(
            questionary,
            "Ligand preparation profile",
            [
                ("engine_aware_full", "Engine-aware full (Open Babel -> Meeko, with AutoDockTools branch when AD4 selected)"),
                ("openbabel_only", "Open Babel only"),
                ("meeko_only", "Meeko only"),
                ("autodocktools_only", "AutoDockTools only"),
                ("openbabel_meeko", "Open Babel -> Meeko"),
                ("openbabel_meeko_autodock", "Open Babel -> Meeko -> AutoDockTools fallback"),
                ("openbabel_autodocktools", "Open Babel -> AutoDockTools"),
            ],
            default="engine_aware_full",
        )
        compatibility = validate_ligand_preparation_profile(ligand_backend, selected_engines)
        if compatibility.is_valid:
            for warning in compatibility.warnings:
                print(f"⚠️  {warning}")
            break
        print("❌ Incompatible preparation profile for selected engines:")
        for issue in compatibility.errors:
            print(f"   - {issue}")
        print("Please choose another ligand preparation profile.")
    autodocktools_prepare_ligand4 = ""
    autodocktools_prepare_receptor4 = ""
    autodocktools_python = ""
    if ligand_backend in {"autodocktools_only", "openbabel_autodocktools", "openbabel_meeko_autodock"}:
        autodocktools_default = _suggest_prepare_ligand4_path()
        autodocktools_prepare_ligand4 = _text(
            questionary,
            "Path to prepare_ligand4.py (leave empty to use PATH lookup)",
            default=autodocktools_default,
            required=False,
        )
        autodocktools_prepare_receptor4 = _text(
            questionary,
            "Path to prepare_receptor4.py (optional; used only for AutoDockTools-only receptor prep)",
            default=_suggest_prepare_receptor4_path(autodocktools_prepare_ligand4 or autodocktools_default),
            required=False,
        )
        autodocktools_python = _text(
            questionary,
            "Python executable for AutoDockTools (optional)",
            default="",
            required=False,
        )
    result = run_autodock_prepare(
        action,
        receptors_input or None,
        ligands_input or None,
        receptors_output or None,
        ligands_output or None,
        force_field=force_field,
        ph=ph,
        ligand_preparation_backend=ligand_backend,
        selected_engines=selected_engines,
        autodocktools_prepare_ligand4=autodocktools_prepare_ligand4 or None,
        autodocktools_prepare_receptor4=autodocktools_prepare_receptor4 or None,
        autodocktools_python=autodocktools_python or None,
    )
    print(f"Preparation status: {result.status}")
    ligand_count = result.outputs.get("ligand_prepared_count")
    receptor_count = result.outputs.get("receptor_prepared_count")
    if ligand_count is not None or receptor_count is not None:
        print(f"Prepared assets: ligands={ligand_count or 0}, receptors={receptor_count or 0}")
    for note in result.notes:
        print(f" - {note}")


def _build_pairlist(questionary, project_root: Path) -> None:
    _header("Build Pairlist")
    layout = _project_layout(project_root)
    excel_file = layout["raw_proteins"] / "multi_pdb_analysis.xlsx"
    if not excel_file.exists():
        nested_candidates = sorted(project_root.glob("*/2-Raw_Protien/multi_pdb_analysis.xlsx"))
        print(f"❌ Required workbook not found: {excel_file}")
        if nested_candidates:
            print("   The selected root looks like a parent directory, not the actual docking project root.")
            if len(nested_candidates) == 1:
                suggested_root = nested_candidates[0].parent.parent
                print(f"   Suggested project root: {suggested_root}")
            else:
                print("   Multiple candidate projects were found under this parent:")
                for candidate in nested_candidates:
                    print(f"    - {candidate.parent.parent}")
        else:
            print("   Run the 'Fetch proteins from PDB' stage first, or switch to the correct project root.")
        return
    workflow_mode = _select(
        questionary,
        "How do you want to handle pair-curation state?",
        [
            ("direct", "Build pairlist directly"),
            ("iterative_only", "Save/update a reusable pair-curation round"),
            ("iterative_freeze", "Save/update a round and freeze it into pairlist.csv"),
            ("freeze_existing", "Freeze an existing saved round into pairlist.csv"),
        ],
        default="direct",
    )
    try:
        latest_round = str(load_manifest(project_root).get("latest_pair_round") or "").strip()
    except FileNotFoundError:
        latest_round = ""
    mode: Optional[str] = None
    if workflow_mode != "freeze_existing":
        mode = _select(
        questionary,
        "How do you want to build the pairlist?",
        [
            ("cocrystal_only", "Cocrystal benchmark only"),
            ("cocrystal_plus_all", "Cocrystal benchmark + all proteins to all ligands"),
            ("cocrystal_plus_nonreference", "Own cocrystal per protein + all non-reference ligands"),
            ("curated_cartesian", "Curated protein and ligand selection"),
            ("curated_per_protein", "Curated per-protein ligand selection"),
        ],
        default="cocrystal_plus_all",
        )
    curated_receptors: Optional[List[str]] = None
    curated_ligands: Optional[List[str]] = None
    curated_mapping: Optional[Dict[str, List[str]]] = None

    if mode == "curated_cartesian":
        receptor_names = list_prepared_asset_names(layout["prepared_proteins"], "receptor")
        ligand_names = list_prepared_asset_names(layout["prepared_ligands"], "ligand")
        curated_receptors = _checkbox(questionary, "Select proteins", receptor_names)
        curated_ligands = _checkbox(questionary, "Select ligands", ligand_names)
    elif mode == "curated_per_protein":
        receptor_names = _checkbox(questionary, "Select proteins", list_prepared_asset_names(layout["prepared_proteins"], "receptor"))
        ligand_pool = list_prepared_asset_names(layout["prepared_ligands"], "ligand")
        curated_mapping = {}
        for receptor_name in receptor_names:
            selected_ligands = _checkbox(questionary, f"Select ligands for {receptor_name}", ligand_pool)
            if selected_ligands:
                curated_mapping[receptor_name] = selected_ligands

    round_id: Optional[str] = None
    iterative = workflow_mode in {"iterative_only", "iterative_freeze"}
    freeze = workflow_mode in {"iterative_freeze", "freeze_existing"}
    if iterative or freeze:
        round_default = latest_round or "round_001"
        round_label = "Pair-curation round to update" if iterative else "Saved pair-curation round to freeze"
        round_id = _text(questionary, round_label, default=round_default, required=True)

    prompt_protein_aliases = _confirm(
        questionary,
        "Prompt for persistent protein display aliases before writing pair metadata?",
        default=True,
    )
    prompt_ligand_aliases = _confirm(
        questionary,
        "Prompt for persistent ligand display aliases before writing pair metadata?",
        default=False,
    )

    result = run_prepare_pairlist(
        project_dir=str(project_root),
        mode=mode,
        prepared_proteins=str(layout["prepared_proteins"]),
        prepared_ligands=str(layout["prepared_ligands"]),
        excel_path=str(excel_file),
        curated_receptors=curated_receptors,
        curated_ligands=curated_ligands,
        curated_mapping=curated_mapping,
        iterative=iterative,
        round_id=round_id,
        freeze=freeze,
        prompt_protein_aliases=prompt_protein_aliases,
        prompt_ligand_aliases=prompt_ligand_aliases,
    )
    print(f"Pairlist rows: {result.outputs.get('pair_count', 0)}")
    if result.outputs.get("protein_alias_file"):
        print(f"Protein aliases: {result.outputs.get('protein_alias_file')}")
    if result.outputs.get("ligand_alias_file"):
        print(f"Ligand aliases: {result.outputs.get('ligand_alias_file')}")
    if result.notes:
        print("Pairlist notes:")
        for note in result.notes:
            print(f" - {note}")


def _materialize_project(questionary, project_root: Path) -> None:
    _header("Prepare Docking Folders")
    layout = _project_layout(project_root)
    prepared_proteins = list(
        path
        for path in Path(layout["prepared_proteins"]).glob("*")
        if path.is_file() and path.suffix.lower() in {".pdb", ".pdbqt"}
    )
    prepared_ligands = list(
        path
        for path in Path(layout["prepared_ligands"]).glob("*")
        if path.is_file() and path.suffix.lower() in {".pdb", ".pdbqt", ".sdf", ".mol2"}
    )
    if not prepared_proteins or not prepared_ligands:
        print("❌ Cannot prepare docking folders yet: prepared assets are missing.")
        if not prepared_proteins:
            print(f"   Missing prepared proteins in: {layout['prepared_proteins']}")
            print("   Expected at least one .pdb or .pdbqt file.")
        if not prepared_ligands:
            print(f"   Missing prepared ligands in: {layout['prepared_ligands']}")
            print("   Expected at least one .pdbqt/.sdf/.mol2/.pdb file.")
        print("   Run the 'Prepare raw proteins and ligands' stage first, then retry this step.")
        return

    manifest = load_manifest(project_root)
    engine_selection = _select(
        questionary,
        "Prepare docking folders for which panels?",
        ENGINE_CHOICES,
        default="all" if set(manifest.get("enabled_panels", [])) == {"gnina", "vina", "smina", "autodock4"} else (manifest.get("enabled_panels", ["gnina"])[0] if manifest.get("enabled_panels") else "all"),
    )
    argv = [
        "--prepared-proteins",
        str(layout["prepared_proteins"]),
        "--prepared-ligands",
        str(layout["prepared_ligands"]),
        "--excel",
        str(layout["raw_proteins"] / "multi_pdb_analysis.xlsx"),
        "--output",
        str(project_root),
        "--layout-profile",
        "docking_legacy",
        "--asset-mode",
        "copy",
        "--engines",
        ",".join(_resolve_engine_selection(engine_selection)),
        "--pairlist-file",
        str(pairlist_path(project_root, "docking_legacy")),
    ]
    if pair_intent_path(project_root, "docking_legacy").exists():
        argv.extend(["--pair-intent", str(pair_intent_path(project_root, "docking_legacy"))])
    result = run_prepare_project(argv)
    print(f"Project materialization status: {result.status}")


def _run_docking_stage(questionary, project_root: Path) -> None:
    _header("Docking Execution")
    manifest = load_manifest(project_root)
    layout = _project_layout(project_root)
    engines = _checkbox(
        questionary,
        "Select engines to run (space to toggle, enter to confirm)",
        manifest.get("enabled_panels", ["gnina", "vina", "smina", "autodock4"]),
    )
    if not engines:
        print("No engines selected.")
        return
    execution_environment = _select(
        questionary,
        "Execution environment",
        [
            ("local_cpu", "Local CPU"),
            ("local_gpu", "Local GPU"),
            ("conda_env", "Conda environment"),
            ("container", "Containerized"),
            ("remote_hpc", "Remote/HPC planning mode"),
        ],
        default="local_cpu",
    )
    execution_workdir = _text(
        questionary,
        "Execution working directory (optional)",
        default=str(project_root),
        required=False,
    )
    prerequisites_dir = _text(
        questionary,
        "Prerequisites directory (optional)",
        default="",
        required=False,
    )
    required_files_dir = _text(
        questionary,
        "Required-files directory (optional)",
        default="",
        required=False,
    )
    parameter_mode = _select(
        questionary,
        "Docking parameter mode",
        [("basic", "Basic (preset-driven)"), ("advanced", "Advanced (manual controls)")],
        default="basic",
    )
    parameter_preset = "balanced"
    exhaustiveness = 16
    num_modes = 20
    seed = ""
    box_scale = "1.0"
    box_padding = "0.0"
    if parameter_mode == "basic":
        parameter_preset = _select(
            questionary,
            "Basic parameter preset",
            [
                ("screening_fast", "Screening Fast (8/10)"),
                ("balanced", "Balanced (16/20)"),
                ("exhaustive", "Exhaustive (32/40)"),
            ],
            default="balanced",
        )
    else:
        exhaustiveness = _int_input(questionary, "Exhaustiveness", default=16, minimum=1)
        num_modes = _int_input(questionary, "Number of modes", default=20, minimum=1)
        seed = _text(questionary, "Seed (optional, Enter skips)", default="", required=False)
        box_scale = str(_float_input(questionary, "Box scale (>0)", default=1.0))
        box_padding = str(_float_input(questionary, "Box padding (>=0)", default=0.0))

    ligand_qc_gate_enabled = _confirm(
        questionary,
        "Enable ligand QC gate before docking?",
        default=True,
    )
    ligand_admet_filters_enabled = True
    ligand_qc_max_lipinski_violations = 1
    ligand_qc_max_molecular_weight = 650.0
    ligand_qc_max_logp = 6.0
    ligand_qc_max_tpsa = 180.0
    ligand_qc_max_rotatable_bonds = 15
    ligand_qc_max_formal_charge_abs = 2
    ligand_qc_min_heavy_atom_count = 6
    ligand_qc_block_pains = True
    ligand_qc_block_brenk = True
    ligand_qc_block_reactive = True
    if ligand_qc_gate_enabled:
        ligand_admet_filters_enabled = _confirm(
            questionary,
            "Enable ADMET filters inside ligand QC gate?",
            default=True,
        )
        if ligand_admet_filters_enabled and _confirm(
            questionary,
            "Customize ligand ADMET thresholds?",
            default=False,
        ):
            ligand_qc_max_lipinski_violations = _int_input(
                questionary,
                "Maximum Lipinski violations",
                default=1,
                minimum=0,
            )
            ligand_qc_max_molecular_weight = _float_input(
                questionary,
                "Maximum molecular weight (Da)",
                default=650.0,
            )
            ligand_qc_max_logp = _float_input(
                questionary,
                "Maximum LogP",
                default=6.0,
            )
            ligand_qc_max_tpsa = _float_input(
                questionary,
                "Maximum TPSA",
                default=180.0,
            )
            ligand_qc_max_rotatable_bonds = _int_input(
                questionary,
                "Maximum rotatable bonds",
                default=15,
                minimum=0,
            )
            ligand_qc_max_formal_charge_abs = _int_input(
                questionary,
                "Maximum absolute formal charge",
                default=2,
                minimum=0,
            )
            ligand_qc_min_heavy_atom_count = _int_input(
                questionary,
                "Minimum heavy-atom count",
                default=6,
                minimum=0,
            )
            ligand_qc_block_pains = _confirm(
                questionary,
                "Block ligands with PAINS alerts?",
                default=True,
            )
            ligand_qc_block_brenk = _confirm(
                questionary,
                "Block ligands with Brenk alerts?",
                default=True,
            )
            ligand_qc_block_reactive = _confirm(
                questionary,
                "Block ligands with reactive SMARTS alerts?",
                default=True,
            )
    receptor_qc_gate_enabled = _confirm(
        questionary,
        "Enable receptor QC gate before docking?",
        default=True,
    )
    receptor_qc_min_atom_count = 100
    receptor_qc_min_heavy_atom_count = 60
    receptor_qc_min_chain_count = 1
    receptor_qc_max_coordinate_span = 500.0
    if receptor_qc_gate_enabled:
        receptor_qc_min_atom_count = _int_input(
            questionary,
            "Receptor QC minimum atom count",
            default=100,
            minimum=1,
        )
        receptor_qc_min_heavy_atom_count = _int_input(
            questionary,
            "Receptor QC minimum heavy-atom count",
            default=60,
            minimum=1,
        )
        receptor_qc_min_chain_count = _int_input(
            questionary,
            "Receptor QC minimum chain count",
            default=1,
            minimum=1,
        )
        receptor_qc_max_coordinate_span = _float_input(
            questionary,
            "Receptor QC maximum coordinate span (A)",
            default=500.0,
        )

    shared_conda_env = ""
    if execution_environment == "conda_env":
        shared_conda_env = _text(
            questionary,
            "Shared conda environment for Vina/Smina (optional)",
            default="",
            required=False,
        )
    container_image = ""
    if execution_environment == "container":
        container_image = _text(
            questionary,
            "Container image path (optional fallback for GNINA)",
            default="",
            required=False,
        )

    dry_run = _confirm(questionary, "Dry run only?", default=True)
    run_in_background = _confirm(
        questionary,
        "Run docking in background so you can continue using other panels?",
        default=True,
    )
    favorite_engine = _select(
        questionary,
        "Default engine for downstream analysis after docking (optional)",
        [("", "Skip")] + [(engine, engine.upper()) for engine in engines],
        default="",
    )
    argv = ["--project-dir", str(project_root), "--engines", ",".join(engines)]
    argv.extend(["--execution-environment", execution_environment])
    if execution_workdir:
        argv.extend(["--execution-workdir", execution_workdir])
    if prerequisites_dir:
        argv.extend(["--prerequisites-dir", prerequisites_dir])
    if required_files_dir:
        argv.extend(["--required-files-dir", required_files_dir])
    argv.extend(["--parameter-mode", parameter_mode, "--parameter-preset", parameter_preset])
    if parameter_mode == "advanced":
        argv.extend(
            [
                "--exhaustiveness",
                str(exhaustiveness),
                "--num-modes",
                str(num_modes),
                "--box-scale",
                str(box_scale),
                "--box-padding",
                str(box_padding),
            ]
        )
        if seed.strip():
            argv.extend(["--seed", seed.strip()])
    if shared_conda_env:
        argv.extend(["--shared-conda-env", shared_conda_env])
    if container_image:
        argv.extend(["--container-image", container_image])
    if not ligand_qc_gate_enabled:
        argv.append("--no-ligand-qc-gate")
    if not ligand_admet_filters_enabled:
        argv.append("--no-ligand-admet-filters")
    argv.extend(
        [
            "--ligand-qc-max-lipinski-violations",
            str(ligand_qc_max_lipinski_violations),
            "--ligand-qc-max-molecular-weight",
            str(ligand_qc_max_molecular_weight),
            "--ligand-qc-max-logp",
            str(ligand_qc_max_logp),
            "--ligand-qc-max-tpsa",
            str(ligand_qc_max_tpsa),
            "--ligand-qc-max-rotatable-bonds",
            str(ligand_qc_max_rotatable_bonds),
            "--ligand-qc-max-formal-charge-abs",
            str(ligand_qc_max_formal_charge_abs),
            "--ligand-qc-min-heavy-atom-count",
            str(ligand_qc_min_heavy_atom_count),
        ]
    )
    if not ligand_qc_block_pains:
        argv.append("--ligand-qc-allow-pains")
    if not ligand_qc_block_brenk:
        argv.append("--ligand-qc-allow-brenk")
    if not ligand_qc_block_reactive:
        argv.append("--ligand-qc-allow-reactive")
    if not receptor_qc_gate_enabled:
        argv.append("--no-receptor-qc-gate")
    argv.extend(
        [
            "--receptor-qc-min-atom-count",
            str(receptor_qc_min_atom_count),
            "--receptor-qc-min-heavy-atom-count",
            str(receptor_qc_min_heavy_atom_count),
            "--receptor-qc-min-chain-count",
            str(receptor_qc_min_chain_count),
            "--receptor-qc-max-coordinate-span",
            str(receptor_qc_max_coordinate_span),
        ]
    )
    if dry_run:
        argv.append("--dry-run")
    if favorite_engine:
        argv.extend(["--favorite-engine", favorite_engine])
    if run_in_background:
        task_id = _TASK_MANAGER.submit(
            root=project_root,
            label="Docking execution",
            target="dock.run",
            fn=run_docking,
            argv=argv,
        )
        print(f"Docking queued in background (task id: {task_id}).")
        return

    result = run_docking(argv)
    print(f"Docking status: {result.status}")


def _run_hpc_stage(questionary, project_root: Path) -> None:
    _header("HPC Deployment")
    manifest = load_manifest(project_root)
    action = _select(
        questionary,
        "What do you want to do on the HPC?",
        [
            ("deploy", "Generate portable HPC deployment bundle"),
            ("sync", "Sync project from this Mac to the HPC"),
            ("submit", "Submit an already synced deployment on the HPC"),
            ("deploy_sync", "Generate deployment and sync project"),
            ("deploy_sync_submit", "Generate deployment, sync project, and submit jobs"),
        ],
        default="deploy_sync",
    )
    run_in_background = False
    if action in {"deploy", "sync", "submit"}:
        run_in_background = _confirm(
            questionary,
            "Run this HPC action in background so you can keep navigating?",
            default=True,
        )

    selected_system: Dict[str, str] = {}
    systems, systems_source = _load_optional_ssh_systems(project_root)
    if systems:
        system_default = systems[0]["name"]
        selected_system_name = _select(
            questionary,
            f"SSH system preset ({systems_source})",
            [
                (
                    entry["name"],
                    (
                        f"{entry['name']} "
                        f"[{entry['kind'] or 'ssh'}]"
                        f"{' -> ' + entry['ssh_target'] if entry.get('ssh_target') else ''}"
                    ),
                )
                for entry in systems
            ],
            default=system_default,
        )
        selected_system = next((entry for entry in systems if entry["name"] == selected_system_name), {})

    profile_default = str(
        selected_system.get("hpc_profile")
        or manifest.get("hpc_profile", {}).get("name")
        or "alex-bibalex"
    )
    profile_file = str(selected_system.get("hpc_profile_file", "") or "").strip()
    profile_name = _text(questionary, "HPC profile name", default=profile_default, required=True)

    try:
        hpc_profile = load_hpc_profile(
            project_root,
            profile_name=profile_name,
            profile_file=profile_file,
        )
    except ValueError as exc:
        print(f"❌ {exc}")
        return

    remote_project_default = resolve_remote_project_dir(project_root, hpc_profile) or ""
    ssh_target_default = resolve_ssh_target(hpc_profile) or ""
    if selected_system:
        remote_project_default = str(selected_system.get("remote_project_dir") or remote_project_default)
        ssh_target_default = str(selected_system.get("ssh_target") or ssh_target_default)
    remote_project_dir = _text(
        questionary,
        "Remote project directory on the HPC",
        default=remote_project_default,
        required=True,
    )
    ssh_target = _text(
        questionary,
        "SSH target (user@host)",
        default=ssh_target_default,
        required=action in {"sync", "submit", "deploy_sync", "deploy_sync_submit"},
    )
    if ssh_target and selected_system.get("ssh_target"):
        ok, probe_message = _probe_ssh_connectivity(ssh_target)
        channel = str(selected_system.get("kind") or "").strip().lower() or "ssh"
        channel_label = "NMRbox SSH" if channel == "nmrbox" else "SSH"
        if ok:
            print(f"✅ {channel_label} connectivity probe passed for {ssh_target}")
        else:
            print(f"⚠️ {channel_label} connectivity probe failed for {ssh_target}: {probe_message}")
    profile_cli_args = (
        ["--hpc-profile-file", profile_file]
        if profile_file
        else ["--hpc-profile", profile_name]
    )
    condor_fit_cli_args: List[str] = []
    condor_fit_payload: Dict[str, object] = {}
    scheduler_name = str((hpc_profile or {}).get("scheduler") or "slurm").strip().lower()
    selected_system_kind = str(selected_system.get("kind") or "").strip().lower()
    deploy_like_action = action in {"deploy", "deploy_sync", "deploy_sync_submit"}
    if selected_system_kind == "nmrbox" and scheduler_name == "condor":
        auto_fit = _confirm(
            questionary,
            "Auto-fit Condor resources using live NMRbox pool health?",
            default=True,
        )
        if auto_fit:
            ok_fit, fit_payload, fit_message = _probe_nmrbox_condor_fit(ssh_target)
            if ok_fit:
                print(
                    "✅ NMRbox Condor fit: "
                    f"idle_slots={fit_payload.get('idle_slots_total', 0)} "
                    f"idle_gpu_slots={fit_payload.get('idle_gpu_slots_total', 0)} "
                    f"cpus={fit_payload.get('recommended_cpus')} "
                    f"memory={fit_payload.get('recommended_memory')} "
                    f"disk={fit_payload.get('recommended_disk')}"
                )
                condor_fit_cli_args.extend(
                    [
                        "--condor-cpus",
                        str(fit_payload.get("recommended_cpus", 8)),
                        "--condor-memory",
                        str(fit_payload.get("recommended_memory", "32GB")),
                        "--condor-disk",
                        str(fit_payload.get("recommended_disk", "20GB")),
                    ]
                )
                condor_fit_payload = dict(fit_payload)
            else:
                print(f"⚠️ NMRbox Condor fit probe failed: {fit_message}")

    enabled_engines = manifest.get("enabled_panels") or manifest.get("engines", ["gnina", "vina", "smina", "autodock4"])
    if not enabled_engines:
        enabled_engines = ["gnina", "vina", "smina", "autodock4"]
    deploy_engines = list(enabled_engines)
    submit_engines = list(enabled_engines)
    inferred_round, inferred_mode = infer_round_mode_from_deployment_root(str(manifest.get("deployment_root", "") or ""))
    round_default = inferred_round or "hpc_screen_001"
    mode_default = inferred_mode or "screen"
    selected_mode = mode_default
    selected_round = round_default
    rerun_manifest_path = ""
    allow_full_exhaustive = False

    if action in {"deploy", "deploy_sync", "deploy_sync_submit"}:
        deploy_engines = _checkbox(
            questionary,
            "Select engines for HPC deployment (space to toggle, enter to confirm)",
            list(enabled_engines),
            preselected=list(enabled_engines),
        )
        if not deploy_engines:
            use_all = _confirm(
                questionary,
                "No engines selected for deployment. Use all enabled engines instead?",
                default=True,
            )
            if use_all:
                deploy_engines = list(enabled_engines)
            else:
                print("No engines selected for deployment.")
                return
        selected_mode = _select(
            questionary,
            "Deployment mode",
            [
                ("screen", "Screening (broad first pass)"),
                ("exhaustive", "Exhaustive rerun (normally from comparative rerun manifest)"),
            ],
            default=mode_default,
        )
        selected_round = _text(
            questionary,
            "Deployment round label",
            default=_suggest_hpc_round_label(project_root, selected_mode, inferred_round, inferred_mode),
            required=True,
        )
        if selected_mode == "exhaustive":
            rerun_manifest_path = _optional_text(
                questionary,
                "Rerun manifest CSV for exhaustive mode",
                default=_default_rerun_manifest(project_root, manifest),
            )
            if not rerun_manifest_path:
                print("No rerun manifest was selected for exhaustive mode.")
                allow_full_exhaustive = _confirm(
                    questionary,
                    "Continue with a full-pair exhaustive deployment anyway? This is usually the wrong first pass.",
                    default=False,
                )
                if not allow_full_exhaustive:
                    print("Aborted. Run screening first or generate a rerun manifest from comparative analysis.")
                    return

    if action in {"submit", "deploy_sync_submit"}:
        submit_engines = _checkbox(
            questionary,
            "Select engines to submit on the HPC (space to toggle, enter to confirm)",
            list(enabled_engines),
            preselected=list(enabled_engines),
        )
        if not submit_engines:
            use_all = _confirm(
                questionary,
                "No engines selected for submission. Use all enabled engines instead?",
                default=True,
            )
            if use_all:
                submit_engines = list(enabled_engines)
            else:
                print("No engines selected for submission.")
                return
        if action == "submit":
            selected_mode = _select(
                questionary,
                "Deployment mode to submit",
                [
                    ("screen", "Screening (broad first pass)"),
                    ("exhaustive", "Exhaustive rerun (normally from comparative rerun manifest)"),
                ],
                default=mode_default,
            )
            selected_round = _text(
                questionary,
                "Deployment round label",
                default=inferred_round if inferred_round and inferred_mode == selected_mode else f"hpc_{selected_mode}_001",
                required=True,
            )
        if selected_mode == "exhaustive":
            deployment_manifest = _load_local_deployment_manifest(project_root, selected_round, selected_mode)
            rerun_manifest_record = str(deployment_manifest.get("rerun_manifest_file", "") or "").strip()
            if not rerun_manifest_record:
                allow_full_exhaustive = _confirm(
                    questionary,
                    "This exhaustive bundle does not record a rerun manifest. Submit a full-pair exhaustive campaign anyway?",
                    default=False,
                )
                if not allow_full_exhaustive:
                    print("Aborted. Submit a screening bundle or regenerate exhaustive from a rerun manifest.")
                    return

    if condor_fit_payload and deploy_like_action and "gnina" in deploy_engines:
        if set(deploy_engines) == {"gnina"}:
            if int(condor_fit_payload.get("recommended_gpus", 0) or 0) > 0:
                condor_fit_cli_args.extend(["--condor-gpus", "1"])
                gpu_requirements = str(condor_fit_payload.get("recommended_requirements", "") or "").strip()
                if gpu_requirements:
                    condor_fit_cli_args.extend(["--condor-requirements", gpu_requirements])
                print(
                    "✅ GNINA-only deploy: applied GPU-aware Condor fit "
                    f"({gpu_requirements or 'GPUs >= 1'})"
                )
            else:
                print("⚠️ No idle GPU slots detected now; keeping GPU requirement unchanged.")
        else:
            print(
                "ℹ️ Mixed-engine deploy detected; applied CPU/memory/disk fit globally. "
                "GPU requirement was not forced globally to avoid constraining CPU-only engines."
            )

    if action in {"sync", "deploy_sync", "deploy_sync_submit"}:
        delete_remote = _confirm(questionary, "Delete remote files that no longer exist locally during sync?", default=False)
    else:
        delete_remote = False

    ligand_qc_gate_enabled = True
    ligand_admet_filters_enabled = True
    receptor_qc_gate_enabled = True
    if action in {"deploy", "submit", "deploy_sync", "deploy_sync_submit"}:
        ligand_qc_gate_enabled = _confirm(
            questionary,
            "Enforce ligand QC gate before HPC action?",
            default=True,
        )
        if ligand_qc_gate_enabled:
            ligand_admet_filters_enabled = _confirm(
                questionary,
                "Enable ADMET filters inside ligand QC gate for HPC action?",
                default=False,
            )
        receptor_qc_gate_enabled = _confirm(
            questionary,
            "Enforce receptor QC gate before HPC action?",
            default=True,
        )

    if action in {"deploy", "submit", "deploy_sync", "deploy_sync_submit"} and (ligand_qc_gate_enabled or receptor_qc_gate_enabled):
        print("Running ligand/receptor preflight validation before any HPC action...")
        df = load_pairlist(project_root)
        rows = [
            PairlistRow(
                receptor=str(row["receptor"]),
                site_id=str(row["site_id"]),
                ligand=str(row["ligand"]),
                center_x=float(row["center_x"]),
                center_y=float(row["center_y"]),
                center_z=float(row["center_z"]),
                size_x=float(row["size_x"]),
                size_y=float(row["size_y"]),
                size_z=float(row["size_z"]),
            )
            for _, row in df.iterrows()
        ]
        metadata_dir = Path(load_manifest(project_root).get("metadata_dir") or project_root / "metadata")
        if ligand_qc_gate_enabled:
            ligand_report_path = metadata_dir / "ligand_validation_report.json"
            ligand_audit = audit_project_ligands(
                project_root,
                rows,
                enable_admet_filters=ligand_admet_filters_enabled,
                report_path=ligand_report_path,
            )
            if ligand_audit.get("issue_count"):
                print("❌ HPC preflight blocked")
                print(f"   Ligand report: {ligand_audit.get('report_path', ligand_report_path)}")
                print(f"   Invalid ligands: {', '.join(ligand_audit.get('affected_ligands', []))}")
                print("   Repair or remove those ligands before deploying or submitting jobs.")
                return
        if receptor_qc_gate_enabled:
            receptor_report_path = metadata_dir / "receptor_validation_report.json"
            receptor_audit = audit_project_receptors(
                project_root=project_root,
                rows=rows,
                min_atom_count=100,
                min_heavy_atom_count=60,
                min_chain_count=1,
                max_coordinate_span=500.0,
                report_path=receptor_report_path,
            )
            if receptor_audit.get("issue_count"):
                print("❌ HPC preflight blocked")
                print(f"   Receptor report: {receptor_audit.get('report_path', receptor_report_path)}")
                print(f"   Invalid receptors: {', '.join(receptor_audit.get('affected_receptors', []))}")
                thresholds = receptor_audit.get("thresholds") or {}
                if thresholds:
                    print(
                        "   Thresholds: "
                        f"min_atom_count={thresholds.get('min_atom_count')} "
                        f"min_heavy_atom_count={thresholds.get('min_heavy_atom_count')} "
                        f"min_chain_count={thresholds.get('min_chain_count')} "
                        f"max_coordinate_span={thresholds.get('max_coordinate_span')}"
                    )
                print("   Repair receptor preparation artifacts or relax QC thresholds in CLI mode before HPC submission.")
                return
        print("✅ Ligand/receptor preflight passed")
    elif action in {"deploy", "submit", "deploy_sync", "deploy_sync_submit"}:
        print("⚠️ HPC preflight skipped because both ligand and receptor QC gates were disabled for this action.")

    if action == "deploy":
        deploy_argv = [
            "--project-dir",
            str(project_root),
            "--engines",
            ",".join(deploy_engines),
            "--mode",
            selected_mode,
            "--round",
            selected_round,
            *profile_cli_args,
            *condor_fit_cli_args,
            "--remote-project-dir",
            remote_project_dir,
            *(["--from-rerun-manifest", rerun_manifest_path] if rerun_manifest_path else []),
            *(["--allow-full-exhaustive"] if allow_full_exhaustive else []),
            *(["--no-ligand-qc-gate"] if not ligand_qc_gate_enabled else []),
            *(["--no-ligand-admet-filters"] if ligand_qc_gate_enabled and not ligand_admet_filters_enabled else []),
            *(["--no-receptor-qc-gate"] if not receptor_qc_gate_enabled else []),
        ]
        if run_in_background:
            task_id = _TASK_MANAGER.submit(
                root=project_root,
                label="HPC deploy",
                target="dock.deploy",
                fn=run_docking_deploy,
                argv=deploy_argv,
            )
            print(f"HPC deploy queued in background (task id: {task_id}).")
            return
        result = run_docking_deploy(deploy_argv)
        print(f"HPC deployment generation: {result.status}")
        return

    if action == "sync":
        sync_argv = [
            "--project-dir",
            str(project_root),
            *profile_cli_args,
            "--ssh-target",
            ssh_target,
            "--remote-project-dir",
            remote_project_dir,
            *(["--delete"] if delete_remote else []),
        ]
        if run_in_background:
            task_id = _TASK_MANAGER.submit(
                root=project_root,
                label="HPC sync",
                target="dock.sync",
                fn=run_docking_sync,
                argv=sync_argv,
            )
            print(f"HPC sync queued in background (task id: {task_id}).")
            return
        result = run_docking_sync(sync_argv)
        print(f"HPC project sync: {result.status}")
        return

    if action == "submit":
        confirm_submit = _confirm(
            questionary,
            (
                f"Submit {len(submit_engines)} engine deployment(s) on the HPC now? "
                "Each engine submits one job script; the active profile decides whether that's a single sequential job or an array."
            ),
            default=False,
        )
        submit_argv = [
            "--project-dir",
            str(project_root),
            "--engines",
            ",".join(submit_engines),
            "--round",
            selected_round,
            "--mode",
            selected_mode,
            *profile_cli_args,
            *condor_fit_cli_args,
            "--ssh-target",
            ssh_target,
            "--remote-project-dir",
            remote_project_dir,
            *(["--allow-full-exhaustive"] if allow_full_exhaustive else []),
            *(["--no-ligand-qc-gate"] if not ligand_qc_gate_enabled else []),
            *(["--no-ligand-admet-filters"] if ligand_qc_gate_enabled and not ligand_admet_filters_enabled else []),
            *(["--no-receptor-qc-gate"] if not receptor_qc_gate_enabled else []),
            *(["--dry-run"] if not confirm_submit else []),
        ]
        if run_in_background:
            task_id = _TASK_MANAGER.submit(
                root=project_root,
                label="HPC submit",
                target="dock.submit",
                fn=run_docking_submit,
                argv=submit_argv,
            )
            print(f"HPC submit queued in background (task id: {task_id}).")
            return
        result = run_docking_submit(submit_argv)
        print(f"HPC submit: {result.status}")
        return

    deploy_result = run_docking_deploy(
        [
            "--project-dir",
            str(project_root),
            "--engines",
            ",".join(deploy_engines),
            "--mode",
            selected_mode,
            "--round",
            selected_round,
            *profile_cli_args,
            *condor_fit_cli_args,
            "--remote-project-dir",
            remote_project_dir,
            *(["--from-rerun-manifest", rerun_manifest_path] if rerun_manifest_path else []),
            *(["--allow-full-exhaustive"] if allow_full_exhaustive else []),
            *(["--no-ligand-qc-gate"] if not ligand_qc_gate_enabled else []),
            *(["--no-ligand-admet-filters"] if ligand_qc_gate_enabled and not ligand_admet_filters_enabled else []),
            *(["--no-receptor-qc-gate"] if not receptor_qc_gate_enabled else []),
        ]
    )
    print(f"HPC deployment generation: {deploy_result.status}")
    if deploy_result.status != "completed":
        return

    sync_result = run_docking_sync(
        [
            "--project-dir",
            str(project_root),
            *profile_cli_args,
            "--ssh-target",
            ssh_target,
            "--remote-project-dir",
            remote_project_dir,
            *(["--delete"] if delete_remote else []),
        ]
    )
    print(f"HPC project sync: {sync_result.status}")
    if sync_result.status != "completed":
        return

    if action == "deploy_sync":
        return

    confirm_submit = _confirm(
        questionary,
        (
            f"Submit {len(submit_engines)} engine deployment(s) on the HPC now? "
            "Each engine submits one job script; the active profile decides whether that's a single sequential job or an array."
        ),
        default=False,
    )
    submit_result = run_docking_submit(
        [
            "--project-dir",
            str(project_root),
            "--engines",
            ",".join(submit_engines),
            "--round",
            selected_round,
            "--mode",
            selected_mode,
            *profile_cli_args,
            "--ssh-target",
            ssh_target,
            "--remote-project-dir",
            remote_project_dir,
            *(["--allow-full-exhaustive"] if allow_full_exhaustive else []),
            *(["--no-ligand-qc-gate"] if not ligand_qc_gate_enabled else []),
            *(["--no-ligand-admet-filters"] if ligand_qc_gate_enabled and not ligand_admet_filters_enabled else []),
            *(["--no-receptor-qc-gate"] if not receptor_qc_gate_enabled else []),
            *(["--dry-run"] if not confirm_submit else []),
        ]
    )
    print(f"HPC submit: {submit_result.status}")


def _edit_display_names(questionary, project_root: Path) -> None:
    _header("Edit Protein And Ligand Display Names")
    from post_docking_analysis.ligand_naming import (
        build_ligand_name_mapping,
        ligand_mapping_to_dataframe,
    )
    from post_docking_analysis.protein_naming import (
        build_protein_name_mapping,
        extract_pdb_code,
        mapping_to_dataframe,
        resolve_protein_display_name,
    )

    profile = detect_layout_profile(project_root)
    naming_dir = _shared_naming_dir(project_root)
    protein_override_file = naming_dir / "protein_name_overrides.csv"
    ligand_override_file = naming_dir / "ligand_name_overrides.csv"

    pairlist_file = pairlist_path(project_root, profile)
    pairlist_df = pd.DataFrame()
    if pairlist_file.exists():
        try:
            pairlist_df = pd.read_csv(pairlist_file)
        except Exception:
            pairlist_df = pd.DataFrame()

    receptors_dir = shared_receptors_dir(project_root, profile)
    receptor_files = sorted(path for path in receptors_dir.glob("*") if path.is_file()) if receptors_dir.exists() else []

    protein_map = build_protein_name_mapping(
        receptor_files,
        pairlist_df=pairlist_df,
        overrides_file=protein_override_file if protein_override_file.exists() else None,
    )
    ligand_identifiers: List[str] = []
    if not pairlist_df.empty:
        for column in ("ligand", "ligand_name", "cocrystal_ligand_name"):
            if column in pairlist_df.columns:
                ligand_identifiers.extend(pairlist_df[column].dropna().astype(str).tolist())
    ligand_map = build_ligand_name_mapping(
        ligand_identifiers,
        pairlist_df=pairlist_df,
        overrides_file=ligand_override_file if ligand_override_file.exists() else None,
    )

    code_map = protein_map.get("pdb_code_to_name", {})
    receptor_map = protein_map.get("receptor_to_name", {})
    prompt_targets: List[tuple[str, str, str]] = []
    seen_codes = set()
    for receptor_file in receptor_files:
        receptor_name = str(receptor_file.name)
        receptor_stem = str(receptor_file.stem)
        code = extract_pdb_code(receptor_name) or extract_pdb_code(receptor_stem)
        if not code or code in seen_codes:
            continue
        seen_codes.add(code)
        current = str(code_map.get(code, "")).strip() or resolve_protein_display_name(receptor_name, protein_map)
        prompt_targets.append(("code", code, current or f"Protein {code}"))

    if not prompt_targets:
        for code in sorted(code_map):
            prompt_targets.append(("code", code, str(code_map.get(code, "")).strip() or f"Protein {code}"))

    if prompt_targets:
        print("Edit protein display names (Enter keeps current).")
        for target_type, target_id, current in prompt_targets:
            updated = _optional_text(questionary, f"Protein name for PDB {target_id} [{current}]", default="")
            final_name = updated or current
            if target_type == "code":
                code_map[target_id] = final_name
                for receptor_key in list(receptor_map.keys()):
                    key_code = extract_pdb_code(receptor_key)
                    if key_code and key_code == target_id:
                        receptor_map[receptor_key] = final_name

    ligand_prompts: List[tuple[str, str]] = []
    seen_stems = set()
    for key in sorted(ligand_map):
        stem = Path(str(key)).stem
        if key != stem or stem in seen_stems:
            continue
        seen_stems.add(stem)
        ligand_prompts.append((stem, str(ligand_map.get(stem, "")).strip() or stem))

    if ligand_prompts:
        print("Edit ligand display names (Enter keeps current).")
        for ligand_id, current in ligand_prompts:
            updated = _optional_text(questionary, f"Ligand name for {ligand_id} [{current}]", default="")
            final_name = updated or current
            ligand_map[ligand_id] = final_name
            ligand_map[f"{ligand_id}.pdbqt"] = final_name

    mapping_to_dataframe(protein_map).to_csv(naming_dir / "protein_name_mapping.csv", index=False)
    ligand_mapping_to_dataframe(ligand_map).to_csv(naming_dir / "ligand_name_mapping.csv", index=False)

    protein_rows = [
        {"pdb_code": code, "receptor": "", "display_name": display}
        for code, display in sorted(code_map.items())
    ]
    pd.DataFrame(protein_rows, columns=["pdb_code", "receptor", "display_name"]).to_csv(
        protein_override_file,
        index=False,
    )

    ligand_rows = []
    seen_ligands = set()
    for key, display in sorted(ligand_map.items()):
        ligand = Path(str(key)).stem
        if ligand in seen_ligands:
            continue
        seen_ligands.add(ligand)
        ligand_rows.append({"ligand": ligand, "display_name": display})
    pd.DataFrame(ligand_rows, columns=["ligand", "display_name"]).to_csv(
        ligand_override_file,
        index=False,
    )

    print(f"Saved protein names: {naming_dir / 'protein_name_mapping.csv'}")
    print(f"Saved ligand names:  {naming_dir / 'ligand_name_mapping.csv'}")


def _analysis_progress_file(session_root: Path) -> Path:
    return session_root / "analysis_progress.json"


def _read_analysis_progress(progress_file: Path) -> Dict[str, object]:
    if not progress_file.exists():
        return {}
    try:
        payload = json.loads(progress_file.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _completed_target_status_map(completed_targets: List[Dict[str, str]]) -> Dict[str, str]:
    status_map: Dict[str, str] = {}
    for row in completed_targets:
        target = str(row.get("target", "")).strip()
        if not target:
            continue
        status_map[target] = str(row.get("status", "")).strip()
    return status_map


def _write_analysis_progress(
    progress_file: Path,
    session_root: Path,
    selected_targets: List[str],
    completed_targets: List[Dict[str, str]],
    running_target: str = "",
    status: str = "pending",
    note: str = "",
) -> None:
    payload = {
        "updated_at": datetime.now().isoformat(timespec="seconds"),
        "session_root": str(session_root),
        "total_targets": len(selected_targets),
        "running_target": running_target,
        "status": status,
        "note": note,
        "selected_targets": selected_targets,
        "completed_targets": completed_targets,
        "completed_count": len(completed_targets),
        "remaining_count": max(len(selected_targets) - len(completed_targets), 0),
    }
    progress_file.parent.mkdir(parents=True, exist_ok=True)
    progress_file.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _show_analysis_progress(progress_file: Path) -> None:
    if not progress_file.exists():
        print("No analysis progress has been recorded for this session yet.")
        return
    payload = _read_analysis_progress(progress_file)
    if not payload:
        print(f"Could not read progress file: {progress_file}")
        return

    total = int(payload.get("total_targets", 0))
    completed = int(payload.get("completed_count", 0))
    running = str(payload.get("running_target", "")).strip()
    status = str(payload.get("status", "unknown"))
    updated = str(payload.get("updated_at", ""))
    note = str(payload.get("note", "")).strip()
    completed_rows = payload.get("completed_targets") or []

    percent = int(round((completed / total) * 100)) if total > 0 else 0
    print(f"Progress: {completed}/{total} completed ({percent}%) | status={status} | updated={updated}")
    if running:
        print(f"Current target: {ANALYSIS_LABELS.get(running, running)}")
    if note:
        print(f"Note: {note}")
    if completed_rows:
        print("Completed targets:")
        for row in completed_rows:
            target = str(row.get("target", ""))
            state = str(row.get("status", ""))
            print(f"  - {ANALYSIS_LABELS.get(target, target)}: {state}")


def _refresh_latest_analysis_shortcuts(project_root: Path, session_root: Path) -> None:
    """
    Maintain stable entry points to the latest analysis session and prune empty dirs.
    """
    try:
        analysis_root = post_docking_root(project_root)
        latest_session_link = analysis_root / "LATEST_SESSION"
        latest_stage_link = analysis_root / "LATEST_STAGE_TARGETS"
        stage_targets = session_root / "stage_targets"

        for link_path, target in (
            (latest_session_link, session_root),
            (latest_stage_link, stage_targets),
        ):
            try:
                if link_path.exists() or link_path.is_symlink():
                    if link_path.is_dir() and not link_path.is_symlink():
                        shutil.rmtree(link_path)
                    else:
                        link_path.unlink()
                if target.exists():
                    link_path.symlink_to(target)
            except Exception:
                continue

        alias_root = project_root / "5-Analysis"
        alias_root.mkdir(parents=True, exist_ok=True)
        alias_readme = alias_root / "README_FIRST.txt"
        lines = [
            "This is the canonical post-docking analysis root.",
            f"Primary root: {analysis_root}",
            f"Latest session: {session_root}",
            f"Latest stage targets: {stage_targets}",
        ]
        alias_readme.write_text("\n".join(lines) + "\n", encoding="utf-8")
        alias_latest = alias_root / "LATEST_SESSION"
        try:
            if alias_latest.exists() or alias_latest.is_symlink():
                alias_latest.unlink()
            alias_latest.symlink_to(session_root)
        except Exception:
            pass

        # Remove only truly empty directories under the analysis root.
        # IMPORTANT: never prune active session roots/sessions container.
        sessions_root = analysis_root / "sessions"
        directories = sorted(
            (path for path in analysis_root.rglob("*") if path.is_dir()),
            key=lambda path: len(path.parts),
            reverse=True,
        )
        for directory in directories:
            try:
                if directory == analysis_root:
                    continue
                if directory == sessions_root or directory == session_root:
                    continue
                if sessions_root in directory.parents:
                    continue
                if any(directory.iterdir()):
                    continue
                directory.rmdir()
            except Exception:
                continue
    except Exception:
        pass


def _run_analysis_stage(questionary, project_root: Path) -> None:
    _header("Post-Docking Analysis")
    manifest = load_manifest(project_root)
    from post_docking_analysis.engine_detector import detect_engines
    targets = _checkbox(
        questionary,
        "Select one or more post-docking functions",
        [(target, label) for target, label in ANALYSIS_LABELS.items()],
    )
    selected_targets, collapsed_alias_targets = _collapse_interaction_alias_targets(list(targets or []))
    if not selected_targets:
        print("No analysis targets selected.")
        return
    clean_target_selected = "analyze.interactions.clean" in selected_targets
    if collapsed_alias_targets:
        collapsed_labels = ", ".join(ANALYSIS_LABELS.get(target, target) for target in collapsed_alias_targets)
        print(
            "ℹ️  Standalone interaction selections now follow the clean-interaction contract.\n"
            f"   Routed to: {ANALYSIS_LABELS.get('analyze.interactions.clean', 'analyze.interactions.clean')}\n"
            f"   Collapsed targets: {collapsed_labels}"
        )

    # Sessionize outputs to avoid scattering artifacts across sibling analyze_* roots.
    analysis_root = _project_layout(project_root)["analysis"]
    sessions_root = analysis_root / "sessions"
    sessions_root.mkdir(parents=True, exist_ok=True)
    existing_sessions = sorted(
        [path for path in sessions_root.glob("interactive_*") if path.is_dir()],
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    session_mode = "new"
    if existing_sessions:
        session_mode = _select(
            questionary,
            "Analysis session mode",
            [
                ("new", "Start a new analysis session"),
                ("resume_latest", f"Resume latest session ({existing_sessions[0].name})"),
            ],
            default="new",
        )
    if session_mode == "resume_latest":
        session_root = existing_sessions[0]
    else:
        session_id = format_run_timestamp()
        session_root = sessions_root / f"interactive_{session_id}"
        session_root.mkdir(parents=True, exist_ok=True)
    shared_stage_output = session_root / "stage_targets"
    shared_stage_output.mkdir(parents=True, exist_ok=True)
    naming_dir = _shared_naming_dir(project_root)
    progress_file = _analysis_progress_file(session_root)
    existing_progress = _read_analysis_progress(progress_file)
    existing_completed_targets = list(existing_progress.get("completed_targets") or [])
    _refresh_latest_analysis_shortcuts(project_root, session_root)

    print(f"Analysis session output: {session_root}")
    if session_mode == "resume_latest" and existing_progress:
        previous_total = int(existing_progress.get("total_targets", 0))
        previous_done = int(existing_progress.get("completed_count", 0))
        print(f"Resuming previous session progress: {previous_done}/{previous_total} completed")

    declared_engines = list(manifest.get("engines", ["gnina"]))
    detection = detect_engines(project_root)
    detected_valid_engines = [
        str(engine).strip().lower()
        for engine in (detection.get("valid_engines") or [])
        if str(engine).strip()
    ]
    if detected_valid_engines:
        engines = detected_valid_engines
    else:
        print("⚠️  No valid docking-engine outputs were auto-detected for this project.")
        print("   Review 4-Working/metadata/engine_detection_report.json for details.")
        allow_manual_fallback = _confirm(
            questionary,
            "Use manual engine fallback selection anyway?",
            default=False,
        )
        if not allow_manual_fallback:
            print("Analysis launch cancelled.")
            return
        engines = declared_engines
    default_engine = str(manifest.get("favorite_engine") or engines[0])
    selected_analysis_engines = list(engines)
    selected_engine_preset = ""
    stage_targets = [t for t in selected_targets if t not in {"analyze.comparative", "analyze.favorite_engine"}]
    stage_engine = None
    run_comparative = "analyze.comparative" in selected_targets

    prompt_protein_names = False
    prompt_ligand_names = False
    complex_query = ""
    exclude_problematic_ligands = False
    positive_affinity_threshold = 0.0
    minimum_pose_count = 1
    rmsd_workers = 0
    consensus_mode = "dockbox_geometric"
    rescoring_scope = "top_n_per_protein"
    rescoring_top_n = 3
    analysis_scope = "full"
    normalization_method = "per_engine_rank"
    biology_file = ""
    biology_mapping_mode = "auto"
    hit_class_policy = "target_percentile"
    hit_class_strong_percentile = 0.10
    hit_class_moderate_percentile = 0.35
    top_pose_policy = str(manifest.get("top_pose_selection_policy", "best_affinity") or "best_affinity").strip().lower()
    if top_pose_policy not in {"best_affinity", "best_consensus", "hybrid"}:
        top_pose_policy = "best_affinity"
    top_pose_aggregation = str(manifest.get("top_pose_global_aggregation", "best_target") or "best_target").strip().lower()
    if top_pose_aggregation not in {"best_target"}:
        top_pose_aggregation = "best_target"
    speed_profile = "standard"
    skip_completed_targets = True

    if run_comparative and stage_targets:
            print("ℹ️  Stage-level engine applies to stage/interactions/visual targets.")
            print("ℹ️  Favorite-engine continuation reuses stage-level engine when available.")
    if len(engines) == 1:
            print(f"ℹ️  Auto-detected engine context: {engines[0].upper()} only; engine prompts will auto-fill.")

    configure_required = True
    while True:
        if configure_required:
            prompt_protein_names = _confirm(
                questionary,
                "Prompt protein display names before running selected targets?",
                default=False,
            )
            prompt_ligand_names = _confirm(
                questionary,
                "Prompt ligand display names before running selected targets?",
                default=False,
            )
            complex_query = _optional_text(
                questionary,
                "Optional complex filter query (example: protein=2FVD;ligand=Sorafenib or exclude_ligand=Acetazolamide; Enter/0 keeps all)",
                default=complex_query or "",
            )
            if complex_query == "0":
                complex_query = ""
            analysis_scope = _select(
                questionary,
                "Analysis scope",
                [
                    ("full", "Full analysis (recommended)"),
                    ("comparison_only", "Comparison only"),
                    ("rescoring_only", "Rescoring shortlist only"),
                    ("top_pose_only", "Top-pose atlas only"),
                    ("qc_only", "QC-only"),
                    ("report_only", "Reports only (reuse score artifacts)"),
                ],
                default=analysis_scope,
            )
            normalization_method = _select(
                questionary,
                "Cross-engine score normalization method",
                [
                    ("per_engine_rank", "Per-engine rank percentile (recommended)"),
                    ("per_engine_minmax", "Per-engine min-max"),
                    ("per_engine_zscore", "Per-engine z-score"),
                ],
                default=normalization_method,
            )
            biology_file = _optional_text(
                questionary,
                "Optional biology annotation file (CSV/TSV/JSON; Enter to skip)",
                default=biology_file or "",
            )
            if biology_file == "0":
                biology_file = ""
            biology_mapping_mode = _select(
                questionary,
                "Biology mapping mode",
                [
                    ("auto", "Auto detect"),
                    ("tag", "Map by tag"),
                    ("protein_ligand", "Map by protein + ligand"),
                    ("protein", "Map by protein"),
                    ("ligand", "Map by ligand"),
                ],
                default=biology_mapping_mode,
            )
            hit_class_policy = _select(
                questionary,
                "Hit classification policy",
                [
                    ("target_percentile", "Target percentile (default)"),
                    ("reference_anchor", "Reference-anchor (validated redocking; use when references exist)"),
                ],
                default=hit_class_policy,
            )
            hit_class_strong_percentile = _float_input(
                questionary,
                "Strong-hit percentile threshold (per target)",
                default=hit_class_strong_percentile,
            )
            hit_class_moderate_percentile = _float_input(
                questionary,
                "Moderate-hit percentile threshold (per target)",
                default=hit_class_moderate_percentile,
            )
            top_pose_policy = _select(
                questionary,
                "Top-pose atlas selection policy",
                [
                    ("best_affinity", "Best affinity first (recommended)"),
                    ("hybrid", "Hybrid (consensus then affinity)"),
                    ("best_consensus", "Best consensus first"),
                ],
                default=top_pose_policy,
            )
            top_pose_aggregation = _select(
                questionary,
                "Top-pose global aggregation mode",
                [("best_target", "Best target (recommended)")],
                default=top_pose_aggregation,
            )

            if len(engines) > 1:
                selected_analysis_engines, selected_engine_preset = _analysis_engine_scope_picker(
                    questionary,
                    manifest,
                    detection,
                    preselected=selected_analysis_engines,
                    preset_name=selected_engine_preset,
                )
            else:
                selected_analysis_engines = list(engines)
                selected_engine_preset = ""
            scoped_default_engine = (
                default_engine if default_engine in selected_analysis_engines else selected_analysis_engines[0]
            )

            if stage_targets and len(selected_analysis_engines) > 1:
                stage_engine = _select(
                    questionary,
                    "Engine for selected stage-level targets (not the same as Favorite-engine continuation)",
                    [(engine_name, engine_name.upper()) for engine_name in selected_analysis_engines],
                    default=stage_engine or scoped_default_engine,
                )
            elif stage_targets and len(selected_analysis_engines) == 1:
                stage_engine = selected_analysis_engines[0]

            if run_comparative:
                consensus_mode = _select(
                    questionary,
                    "Consensus mode for comparative hit ranking",
                    [
                        ("dockbox_geometric", "DockBox geometric consensus (recommended)"),
                        ("weighted_hybrid", "Weighted hybrid"),
                        ("strict_consensus", "Strict consensus (requires >=2 engines per pair)"),
                        ("favorite_guardrails", "Favorite + guardrails"),
                    ],
                    default=consensus_mode,
                )
                rescoring_scope = _select(
                    questionary,
                    "Rescoring scope for shortlist generation",
                    [
                        ("top_n_per_protein", "Top-N per protein (recommended)"),
                        ("top_n_global", "Top-N global"),
                    ],
                    default=rescoring_scope,
                )
                rescoring_top_n = _int_input(
                    questionary,
                    "Rescoring Top-N candidates",
                    default=rescoring_top_n,
                    minimum=1,
                )

            exclude_problematic_ligands = _confirm(
                questionary,
                "Exclude problematic ligands before analysis (positive affinity / missing poses)?",
                default=exclude_problematic_ligands,
            )
            if exclude_problematic_ligands:
                positive_affinity_threshold = _float_input(
                    questionary,
                    "Exclude ligands when best affinity is greater than this threshold (kcal/mol)",
                    default=positive_affinity_threshold,
                )
                minimum_pose_count = _int_input(
                    questionary,
                    "Minimum required poses per ligand",
                    default=minimum_pose_count,
                    minimum=1,
                )
            rmsd_workers = _rmsd_workers_input(questionary, default=rmsd_workers)
            speed_profile = _select(
                questionary,
                "Analysis speed profile",
                [
                    ("standard", "Standard (run all selected stages)"),
                    ("fast", "Fast mode (skip optional heavy stages where supported)"),
                ],
                default=speed_profile,
            )
            skip_completed_targets = _confirm(
                questionary,
                "Skip targets already completed in this session?",
                default=skip_completed_targets,
            )
            configure_required = False

        action = _select(
            questionary,
            "Analysis control",
            [
                ("run", "Run selected targets now"),
                ("progress", "Show current session progress"),
                ("edit", "Edit analysis settings"),
                ("cancel", "Back to workflow menu"),
                ("exit", "Exit wizard"),
            ],
            default="run",
        )
        if action == "progress":
            _show_analysis_progress(progress_file)
            continue
        if action == "edit":
            configure_required = True
            continue
        if action == "exit":
            raise _PromptNavigation("exit")
        if action == "cancel":
            print("Analysis launch cancelled.")
            return
        break

    remaining_prompt_protein = prompt_protein_names
    remaining_prompt_ligand = prompt_ligand_names
    completed_targets: List[Dict[str, str]] = list(existing_completed_targets)
    _write_analysis_progress(
        progress_file,
        session_root,
        selected_targets,
        completed_targets,
        running_target="",
        status="started",
        note="Session initialized" if session_mode == "new" else "Session resumed",
    )

    completed_status_map = _completed_target_status_map(completed_targets)
    rmsd_scopes_setting = "per_complex"
    if speed_profile == "fast":
        print("⚡ Fast mode enabled: optional heavy stages are minimized where supported.")
    total_targets = len(selected_targets)
    for index, target in enumerate(selected_targets, start=1):
        label = ANALYSIS_LABELS.get(target, target)
        if skip_completed_targets and completed_status_map.get(target) == "completed":
            print(f"[{index}/{total_targets}] Skipping {label} (already completed in this session)")
            _write_analysis_progress(
                progress_file,
                session_root,
                selected_targets,
                completed_targets,
                running_target="",
                status="running",
                note=f"Skipped {label}: already completed",
            )
            continue
        print(f"[{index}/{total_targets}] Starting {label} ...")
        _write_analysis_progress(
            progress_file,
            session_root,
            selected_targets,
            completed_targets,
            running_target=target,
            status="running",
            note=f"Running {label}",
        )
        if target == "analyze.comparative":
            comparative_output = session_root / "comparative_all_engines"
            promote_exhaustive = _confirm(
                questionary,
                "Generate a follow-up rerun manifest from comparative results? (creates list only, does not run docking)",
                default=False,
            )
            rerun_engine = None
            top_per_protein = 1
            winner_only = False
            min_affinity_advantage = 0.0
            max_rerun_pairs = 0
            pair_allowlist = None
            if promote_exhaustive:
                rerun_engine = _select(
                    questionary,
                    "Engine to promote into exhaustive follow-up",
                    [(engine_name, engine_name.upper()) for engine_name in selected_analysis_engines],
                    default=scoped_default_engine,
                )
                top_per_protein = int(_text(questionary, "Top ligands per protein to consider", default="1", required=True))
                winner_only = _confirm(questionary, "Only promote pairs won by the selected engine?", default=False)
                min_affinity_advantage = float(
                    _text(questionary, "Minimum affinity advantage required (kcal/mol)", default="0.0", required=True)
                )
                max_rerun_pairs = int(_text(questionary, "Maximum promoted pairs (0 = no limit)", default="0", required=True))
                pair_allowlist = _optional_text(questionary, "Optional pair allowlist file", default="")
            result = run_analysis_comparative(
                str(project_root),
                output_dir=str(comparative_output),
                engines=selected_analysis_engines,
                engine_preset=selected_engine_preset or None,
                promote_exhaustive=promote_exhaustive,
                rerun_engine=rerun_engine,
                top_per_protein=top_per_protein,
                winner_only=winner_only,
                min_affinity_advantage=min_affinity_advantage,
                max_rerun_pairs=max_rerun_pairs,
                pair_allowlist=pair_allowlist or None,
                consensus_mode=consensus_mode,
                rescoring_scope=rescoring_scope,
                rescoring_top_n=rescoring_top_n,
                prompt_protein_names=remaining_prompt_protein,
                prompt_ligand_names=remaining_prompt_ligand,
                complex_query=complex_query or None,
                shared_mapping_dir=str(naming_dir),
                exclude_problematic_ligands=exclude_problematic_ligands,
                positive_affinity_threshold=positive_affinity_threshold,
                minimum_pose_count=minimum_pose_count,
                rmsd_workers=rmsd_workers,
                analysis_scope=analysis_scope,
                normalization_method=normalization_method,
                biology_file=biology_file or None,
                biology_mapping_mode=biology_mapping_mode,
                hit_class_policy=hit_class_policy,
                hit_class_strong_percentile=hit_class_strong_percentile,
                hit_class_moderate_percentile=hit_class_moderate_percentile,
                top_pose_selection_policy=top_pose_policy,
                top_pose_global_aggregation=top_pose_aggregation,
            )
            if rerun_manifest := result.outputs.get("rerun_manifest_file"):
                print(f"Follow-up rerun manifest generated: {rerun_manifest}")
                print("Note: Docking was not re-run automatically. Use this manifest in your next docking round.")
            if consensus_file := result.outputs.get("consensus_ranked_hits_file"):
                print(f"Consensus ranked hits: {consensus_file}")
            if rescoring_file := result.outputs.get("rescoring_candidates_file"):
                print(f"Rescoring candidates: {rescoring_file}")
            if explain_file := result.outputs.get("consensus_explainability_file"):
                print(f"Consensus explainability: {explain_file}")
            if corr_file := result.outputs.get("engine_rank_correlation_per_protein_file"):
                print(f"Per-protein rank correlations: {corr_file}")
            if corr_global_file := result.outputs.get("engine_rank_correlation_global_file"):
                print(f"Global rank correlations: {corr_global_file}")
            if canonical_scores := result.outputs.get("canonical_scores_root"):
                print(f"Canonical score outputs: {canonical_scores}")
            if classes_file := result.outputs.get("consensus_ranked_hits_with_classes_file"):
                print(f"Hit classes file: {classes_file}")
            if bio_map := result.outputs.get("biology_mapping_report_file"):
                print(f"Biology mapping report: {bio_map}")
            if bio_corr := result.outputs.get("biology_correlation_global_file"):
                print(f"Biology global correlations: {bio_corr}")
            if top_pose_global := result.outputs.get("top_pose_per_ligand_global_file"):
                print(f"Top-pose global atlas: {top_pose_global}")
            if top_pose_summary := result.outputs.get("top_pose_summary_file"):
                print(f"Top-pose summary: {top_pose_summary}")
            if top_pose_canonical := result.outputs.get("top_pose_canonical_root"):
                print(f"Top-pose canonical mirror: {top_pose_canonical}")
        elif target == "analyze.favorite_engine":
            favorite_output = session_root / "favorite_engine"
            if stage_engine and stage_engine in selected_analysis_engines:
                favorite_engine = stage_engine
                print(f"Using stage-level engine for favorite continuation: {favorite_engine.upper()}")
            elif len(selected_analysis_engines) == 1:
                favorite_engine = selected_analysis_engines[0]
                print(f"Using only available engine for favorite continuation: {favorite_engine.upper()}")
            else:
                favorite_engine = _select(
                    questionary,
                    "Favorite engine for this continuation target",
                    [(engine, engine.upper()) for engine in selected_analysis_engines],
                    default=scoped_default_engine,
                )
            result = run_analysis_favorite(
                str(project_root),
                favorite_engine,
                output_dir=str(favorite_output),
                engines=selected_analysis_engines,
                engine_preset=selected_engine_preset or None,
                prompt_protein_names=remaining_prompt_protein,
                prompt_ligand_names=remaining_prompt_ligand,
                complex_query=complex_query or None,
                rmsd_scopes=rmsd_scopes_setting,
                shared_mapping_dir=str(naming_dir),
                exclude_problematic_ligands=exclude_problematic_ligands,
                positive_affinity_threshold=positive_affinity_threshold,
                minimum_pose_count=minimum_pose_count,
                rmsd_workers=rmsd_workers,
                speed_profile=speed_profile,
                analysis_scope=analysis_scope,
                normalization_method=normalization_method,
                biology_file=biology_file or None,
                biology_mapping_mode=biology_mapping_mode,
                hit_class_policy=hit_class_policy,
                hit_class_strong_percentile=hit_class_strong_percentile,
                hit_class_moderate_percentile=hit_class_moderate_percentile,
                top_pose_selection_policy=top_pose_policy,
                top_pose_global_aggregation=top_pose_aggregation,
            )
            if top_pose_global := result.outputs.get("top_pose_per_ligand_global_file"):
                print(f"Top-pose global atlas: {top_pose_global}")
            if top_pose_summary := result.outputs.get("top_pose_summary_file"):
                print(f"Top-pose summary: {top_pose_summary}")
            if top_pose_canonical := result.outputs.get("top_pose_canonical_root"):
                print(f"Top-pose canonical mirror: {top_pose_canonical}")
            if clean_target_selected:
                print(
                    "ℹ️  Clean interactions pipeline is already selected in this session target list; "
                    "skipping duplicate favorite-flow trigger."
                )
            elif str(result.status) == "completed":
                clean_followup_mode = _select(
                    questionary,
                    "Favorite-engine continuation interaction follow-up",
                    [
                        ("run_clean", "Run clean interactions pipeline now (recommended)"),
                        ("skip_clean", "Skip clean interactions follow-up"),
                    ],
                    default="run_clean",
                )
                if clean_followup_mode == "run_clean":
                    print("Running clean interactions follow-up from favorite-engine continuation ...")
                    clean_followup_result = run_analysis_target(
                        "analyze.interactions.clean",
                        project_dir=str(project_root),
                        output_dir=str(favorite_output / "clean_interactions"),
                        engine=favorite_engine,
                        engines=selected_analysis_engines,
                        engine_preset=selected_engine_preset or None,
                        prompt_protein_names=False,
                        prompt_ligand_names=False,
                        complex_query=complex_query or None,
                        rmsd_scopes=rmsd_scopes_setting,
                        shared_mapping_dir=str(naming_dir),
                        exclude_problematic_ligands=exclude_problematic_ligands,
                        positive_affinity_threshold=positive_affinity_threshold,
                        minimum_pose_count=minimum_pose_count,
                        rmsd_workers=rmsd_workers,
                        speed_profile=speed_profile,
                    )
                    result.outputs["favorite_clean_followup_target"] = str(clean_followup_result.target)
                    result.outputs["favorite_clean_followup_status"] = str(clean_followup_result.status)
                    result.outputs["favorite_clean_followup_output_dir"] = str(
                        clean_followup_result.outputs.get("output_dir", "")
                    )
                    print(f"Favorite clean-interactions follow-up: {clean_followup_result.status}")
                    print("Follow-up run-tracking summary:")
                    _print_run_tracking_summary(clean_followup_result)
                    if str(clean_followup_result.status) != "completed":
                        result.status = "failed"
                        result.notes.append(
                            f"favorite_clean_followup_failed={clean_followup_result.status}"
                        )
            else:
                print("Skipping favorite-flow clean interactions follow-up because favorite-engine run did not complete.")
        else:
            result = run_analysis_target(
                target,
                project_dir=str(project_root),
                output_dir=str(shared_stage_output),
                engine=stage_engine or default_engine,
                engines=selected_analysis_engines,
                engine_preset=selected_engine_preset or None,
                enable_poseview=(target == "analyze.interactions.poseview"),
                prompt_protein_names=remaining_prompt_protein,
                prompt_ligand_names=remaining_prompt_ligand,
                complex_query=complex_query or None,
                rmsd_scopes=rmsd_scopes_setting,
                shared_mapping_dir=str(naming_dir),
                exclude_problematic_ligands=exclude_problematic_ligands,
                positive_affinity_threshold=positive_affinity_threshold,
                minimum_pose_count=minimum_pose_count,
                rmsd_workers=rmsd_workers,
                speed_profile=speed_profile,
            )
        if remaining_prompt_protein or remaining_prompt_ligand:
            remaining_prompt_protein = False
            remaining_prompt_ligand = False
        completed_targets.append(
            {
                "target": target,
                "status": str(result.status),
                "timestamp": datetime.now().isoformat(timespec="seconds"),
            }
        )
        completed_status_map[target] = str(result.status)
        _write_analysis_progress(
            progress_file,
            session_root,
            selected_targets,
            completed_targets,
            running_target="",
            status="running",
            note=f"Completed {label}: {result.status}",
        )
        print(f"[{index}/{total_targets}] {label}: {result.status}")
        _print_run_tracking_summary(result)

    _write_analysis_progress(
        progress_file,
        session_root,
        selected_targets,
        completed_targets,
        running_target="",
        status="completed",
        note="All selected targets finished",
    )
    print(f"Session progress file: {progress_file}")
    _refresh_latest_analysis_shortcuts(project_root, session_root)


def run_interactive_workflow(initial_root: Optional[str] = None) -> int:
    questionary = _questionary()
    print("Omni-DockForge: End-to-End Docking, Consensus by Design.")
    print("================================")

    startup_root = _suggest_startup_project_root(
        Path(initial_root).expanduser().resolve() if initial_root else None
    )
    active_root = Path(initial_root).expanduser().resolve() if initial_root else None
    if active_root is not None and not active_root.exists():
        active_root = None
    if active_root is not None and not (active_root / "4-Docking" / "project_manifest.json").exists() and not (active_root / "project_manifest.json").exists():
        run_workflow_init(str(active_root), layout_profile="docking_legacy")

    while active_root is None:
        try:
            opening_choice = _select(
                questionary,
                "Docking project workflow",
                [
                    ("new", "Create a new docking project"),
                    ("resume", "Resume existing docking project"),
                    ("status", "Inspect existing project status"),
                    ("quit", "Exit"),
                ],
                default="new",
            )
        except _PromptNavigation as nav:
            if nav.action == "exit":
                return 0
            continue
        if opening_choice == "quit":
            return 0
        if opening_choice == "new":
            try:
                active_root = _create_project(questionary)
            except _PromptNavigation as nav:
                if nav.action == "exit":
                    return 0
                continue
            break
        if opening_choice == "resume":
            try:
                active_root = _resume_project(questionary, initial_root=active_root or startup_root)
            except _PromptNavigation as nav:
                if nav.action == "exit":
                    return 0
                continue
            break
        if opening_choice == "status":
            try:
                _show_status(_resume_project(questionary, initial_root=active_root or startup_root))
            except _PromptNavigation as nav:
                if nav.action == "exit":
                    return 0
                continue

    stage_default = "analyze"
    while True:
        _show_workflow_timeline(active_root)
        try:
            stage = _select(
                questionary,
                "Workflow stage",
                [
                    ("pdb_collect", "Fetch proteins from PDB"),
                    ("prepare_assets", "Prepare raw proteins and ligands"),
                    ("build_pairlist", "Build pairlist"),
                    ("materialize_project", "Prepare docking folders"),
                    ("checkpoint", "Checkpoint & Revise workspace clone"),
                    ("dock", "Run docking"),
                    ("hpc", "Deploy, sync, and submit on HPC"),
                    ("tasks", "Show background task progress"),
                    ("edit_names", "Edit protein/ligand display names"),
                    ("analyze", "Post-docking analysis"),
                    ("status", "Show workflow status"),
                    ("switch", "Switch project"),
                    ("quit", "Exit"),
                ],
                default=stage_default,
            )
        except _PromptNavigation as nav:
            if nav.action == "exit":
                return 0
            continue
        stage_default = stage

        if stage == "quit":
            return 0
        try:
            if stage == "switch":
                active_root = _resume_project(questionary, initial_root=active_root)
                continue
            if stage == "status":
                _show_status(active_root)
                continue
            if stage == "pdb_collect":
                _collect_pdbs(questionary, active_root)
                continue
            if stage == "prepare_assets":
                _prepare_assets(questionary, active_root)
                continue
            if stage == "build_pairlist":
                _build_pairlist(questionary, active_root)
                continue
            if stage == "materialize_project":
                _materialize_project(questionary, active_root)
                continue
            if stage == "checkpoint":
                _create_checkpoint_workspace(questionary, active_root)
                continue
            if stage == "dock":
                _run_docking_stage(questionary, active_root)
                continue
            if stage == "hpc":
                _run_hpc_stage(questionary, active_root)
                continue
            if stage == "tasks":
                _show_background_tasks(active_root)
                continue
            if stage == "edit_names":
                _edit_display_names(questionary, active_root)
                continue
            if stage == "analyze":
                _run_analysis_stage(questionary, active_root)
                continue
        except _PromptNavigation as nav:
            if nav.action == "exit":
                return 0
            print("↩️ Returning to workflow menu.")
            continue
