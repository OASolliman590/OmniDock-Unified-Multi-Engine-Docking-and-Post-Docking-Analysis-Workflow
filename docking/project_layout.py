from __future__ import annotations

import json
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from .models import ProjectManifest


PAIRLIST_COLUMNS = [
    "receptor",
    "site_id",
    "ligand",
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
    "protein_display_name",
    "ligand_display_name",
]

LAYOUT_CANONICAL = "canonical"
LAYOUT_DOCKING_LEGACY = "docking_legacy"
LAYOUT_PROFILES = {LAYOUT_CANONICAL, LAYOUT_DOCKING_LEGACY}

DOCKING_LEGACY_DIRS = {
    "raw_ligands": "1-Raw_Ligand",
    # Compatibility alias: normalized SDF copies now live beside the source
    # ligands in the single raw-ligands staging folder.
    "raw_ligands_sdf": "1-Raw_Ligand",
    "raw_proteins": "2-Raw_Protien",
    "preparation": "3-Preparation",
    "prepared_proteins": "3-Preparation/2-Prepared_Protiens",
    "prepared_ligands": "3-Preparation/4-Prepared_Ligand",
    "preparation_logs": "3-Preparation/preparation_logs",
    "docking_root": "4-Docking",
    "receptors": "4-Docking/receptors",
    "ligands": "4-Docking/ligands",
    "metadata": "4-Docking/metadata",
    "post_docking_root": "5-Analysis",
    "analysis": "5-Analysis",
    "raw_data": "5-Analysis/raw_data",
    "reports": "7-Reports",
    "visualizations": "5-Analysis/visualizations",
    "interactions": "5-Analysis/interactions",
    "rmsd_analysis": "5-Analysis/rmsd_analysis",
    "complexes": "5-Analysis/complexes",
    "best_poses": "5-Analysis/best_poses",
}

GNINA_HPC_COMPAT_DIRS = [
    "results",
    "scripts",
    "scripts_hpc",
]

NUMBERED_OUTPUT_TOPOLOGY = {
    "meta": ".meta",
    "input_root": "0-Input",
    "input_receptors": "0-Input/receptors",
    "input_ligands": "0-Input/ligands",
    "input_references": "0-Input/references",
    "preparation_root": "1-Preparation",
    "prep_receptors_prepared": "1-Preparation/receptors/prepared",
    "prep_receptors_qc": "1-Preparation/receptors/qc",
    "prep_ligands_prepared": "1-Preparation/ligands/prepared",
    "prep_ligands_flagged": "1-Preparation/ligands/flagged",
    "gridboxes_root": "2-GridBoxes",
    "grid_definitions": "2-GridBoxes/definitions",
    "grid_validation": "2-GridBoxes/validation",
    "docking_root_numbered": "3-Docking",
    "docking_gnina_configs": "3-Docking/gnina/configs",
    "docking_gnina_poses": "3-Docking/gnina/poses",
    "docking_gnina_logs": "3-Docking/gnina/logs",
    "docking_vina_configs": "3-Docking/vina/configs",
    "docking_vina_poses": "3-Docking/vina/poses",
    "docking_vina_logs": "3-Docking/vina/logs",
    "docking_smina_configs": "3-Docking/smina/configs",
    "docking_smina_poses": "3-Docking/smina/poses",
    "docking_smina_logs": "3-Docking/smina/logs",
    "docking_autodock4_configs": "3-Docking/autodock4/configs",
    "docking_autodock4_poses": "3-Docking/autodock4/poses",
    "docking_autodock4_logs": "3-Docking/autodock4/logs",
    "post_docking_root_numbered": "4-Working",
    "post_metadata": "4-Working/metadata",
    "post_scores_raw": "4-Working/scores/raw",
    "post_scores_unified": "4-Working/scores/unified",
    "post_scores_consensus": "4-Working/scores/consensus",
    "post_filtering": "4-Working/filtering",
    "post_filter_reports": "4-Working/filtering/filter_reports",
    "post_rescoring_candidates": "4-Working/rescoring",
    "post_rescoring_results": "4-Working/rescoring/rescored",
    "post_rerun_manifests": "4-Working/rerun_manifests",
    "post_storage": "4-Working/storage",
    "analysis_root_numbered": "5-Analysis",
    "analysis_sessions_root": "5-Analysis/sessions",
    "analysis_comparative": "5-Analysis/comparative",
    "analysis_raw_data": "5-Analysis/raw_data",
    "analysis_rmsd": "5-Analysis/rmsd",
    "analysis_polypharmacology": "5-Analysis/polypharmacology",
    "analysis_structure_quality": "5-Analysis/structure_quality",
    "analysis_top_pose": "5-Analysis/top_pose_ligand_performance",
    "visualizations_root_numbered": "6-Visualizations",
    "visualizations_2d": "6-Visualizations/2d",
    "visualizations_2d_heatmaps": "6-Visualizations/2d/heatmaps",
    "visualizations_2d_score_plots": "6-Visualizations/2d/score_plots",
    "visualizations_2d_interaction_maps": "6-Visualizations/2d/interaction_maps",
    "visualizations_3d": "6-Visualizations/3d",
    "visualizations_3d_py3dmol": "6-Visualizations/3d/py3dmol",
    "visualizations_3d_pymol": "6-Visualizations/3d/pymol",
    "reports_root_numbered": "7-Reports",
    "reports_per_protein": "7-Reports/per_protein",
    "reports_figures": "7-Reports/figures",
    "reports_tables": "7-Reports/tables",
    "sessions_root_numbered": "sessions",
}


def detect_layout_profile(project_root: Path, requested: Optional[str] = None) -> str:
    if requested:
        if requested not in LAYOUT_PROFILES:
            raise ValueError(f"Unsupported layout profile: {requested}")
        return requested

    root = Path(project_root).expanduser().resolve()
    if (root / DOCKING_LEGACY_DIRS["docking_root"]).exists():
        return LAYOUT_DOCKING_LEGACY
    return LAYOUT_CANONICAL


def docking_root(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    if profile == LAYOUT_DOCKING_LEGACY:
        return root / DOCKING_LEGACY_DIRS["docking_root"]
    return root


def post_docking_root(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    if profile == LAYOUT_DOCKING_LEGACY:
        return root / DOCKING_LEGACY_DIRS["post_docking_root"]
    return root / "5-Analysis"


def _canonical_layout(project_root: Path) -> Dict[str, Path]:
    return {
        "project_root": project_root,
        "raw_ligands": project_root / "ligands_raw",
        "raw_ligands_sdf": project_root / "ligands_raw",
        "raw_proteins": project_root / "proteins_raw",
        "preparation": project_root / "preparation",
        "prepared_proteins": project_root / "prepared_proteins",
        "prepared_ligands": project_root / "prepared_ligands",
        "preparation_logs": project_root / "preparation_logs",
        "docking_root": project_root,
        "receptors": project_root / "receptors",
        "ligands": project_root / "ligands",
        "metadata": project_root / "metadata",
        "post_docking_root": project_root / "5-Analysis",
        "analysis": project_root / "5-Analysis",
        "raw_data": project_root / "5-Analysis" / "raw_data",
        "reports": project_root / "7-Reports",
        "visualizations": project_root / "5-Analysis" / "visualizations",
        "interactions": project_root / "5-Analysis" / "interactions",
        "rmsd_analysis": project_root / "5-Analysis" / "rmsd_analysis",
        "complexes": project_root / "5-Analysis" / "complexes",
        "best_poses": project_root / "5-Analysis" / "best_poses",
    }


def _legacy_layout(project_root: Path) -> Dict[str, Path]:
    layout = {"project_root": project_root}
    for key, rel_path in DOCKING_LEGACY_DIRS.items():
        layout[key] = project_root / rel_path
    return layout


def ensure_project_layout(project_root: Path, layout_profile: str = LAYOUT_CANONICAL) -> Dict[str, Path]:
    project_root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(project_root, layout_profile)
    layout = _legacy_layout(project_root) if profile == LAYOUT_DOCKING_LEGACY else _canonical_layout(project_root)
    layout["layout_profile"] = profile
    for key, path in layout.items():
        if key in {"project_root", "layout_profile"}:
            continue
        path.mkdir(parents=True, exist_ok=True)
    ensure_post_docking_compat_shim(project_root)
    return layout


def ensure_engine_layout(project_root: Path, engine: str, layout_profile: Optional[str] = None) -> Dict[str, Path]:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    dock_root = docking_root(root, profile)

    if profile == LAYOUT_DOCKING_LEGACY:
        if engine == "gnina":
            layout = {
                "root": dock_root,
                "poses": dock_root / "gnina_out",
                "logs": dock_root / "logs",
                "scores": dock_root / "results",
            }
            for path in layout.values():
                path.mkdir(parents=True, exist_ok=True)
            for dirname in GNINA_HPC_COMPAT_DIRS:
                (dock_root / dirname).mkdir(parents=True, exist_ok=True)
            return layout

        engine_root = dock_root / f"{engine}_out"
        layout = {
            "root": engine_root,
            "poses": engine_root / "poses",
            "logs": engine_root / "logs",
            "scores": engine_root / "scores",
        }
        for path in layout.values():
            path.mkdir(parents=True, exist_ok=True)
        return layout

    engine_root = dock_root / "engines" / engine
    layout = {
        "root": engine_root,
        "poses": engine_root / "poses",
        "logs": engine_root / "logs",
        "scores": engine_root / "scores",
    }
    for path in layout.values():
        path.mkdir(parents=True, exist_ok=True)
    return layout


def ensure_gnina_hpc_compat_layout(project_root: Path, layout_profile: Optional[str] = None) -> Dict[str, str]:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    gnina_layout = ensure_engine_layout(root, "gnina", layout_profile=profile)
    dock_root = docking_root(root, profile)

    created_dirs: Dict[str, str] = {
        "gnina_out": str(gnina_layout["poses"]),
        "logs": str(gnina_layout["logs"]),
    }
    for dirname in GNINA_HPC_COMPAT_DIRS:
        path = dock_root / dirname
        path.mkdir(parents=True, exist_ok=True)
        created_dirs[dirname] = str(path)
    return created_dirs


def ensure_numbered_output_layout(project_root: Path, layout_profile: Optional[str] = None) -> Dict[str, Path]:
    """
    Ensure canonical numbered DockForge output folders exist.

    Numbered DockForge folders are always created at project root. This keeps
    human-facing analysis outputs in `5-Analysis/` and intermediate scratch data
    in `4-Working/`, regardless of whether the project started from a legacy dock
    layout.
    """
    root = Path(project_root).expanduser().resolve()
    detect_layout_profile(root, layout_profile)
    numbered_base = root

    layout: Dict[str, Path] = {}
    for key, relative in NUMBERED_OUTPUT_TOPOLOGY.items():
        path = numbered_base / relative
        path.mkdir(parents=True, exist_ok=True)
        layout[key] = path
    layout["numbered_base"] = numbered_base
    ensure_post_docking_compat_shim(root)
    return layout


def ensure_post_docking_compat_shim(project_root: Path) -> Dict[str, object]:
    root = Path(project_root).expanduser().resolve()
    canonical = root / "5-Analysis"
    legacy = root / "5-Post_Docking_Analysis"
    canonical.mkdir(parents=True, exist_ok=True)

    status = "created"
    if legacy.exists() or legacy.is_symlink():
        if legacy.is_symlink():
            status = "existing_symlink"
        else:
            status = "legacy_dir_exists"
            warnings.warn(
                f"Legacy post-docking directory exists at {legacy}. New runs write to {canonical}.",
                RuntimeWarning,
                stacklevel=2,
            )
    else:
        legacy.symlink_to(canonical, target_is_directory=True)
        warnings.warn(
            f"Created compatibility symlink {legacy} -> {canonical}. Use {canonical} for new analysis outputs.",
            RuntimeWarning,
            stacklevel=2,
        )
    return {
        "status": status,
        "canonical_root": canonical,
        "legacy_root": legacy,
    }


def pairlist_path(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return docking_root(project_root, layout_profile) / "pairlist.csv"


def pair_intent_path(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return docking_root(project_root, layout_profile) / "metadata" / "pair_intent.csv"


def pair_curation_state_path(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return docking_root(project_root, layout_profile) / "metadata" / "pair_curation_state.json"


def manifest_path(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return docking_root(project_root, layout_profile) / "project_manifest.json"


def deployment_root(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return docking_root(project_root, layout_profile) / "deployments"


def shared_receptors_dir(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return ensure_project_layout(project_root, detect_layout_profile(project_root, layout_profile))["receptors"]


def shared_ligands_dir(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    return ensure_project_layout(project_root, detect_layout_profile(project_root, layout_profile))["ligands"]


def load_pairlist(project_root: Path, layout_profile: Optional[str] = None) -> pd.DataFrame:
    return pd.read_csv(pairlist_path(project_root, layout_profile))


def load_manifest(project_root: Path, layout_profile: Optional[str] = None) -> Dict[str, object]:
    root = Path(project_root).expanduser().resolve()
    primary = manifest_path(root, layout_profile)
    candidates = [primary]
    root_candidate = root / "project_manifest.json"
    legacy_candidate = root / DOCKING_LEGACY_DIRS["docking_root"] / "project_manifest.json"

    for candidate in (root_candidate, legacy_candidate):
        if candidate not in candidates:
            candidates.append(candidate)

    for candidate in candidates:
        if not candidate.exists():
            continue
        with open(candidate, "r", encoding="utf-8") as handle:
            return json.load(handle)

    raise FileNotFoundError(f"No project manifest found. Checked: {', '.join(str(path) for path in candidates)}")


def save_manifest(project_root: Path, payload: Dict[str, object], layout_profile: Optional[str] = None) -> Path:
    path = manifest_path(project_root, layout_profile)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    return path


def ensure_pairlist_stub(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    path = pairlist_path(project_root, layout_profile)
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(columns=PAIRLIST_COLUMNS).to_csv(path, index=False)
    return path


def bootstrap_project_layout(
    project_root: Path,
    engines: List[str],
    project_name: str = "",
    favorite_engine: str = "",
    asset_mode: str = "symlink",
    source_paths: Optional[Dict[str, str]] = None,
    notes: Optional[List[str]] = None,
    layout_profile: str = LAYOUT_CANONICAL,
    pairlist_mode: str = "",
    has_cocrystal_benchmark_rows: bool = False,
) -> Dict[str, object]:
    project_root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(project_root, layout_profile)
    layout = ensure_project_layout(project_root, profile)
    pairlist_file = ensure_pairlist_stub(project_root, profile)
    existing = load_manifest(project_root, profile) if manifest_path(project_root, profile).exists() else {}

    for engine in engines:
        ensure_engine_layout(project_root, engine, layout_profile=profile)

    compatibility_profiles: List[str] = list(existing.get("compatibility_profiles", []))
    if "gnina" in engines:
        ensure_gnina_hpc_compat_layout(project_root, layout_profile=profile)
        if "gnina_hpc" not in compatibility_profiles:
            compatibility_profiles.append("gnina_hpc")

    merged_engines = sorted(set(existing.get("engines", [])) | set(engines))
    merged_source_paths = dict(existing.get("source_paths", {}))
    merged_source_paths.update(source_paths or {})
    merged_notes = list(existing.get("notes", []))
    for note in notes or []:
        if note not in merged_notes:
            merged_notes.append(note)
    if "gnina" in merged_engines:
        compat_note = "compatibility_profile=gnina_hpc"
        if compat_note not in merged_notes:
            merged_notes.append(compat_note)

    manifest = ProjectManifest(
        project_name=project_name or existing.get("project_name") or project_root.name,
        project_root=str(project_root),
        created_at=existing.get("created_at") or datetime.now(timezone.utc).isoformat(),
        asset_mode=existing.get("asset_mode") or asset_mode,
        engines=merged_engines,
        pairlist_file=str(pairlist_file),
        receptors_dir=str(layout["receptors"]),
        ligands_dir=str(layout["ligands"]),
        metadata_dir=str(layout["metadata"]),
        layout_profile=profile,
        docking_root=str(layout["docking_root"]),
        post_docking_root=str(layout["post_docking_root"]),
        raw_proteins_dir=str(layout["raw_proteins"]),
        raw_ligands_dir=str(layout["raw_ligands"]),
        raw_ligands_sdf_dir=str(layout["raw_ligands_sdf"]),
        prepared_proteins_dir=str(layout["prepared_proteins"]),
        prepared_ligands_dir=str(layout["prepared_ligands"]),
        enabled_panels=list(existing.get("enabled_panels", []) or merged_engines),
        pairlist_mode=pairlist_mode or str(existing.get("pairlist_mode", "")),
        has_cocrystal_benchmark_rows=bool(
            existing.get("has_cocrystal_benchmark_rows", False) or has_cocrystal_benchmark_rows
        ),
        source_paths=merged_source_paths,
        favorite_engine=favorite_engine or existing.get("favorite_engine", ""),
        engine_settings=existing.get("engine_settings", {engine: {} for engine in merged_engines}),
        pair_curation_state_file=str(existing.get("pair_curation_state_file", "") or pair_curation_state_path(project_root, profile)),
        latest_pair_round=existing.get("latest_pair_round", ""),
        deployment_root=str(existing.get("deployment_root", "") or deployment_root(project_root, profile)),
        latest_rerun_manifest_file=existing.get("latest_rerun_manifest_file", ""),
        top_pose_selection_policy=str(existing.get("top_pose_selection_policy", "best_affinity") or "best_affinity"),
        top_pose_global_aggregation=str(existing.get("top_pose_global_aggregation", "best_target") or "best_target"),
        hpc_profile=dict(existing.get("hpc_profile", {})),
        notes=merged_notes,
    )
    manifest.engine_settings = {
        **{engine: {} for engine in merged_engines},
        **dict(existing.get("engine_settings", {})),
        **dict(manifest.engine_settings),
    }
    save_manifest(project_root, manifest.to_dict(), layout_profile=profile)
    return {
        "project_root": str(project_root),
        "docking_root": str(layout["docking_root"]),
        "post_docking_root": str(layout["post_docking_root"]),
        "engines": merged_engines,
        "manifest_path": str(manifest_path(project_root, profile)),
        "pairlist_file": str(pairlist_file),
        "layout_type": profile,
        "compatibility_profiles": compatibility_profiles,
        "enabled_panels": manifest.enabled_panels,
    }
