from __future__ import annotations

import argparse
import json
import sys
from dataclasses import fields
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .deployment import generate_condor_deployment, generate_slurm_deployment
from .engine_registry import build_runner, normalize_engines
from .execution_environment import (
    ExecutionEnvironmentConfig,
    apply_environment_runtime,
    normalize_environment,
    validate_environment_selection,
)
from .hpc_profiles import load_hpc_profile, merge_runtime_with_profile
from .models import PairlistRow
from .parameter_schema import (
    apply_schema_to_runtime,
    resolve_parameter_schema,
    transform_pairlist_rows,
    validate_parameter_schema,
    validate_pairlist_geometry,
)
from .preflight import run_docking_preflight
from .preparation.ligand_quality import LIGAND_ADMET_DEFAULTS, audit_project_ligands
from .preparation.receptor_quality import audit_project_receptors
from .preparation.project_builder import DockingPreparationConfig, DockingProjectBuilder
from .project_layout import deployment_root, load_manifest, load_pairlist, save_manifest
from .remote_ops import (
    infer_round_mode_from_deployment_root,
    resolve_remote_target,
    submit_remote_deployment,
    sync_project_to_remote,
)


PAIRLIST_REQUIRED_COLUMNS = [
    "receptor",
    "site_id",
    "ligand",
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
]

RECEPTOR_QC_MIN_ATOM_DEFAULT = 100
RECEPTOR_QC_MIN_HEAVY_ATOM_DEFAULT = 60
RECEPTOR_QC_MIN_CHAIN_DEFAULT = 1
RECEPTOR_QC_MAX_SPAN_DEFAULT = 500.0
LIGAND_QC_MAX_LIPINSKI_DEFAULT = int(LIGAND_ADMET_DEFAULTS.get("max_lipinski_violations", 1))
LIGAND_QC_MAX_MW_DEFAULT = float(LIGAND_ADMET_DEFAULTS.get("max_molecular_weight", 650.0))
LIGAND_QC_MAX_LOGP_DEFAULT = float(LIGAND_ADMET_DEFAULTS.get("max_logp", 6.0))
LIGAND_QC_MAX_TPSA_DEFAULT = float(LIGAND_ADMET_DEFAULTS.get("max_tpsa", 180.0))
LIGAND_QC_MAX_ROTATABLE_DEFAULT = int(LIGAND_ADMET_DEFAULTS.get("max_rotatable_bonds", 15))
LIGAND_QC_MAX_FORMAL_CHARGE_DEFAULT = int(LIGAND_ADMET_DEFAULTS.get("max_formal_charge_abs", 2))
LIGAND_QC_MIN_HEAVY_ATOM_DEFAULT = int(LIGAND_ADMET_DEFAULTS.get("min_heavy_atom_count", 6))


def _pairlist_rows_from_df(df: pd.DataFrame) -> List[PairlistRow]:
    return [
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
            **{item.name: str(row[item.name]) for item in fields(PairlistRow) if item.name not in PAIRLIST_REQUIRED_COLUMNS and item.name in row and pd.notna(row[item.name])},
        )
        for _, row in df.iterrows()
    ]


def _load_pairlist_rows(project_dir: Path) -> List[PairlistRow]:
    return _pairlist_rows_from_df(load_pairlist(project_dir))


def _load_rerun_manifest_rows(rerun_manifest: Path) -> Tuple[List[PairlistRow], List[str]]:
    df = pd.read_csv(rerun_manifest)
    missing = [column for column in PAIRLIST_REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(
            f"Rerun manifest is missing required columns: {', '.join(missing)}"
        )
    engines = []
    if "engine" in df.columns:
        engines = normalize_engines(df["engine"].dropna().astype(str).tolist())
    return _pairlist_rows_from_df(df), engines


def _build_runtime_by_engine(
    args: argparse.Namespace,
    docking_mode: str,
    round_id: str,
) -> Dict[str, Dict[str, object]]:
    return {
        "gnina": {
            "image": args.gnina_image,
            "binary": args.gnina_binary,
            "device": args.gnina_device,
            "cpu": args.gnina_cpu,
            "cnn_scoring": args.gnina_cnn_scoring,
            "score_only": getattr(args, "mode", "standard") == "score-only",
            "exhaustiveness": args.exhaustiveness,
            "num_modes": args.num_modes,
            "use_gpu": True if args.gnina_device is not None else (False if args.gnina_cpu not in (None, 0) else None),
            "round_id": round_id,
            "docking_mode": docking_mode,
        },
        "vina": {
            "conda_env": args.vina_conda_env,
            "binary": args.vina_binary,
            "cpu": args.vina_cpu,
            "exhaustiveness": args.exhaustiveness,
            "num_modes": args.num_modes,
            "round_id": round_id,
            "docking_mode": docking_mode,
        },
        "smina": {
            "conda_env": args.smina_conda_env,
            "binary": args.smina_binary,
            "cpu": args.smina_cpu,
            "scoring": args.smina_scoring,
            "exhaustiveness": args.exhaustiveness,
            "num_modes": args.num_modes,
            "round_id": round_id,
            "docking_mode": docking_mode,
        },
        "autodock4": {
            "binary": args.autodock4_binary,
            "autodock_binary": args.autodock4_binary,
            "autogrid_binary": args.autogrid4_binary,
            "parameter_file": args.autodock4_parameter_file,
            "spacing": args.autodock4_spacing,
            "ga_pop_size": args.autodock4_ga_pop_size,
            "ga_num_evals": args.autodock4_ga_num_evals,
            "ga_num_generations": args.autodock4_ga_num_generations,
            "ga_run": args.autodock4_ga_run,
            "ls_search_freq": args.autodock4_ls_search_freq,
            "torsdof": args.autodock4_torsdof,
            "autodocktools_python": args.autodocktools_python,
            "autodocktools_prepare_gpf4": args.autodocktools_prepare_gpf4,
            "autodocktools_prepare_dpf4": args.autodocktools_prepare_dpf4,
            "round_id": round_id,
            "docking_mode": docking_mode,
        },
    }


def _restore_saved_protocol(args, manifest, argv):
    """Preserve reviewed scientific settings unless the caller overrides them."""
    supplied = {token.split("=", 1)[0] for token in (sys.argv[1:] if argv is None else argv) if token.startswith("--")}
    schema = manifest.get("parameter_schema") or {}
    common = schema.get("common") or {}
    values = {"parameter_mode": schema.get("mode"), "parameter_preset": schema.get("preset"), **common}
    for name, value in values.items():
        if hasattr(args, name) and value is not None and "--" + name.replace("_", "-") not in supplied:
            setattr(args, name, value)
    if supplied & {"--exhaustiveness", "--num-modes"} and not supplied & {"--parameter-mode", "--parameter-preset"}:
        args.parameter_mode = "advanced"
    mappings = {
        "gnina": {"cnn_scoring": "gnina_cnn_scoring"},
        "smina": {"scoring": "smina_scoring"},
        "autodock4": {key: "autodock4_" + key.replace("parameter_file_path", "parameter_file") for key in ["parameter_file", "spacing", "ga_pop_size", "ga_num_evals", "ga_num_generations", "ga_run", "ls_search_freq", "torsdof"]},
    }
    for engine, mapping in mappings.items():
        saved = (manifest.get("engine_settings") or {}).get(engine) or {}
        for key, name in mapping.items():
            if saved.get(key) is not None and "--" + name.replace("_", "-") not in supplied:
                setattr(args, name, saved[key])


def _deployment_manifest_path(project_dir: Path, round_id: str, mode: str) -> Path:
    return deployment_root(project_dir) / round_id / mode / "deployment_manifest.json"


def _load_deployment_manifest(project_dir: Path, round_id: str, mode: str) -> Dict[str, object]:
    path = _deployment_manifest_path(project_dir, round_id, mode)
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return payload if isinstance(payload, dict) else {}


def _validate_exhaustive_deploy(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    rerun_manifest_file: Optional[Path],
) -> None:
    if args.mode != "exhaustive":
        return
    if rerun_manifest_file:
        return
    if args.allow_full_exhaustive:
        print("⚠️  Full-pair exhaustive deployment override enabled. This bypasses the normal rerun-manifest workflow.")
        return
    parser.error(
        "--mode exhaustive requires --from-rerun-manifest. "
        "Use --allow-full-exhaustive only when you intentionally want to rerun the entire pairlist at exhaustive settings."
    )


def _validate_exhaustive_submit(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    project_dir: Path,
    round_id: str,
    mode: str,
) -> None:
    if mode != "exhaustive":
        return
    deployment_manifest = _load_deployment_manifest(project_dir, round_id, mode)
    rerun_manifest_file = str(deployment_manifest.get("rerun_manifest_file", "") or "").strip()
    if rerun_manifest_file:
        return
    if args.allow_full_exhaustive:
        print("⚠️  Full-pair exhaustive submit override enabled. This bundle does not record a rerun manifest.")
        return
    manifest_path = _deployment_manifest_path(project_dir, round_id, mode)
    if deployment_manifest:
        parser.error(
            f"Refusing to submit exhaustive deployment '{round_id}' because {manifest_path} does not record a rerun manifest. "
            "Regenerate from analyze comparative or pass --allow-full-exhaustive to override."
        )
    parser.error(
        f"Cannot verify exhaustive deployment '{round_id}' because {manifest_path} was not found locally. "
        "Regenerate the bundle locally first or pass --allow-full-exhaustive to override."
    )


def _ligand_validation_report_path(project_dir: Path) -> Path:
    manifest = load_manifest(project_dir)
    metadata_dir = Path(manifest.get("metadata_dir") or (project_dir / "metadata"))
    return metadata_dir / "ligand_validation_report.json"


def _receptor_validation_report_path(project_dir: Path) -> Path:
    manifest = load_manifest(project_dir)
    metadata_dir = Path(manifest.get("metadata_dir") or (project_dir / "metadata"))
    return metadata_dir / "receptor_validation_report.json"


def _print_ligand_validation_failure(payload: Dict[str, object]) -> None:
    report_path = str(payload.get("report_path", "") or "")
    print("❌ Ligand validation failed")
    if report_path:
        print(f"   Report: {report_path}")
    print(f"   Invalid ligands: {', '.join(payload.get('affected_ligands', []))}")
    thresholds = payload.get("admet_thresholds") or {}
    if thresholds and bool(payload.get("admet_filters_enabled", False)):
        print(
            "   ADMET thresholds: "
            f"lipinski<={thresholds.get('max_lipinski_violations')} "
            f"MW<={thresholds.get('max_molecular_weight')} "
            f"logP<={thresholds.get('max_logp')} "
            f"TPSA<={thresholds.get('max_tpsa')} "
            f"rotors<={thresholds.get('max_rotatable_bonds')} "
            f"|charge|<={thresholds.get('max_formal_charge_abs')} "
            f"heavy_atoms>={thresholds.get('min_heavy_atom_count')}"
        )
    for issue in payload.get("issues", [])[:10]:
        affected_pairs = int(issue.get("affected_pairs", 0) or 0)
        print(
            f"   {issue['ligand']}: {issue['reason']} "
            f"(affected_pairs={affected_pairs})"
        )
        if issue.get("details"):
            print(f"      details: {issue['details']}")
        if issue.get("raw_pdb_file"):
            print(f"      raw_pdb: {issue['raw_pdb_file']}")
    remaining = int(payload.get("issue_count", 0) or 0) - min(len(payload.get("issues", [])), 10)
    if remaining > 0:
        print(f"   ... {remaining} additional issue(s) omitted")
    if int(payload.get("warning_count", 0) or 0) > 0:
        print(f"   Warnings: {int(payload.get('warning_count', 0) or 0)} (see report for details)")
    print("   Refusing to continue until the invalid ligand files are repaired or removed from the active pair set.")


def _print_receptor_validation_failure(payload: Dict[str, object]) -> None:
    report_path = str(payload.get("report_path", "") or "")
    print("❌ Receptor QC gate failed")
    if report_path:
        print(f"   Report: {report_path}")
    print(f"   Invalid receptors: {', '.join(payload.get('affected_receptors', []))}")
    thresholds = payload.get("thresholds") or {}
    if thresholds:
        print(
            "   Thresholds: "
            f"min_atom_count={thresholds.get('min_atom_count')} "
            f"min_heavy_atom_count={thresholds.get('min_heavy_atom_count')} "
            f"min_chain_count={thresholds.get('min_chain_count')} "
            f"max_coordinate_span={thresholds.get('max_coordinate_span')}"
        )
    for issue in payload.get("issues", [])[:10]:
        affected_pairs = int(issue.get("affected_pairs", 0) or 0)
        print(
            f"   {issue['receptor']}: {issue['reason']} "
            f"(affected_pairs={affected_pairs})"
        )
        if issue.get("details"):
            print(f"      details: {issue['details']}")
    remaining = int(payload.get("issue_count", 0) or 0) - min(len(payload.get("issues", [])), 10)
    if remaining > 0:
        print(f"   ... {remaining} additional issue(s) omitted")
    print("   Fix receptor preparation or relax QC thresholds before docking.")


def _audit_project_ligands_for_rows(
    project_dir: Path,
    rows: List[PairlistRow],
    *,
    enable_admet_filters: bool,
    max_lipinski_violations: int,
    max_molecular_weight: float,
    max_logp: float,
    max_tpsa: float,
    max_rotatable_bonds: int,
    max_formal_charge_abs: int,
    min_heavy_atom_count: int,
    block_pains: bool,
    block_brenk: bool,
    block_reactive: bool,
) -> Dict[str, object]:
    return audit_project_ligands(
        project_root=project_dir,
        rows=rows,
        enable_admet_filters=enable_admet_filters,
        admet_thresholds={
            "max_lipinski_violations": int(max_lipinski_violations),
            "max_molecular_weight": float(max_molecular_weight),
            "max_logp": float(max_logp),
            "max_tpsa": float(max_tpsa),
            "max_rotatable_bonds": int(max_rotatable_bonds),
            "max_formal_charge_abs": int(max_formal_charge_abs),
            "min_heavy_atom_count": int(min_heavy_atom_count),
            "block_pains": bool(block_pains),
            "block_brenk": bool(block_brenk),
            "block_reactive": bool(block_reactive),
        },
        report_path=_ligand_validation_report_path(project_dir),
    )


def _audit_project_receptors_for_rows(
    project_dir: Path,
    rows: List[PairlistRow],
    *,
    min_atom_count: int,
    min_heavy_atom_count: int,
    min_chain_count: int,
    max_coordinate_span: float,
) -> Dict[str, object]:
    return audit_project_receptors(
        project_root=project_dir,
        rows=rows,
        min_atom_count=min_atom_count,
        min_heavy_atom_count=min_heavy_atom_count,
        min_chain_count=min_chain_count,
        max_coordinate_span=max_coordinate_span,
        report_path=_receptor_validation_report_path(project_dir),
    )


def prepare_docking_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Build a docking-ready multi-engine project")
    parser.add_argument("--prepared-proteins", required=True, help="Directory of prepared receptor files")
    parser.add_argument("--prepared-ligands", required=True, help="Directory of prepared ligand files")
    parser.add_argument("--pair-intent", help="Manual pair-intent CSV (required when --pair-mode manual)")
    parser.add_argument("--pairlist-file", help="Prebuilt pairlist.csv to materialize into the docking project")
    parser.add_argument("--excel", help="multi_pdb_analysis.xlsx path (required unless --pairlist-file is provided)")
    parser.add_argument("--output", required=True, help="Output project directory")
    parser.add_argument("--raw-proteins", help="Optional raw protein directory to record in manifest")
    parser.add_argument("--raw-ligands", help="Optional raw ligand directory to record in manifest")
    parser.add_argument("--default-site-id", default="site_1", help="Default site id when pair intent omits it")
    parser.add_argument("--default-box-size", type=float, default=20.0, help="Default cubic box size")
    parser.add_argument("--asset-mode", choices=["symlink", "copy"], default="symlink", help="How to materialize receptors/ligands into the project")
    parser.add_argument("--engines", default="gnina,vina,smina,autodock4", help="Comma-separated engine list to initialize")
    parser.add_argument("--project-name", help="Optional manifest project name")
    parser.add_argument("--layout-profile", choices=["canonical", "docking_legacy"], default="docking_legacy")
    parser.add_argument(
        "--pair-mode",
        choices=["manual", "protein-based", "cocrystal_only", "cocrystal_plus_all", "curated_cartesian", "curated_per_protein"],
        default="manual",
        help="How to generate pairlist rows",
    )
    parser.add_argument(
        "--dock-all-proteins-to-all-ligands",
        action="store_true",
        help="In protein-based mode, generate the full protein x ligand Cartesian product instead of only reference ligands",
    )

    args = parser.parse_args(argv)
    if args.pair_mode == "manual" and not args.pair_intent and not args.pairlist_file:
        parser.error("--pair-intent is required when --pair-mode manual")
    if not args.excel and not args.pairlist_file:
        parser.error("--excel is required unless --pairlist-file is provided")
    if args.dock_all_proteins_to_all_ligands and args.pair_mode != "protein-based":
        parser.error("--dock-all-proteins-to-all-ligands is only valid with --pair-mode protein-based")

    config = DockingPreparationConfig(
        prepared_proteins=Path(args.prepared_proteins).expanduser(),
        prepared_ligands=Path(args.prepared_ligands).expanduser(),
        pair_intent=Path(args.pair_intent).expanduser() if args.pair_intent else None,
        excel_path=Path(args.excel).expanduser() if args.excel else None,
        output_dir=Path(args.output).expanduser(),
        pairlist_file=Path(args.pairlist_file).expanduser() if args.pairlist_file else None,
        raw_proteins=Path(args.raw_proteins).expanduser() if args.raw_proteins else None,
        raw_ligands=Path(args.raw_ligands).expanduser() if args.raw_ligands else None,
        default_site_id=args.default_site_id,
        default_box_size=args.default_box_size,
        asset_mode=args.asset_mode,
        engines=normalize_engines([args.engines]),
        project_name=args.project_name or "",
        pair_mode=args.pair_mode,
        dock_all_proteins_to_all_ligands=args.dock_all_proteins_to_all_ligands,
        layout_profile=args.layout_profile,
    )
    summary = DockingProjectBuilder(config).build()
    print("✅ Docking project created")
    print(f"   Project root: {summary['project_root']}")
    print(f"   Pair mode: {summary['pair_mode']}")
    print(f"   Pairlist rows: {summary['pair_count']}")
    print(f"   Receptors: {summary['receptor_count']}")
    print(f"   Ligands: {summary['ligand_count']}")
    print(f"   Engines: {', '.join(summary['engines'])}")
    if summary.get("warning_count"):
        print(f"   Warnings: {summary['warning_count']} (see metadata/preparation_warnings.txt)")
    return 0


def dock_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Run GNINA, Vina, Smina, and/or AutoDock4 on a prepared docking project")
    parser.add_argument("--project-dir", required=True, help="Docking project root")
    parser.add_argument("--engines", help="Comma-separated engine list (default: project manifest engines)")
    parser.add_argument("--mode", choices=["standard", "score-only"], default="standard", help="Docking mode")
    parser.add_argument("--favorite-engine", help="Optional favorite engine to record in the manifest")
    parser.add_argument("--dry-run", action="store_true", help="Print commands and write manifests without executing engines")
    parser.add_argument("--skip-completed", action="store_true", help="Skip rows with existing pose output")
    parser.add_argument("--gnina-image", help="Apptainer/Singularity image for GNINA")
    parser.add_argument("--gnina-binary", help="GNINA binary path inside or outside the container")
    parser.add_argument("--gnina-device", help="GNINA device id")
    parser.add_argument("--gnina-cpu", type=int, help="CPU cores for GNINA CPU mode")
    parser.add_argument("--gnina-cnn-scoring", default="rescore", help="GNINA CNN scoring mode")
    parser.add_argument("--vina-conda-env", help="Conda environment name for AutoDock Vina")
    parser.add_argument("--vina-binary", help="Vina binary name/path")
    parser.add_argument("--vina-cpu", type=int, help="CPU cores for Vina")
    parser.add_argument("--smina-conda-env", help="Conda environment name for Smina")
    parser.add_argument("--smina-binary", help="Smina binary name/path")
    parser.add_argument("--smina-cpu", type=int, help="CPU cores for Smina")
    parser.add_argument("--smina-scoring", help="Optional Smina scoring function")
    parser.add_argument("--autodock4-binary", help="AutoDock4 binary name/path")
    parser.add_argument("--autogrid4-binary", help="AutoGrid4 binary name/path")
    parser.add_argument("--autodock4-parameter-file", default="AD4.1_bound.dat", help="AutoDock4 parameter file")
    parser.add_argument("--autodock4-spacing", type=float, default=0.375, help="AutoGrid spacing (A)")
    parser.add_argument("--autodock4-ga-pop-size", type=int, default=150, help="AutoDock4 GA population size")
    parser.add_argument("--autodock4-ga-num-evals", type=int, default=25000000, help="AutoDock4 GA maximum evaluations")
    parser.add_argument("--autodock4-ga-num-generations", type=int, default=27000, help="AutoDock4 GA maximum generations")
    parser.add_argument("--autodock4-ga-run", type=int, default=100, help="AutoDock4 GA run count")
    parser.add_argument("--autodock4-ls-search-freq", type=float, default=0.06, help="AutoDock4 local-search frequency")
    parser.add_argument("--autodock4-torsdof", type=int, default=0, help="Override AutoDock4 torsional degrees of freedom (0 keeps default)")
    parser.add_argument("--autodocktools-python", help="Python executable used to run AutoDockTools prepare_gpf4.py/prepare_dpf4.py")
    parser.add_argument("--autodocktools-prepare-gpf4", help="Path to AutoDockTools prepare_gpf4.py script")
    parser.add_argument("--autodocktools-prepare-dpf4", help="Path to AutoDockTools prepare_dpf4.py script")
    parser.add_argument("--exhaustiveness", type=int, default=16, help="Shared exhaustiveness default")
    parser.add_argument("--num-modes", type=int, default=20, help="Shared num_modes default")
    parser.add_argument("--seed", type=int, help="Shared deterministic seed where supported by selected engines")
    parser.add_argument("--box-scale", type=float, default=1.0, help="Scale factor applied to pairlist box sizes")
    parser.add_argument("--box-padding", type=float, default=0.0, help="Padding added to each side of docking box sizes (Angstrom)")
    parser.add_argument("--parameter-mode", choices=["basic", "advanced"], default="basic", help="Parameter profile mode")
    parser.add_argument(
        "--parameter-preset",
        choices=["screening_fast", "balanced", "exhaustive"],
        default="balanced",
        help="Basic-mode preset for shared docking parameters",
    )
    parser.add_argument(
        "--execution-environment",
        choices=["local_cpu", "local_gpu", "conda_env", "container", "remote_hpc"],
        default="local_cpu",
        help="Execution environment profile for docking runtime orchestration",
    )
    parser.add_argument("--execution-workdir", help="Optional execution working directory for local docking processes")
    parser.add_argument("--prerequisites-dir", help="Optional directory containing prerequisite assets/modules for docking")
    parser.add_argument("--required-files-dir", help="Optional required-file directory that must exist before launch")
    parser.add_argument("--shared-conda-env", help="Shared conda environment for Vina/Smina in conda_env mode")
    parser.add_argument("--container-image", help="Fallback container image for container mode (used when engine image is absent)")
    parser.add_argument("--no-ligand-qc-gate", action="store_true", help="Disable ligand structural QC gate before docking")
    parser.add_argument("--no-receptor-qc-gate", action="store_true", help="Disable receptor structural QC gate before docking")
    parser.add_argument("--no-ligand-admet-filters", action="store_true", help="Disable ADMET-based ligand filters inside ligand QC gate")
    parser.add_argument("--ligand-qc-max-lipinski-violations", type=int, default=LIGAND_QC_MAX_LIPINSKI_DEFAULT, help="Maximum allowed Lipinski violations per ligand")
    parser.add_argument("--ligand-qc-max-molecular-weight", type=float, default=LIGAND_QC_MAX_MW_DEFAULT, help="Maximum allowed molecular weight (Da)")
    parser.add_argument("--ligand-qc-max-logp", type=float, default=LIGAND_QC_MAX_LOGP_DEFAULT, help="Maximum allowed LogP")
    parser.add_argument("--ligand-qc-max-tpsa", type=float, default=LIGAND_QC_MAX_TPSA_DEFAULT, help="Maximum allowed topological polar surface area")
    parser.add_argument("--ligand-qc-max-rotatable-bonds", type=int, default=LIGAND_QC_MAX_ROTATABLE_DEFAULT, help="Maximum allowed rotatable bonds")
    parser.add_argument("--ligand-qc-max-formal-charge-abs", type=int, default=LIGAND_QC_MAX_FORMAL_CHARGE_DEFAULT, help="Maximum allowed absolute formal charge")
    parser.add_argument("--ligand-qc-min-heavy-atom-count", type=int, default=LIGAND_QC_MIN_HEAVY_ATOM_DEFAULT, help="Minimum required heavy-atom count")
    parser.add_argument("--ligand-qc-allow-pains", action="store_true", help="Do not block ligands with PAINS alerts")
    parser.add_argument("--ligand-qc-allow-brenk", action="store_true", help="Do not block ligands with Brenk alerts")
    parser.add_argument("--ligand-qc-allow-reactive", action="store_true", help="Do not block ligands with reactive SMARTS alerts")
    parser.add_argument("--receptor-qc-min-atom-count", type=int, default=100, help="Minimum receptor atom count for QC gate")
    parser.add_argument("--receptor-qc-min-heavy-atom-count", type=int, default=60, help="Minimum receptor heavy-atom count for QC gate")
    parser.add_argument("--receptor-qc-min-chain-count", type=int, default=1, help="Minimum receptor chain count for QC gate")
    parser.add_argument("--receptor-qc-max-coordinate-span", type=float, default=500.0, help="Maximum receptor coordinate span (A) for QC gate")

    args = parser.parse_args(argv)
    project_dir = Path(args.project_dir).expanduser()
    manifest = load_manifest(project_dir)
    _restore_saved_protocol(args, manifest, argv)
    round_id = str(manifest.get("latest_pair_round") or "round_001")
    engines = normalize_engines([args.engines]) if args.engines else normalize_engines(manifest.get("engines", []))
    rows = _load_pairlist_rows(project_dir)
    runtime_by_engine = _build_runtime_by_engine(args, docking_mode="screen", round_id=round_id)

    env_config = ExecutionEnvironmentConfig(
        environment=normalize_environment(args.execution_environment),
        workdir=str(args.execution_workdir or "").strip(),
        prerequisites_dir=str(args.prerequisites_dir or "").strip(),
        required_files_dir=str(args.required_files_dir or "").strip(),
        shared_conda_env=str(args.shared_conda_env or "").strip(),
        container_image=str(args.container_image or "").strip(),
    )
    env_errors, env_warnings = validate_environment_selection(engines, env_config, runtime_by_engine)
    for warning in env_warnings:
        print(f"⚠️  {warning}")
    if env_errors:
        print("❌ Docking environment selection is invalid:")
        for issue in env_errors:
            print(f"   - {issue}")
        return 1
    runtime_by_engine = apply_environment_runtime(engines, runtime_by_engine, env_config)

    parameter_schema = resolve_parameter_schema(
        mode=args.parameter_mode,
        preset=args.parameter_preset,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        seed=args.seed,
        box_scale=args.box_scale,
        box_padding=args.box_padding,
        runtime_by_engine=runtime_by_engine,
    )
    schema_errors, schema_warnings = validate_parameter_schema(parameter_schema, engines)
    for warning in schema_warnings:
        print(f"⚠️  {warning}")
    if schema_errors:
        print("❌ Docking parameter schema is invalid:")
        for issue in schema_errors:
            print(f"   - {issue}")
        return 1
    runtime_by_engine = apply_schema_to_runtime(parameter_schema, runtime_by_engine)
    transformed_rows = transform_pairlist_rows(rows, parameter_schema)

    preflight = run_docking_preflight(
        project_dir=project_dir,
        engines=engines,
        pairlist_rows=transformed_rows,
        runtime_by_engine=runtime_by_engine,
        execution_workdir=env_config.workdir,
        prerequisites_dir=env_config.prerequisites_dir,
        required_files_dir=env_config.required_files_dir,
        enable_ligand_qc=not bool(args.no_ligand_qc_gate),
        enable_receptor_qc=not bool(args.no_receptor_qc_gate),
        ligand_qc_report_path=str(_ligand_validation_report_path(project_dir)),
        receptor_qc_report_path=str(_receptor_validation_report_path(project_dir)),
        ligand_admet_filters_enabled=not bool(args.no_ligand_admet_filters),
        ligand_max_lipinski_violations=int(args.ligand_qc_max_lipinski_violations),
        ligand_max_molecular_weight=float(args.ligand_qc_max_molecular_weight),
        ligand_max_logp=float(args.ligand_qc_max_logp),
        ligand_max_tpsa=float(args.ligand_qc_max_tpsa),
        ligand_max_rotatable_bonds=int(args.ligand_qc_max_rotatable_bonds),
        ligand_max_formal_charge_abs=int(args.ligand_qc_max_formal_charge_abs),
        ligand_min_heavy_atom_count=int(args.ligand_qc_min_heavy_atom_count),
        ligand_block_pains=not bool(args.ligand_qc_allow_pains),
        ligand_block_brenk=not bool(args.ligand_qc_allow_brenk),
        ligand_block_reactive=not bool(args.ligand_qc_allow_reactive),
        receptor_min_atom_count=int(args.receptor_qc_min_atom_count),
        receptor_min_heavy_atom_count=int(args.receptor_qc_min_heavy_atom_count),
        receptor_min_chain_count=int(args.receptor_qc_min_chain_count),
        receptor_max_coordinate_span=float(args.receptor_qc_max_coordinate_span),
        dry_run=bool(args.dry_run),
    )
    for warning in preflight.get("warnings", []):
        print(f"⚠️  {warning}")
    if not preflight.get("ok", False):
        print("❌ Docking preflight blocked launch")
        preflight_details = preflight.get("details", {}) or {}
        ligand_qc = preflight_details.get("ligand_qc")
        receptor_qc = preflight_details.get("receptor_qc")
        if isinstance(ligand_qc, dict) and int(ligand_qc.get("issue_count", 0) or 0) > 0:
            _print_ligand_validation_failure(ligand_qc)
        if isinstance(receptor_qc, dict) and int(receptor_qc.get("issue_count", 0) or 0) > 0:
            _print_receptor_validation_failure(receptor_qc)
        for issue in preflight.get("errors", []):
            print(f"   - {issue}")
        return 1
    print(f"✅ Docking preflight passed for {len(transformed_rows)} pair(s)")

    results = {}
    for engine in engines:
        runner = build_runner(engine, project_dir, runtime_by_engine[engine])
        results[engine] = runner.run(transformed_rows, dry_run=args.dry_run, skip_completed=args.skip_completed)
        completed = sum(1 for job in results[engine]["jobs"] if job["status"] == "completed")
        dry = sum(1 for job in results[engine]["jobs"] if job["status"] == "dry_run")
        skipped = sum(1 for job in results[engine]["jobs"] if job["status"] == "skipped")
        failed = sum(1 for job in results[engine]["jobs"] if job["status"] == "failed")
        print(f"🔧 {engine}: completed={completed} dry_run={dry} skipped={skipped} failed={failed}")

    if args.favorite_engine:
        manifest["favorite_engine"] = args.favorite_engine
    manifest.setdefault("engine_settings", {}).update(runtime_by_engine)
    manifest["execution_settings"] = env_config.to_dict()
    manifest["parameter_schema"] = parameter_schema.to_dict()
    manifest["qc_settings"] = {
        "ligand_qc_gate_enabled": not bool(args.no_ligand_qc_gate),
        "ligand_admet_filters_enabled": not bool(args.no_ligand_admet_filters),
        "ligand_max_lipinski_violations": int(args.ligand_qc_max_lipinski_violations),
        "ligand_max_molecular_weight": float(args.ligand_qc_max_molecular_weight),
        "ligand_max_logp": float(args.ligand_qc_max_logp),
        "ligand_max_tpsa": float(args.ligand_qc_max_tpsa),
        "ligand_max_rotatable_bonds": int(args.ligand_qc_max_rotatable_bonds),
        "ligand_max_formal_charge_abs": int(args.ligand_qc_max_formal_charge_abs),
        "ligand_min_heavy_atom_count": int(args.ligand_qc_min_heavy_atom_count),
        "ligand_block_pains": not bool(args.ligand_qc_allow_pains),
        "ligand_block_brenk": not bool(args.ligand_qc_allow_brenk),
        "ligand_block_reactive": not bool(args.ligand_qc_allow_reactive),
        "receptor_qc_gate_enabled": not bool(args.no_receptor_qc_gate),
        "receptor_min_atom_count": int(args.receptor_qc_min_atom_count),
        "receptor_min_heavy_atom_count": int(args.receptor_qc_min_heavy_atom_count),
        "receptor_min_chain_count": int(args.receptor_qc_min_chain_count),
        "receptor_max_coordinate_span": float(args.receptor_qc_max_coordinate_span),
    }
    manifest["latest_pair_round"] = round_id
    save_manifest(project_dir, manifest)
    return 1 if any(job["status"] == "failed" for result in results.values() for job in result["jobs"]) else 0


def deploy_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Generate Slurm deployment assets for a prepared docking project")
    parser.add_argument("--project-dir", required=True, help="Docking project root")
    parser.add_argument("--engines", help="Comma-separated engine list (default: project manifest engines or rerun manifest engine)")
    parser.add_argument("--round", help="Logical round identifier for the deployment bundle")
    parser.add_argument("--mode", choices=["screen", "exhaustive"], default="screen", help="Deployment stage mode")
    parser.add_argument("--from-rerun-manifest", help="Optional rerun-manifest CSV generated by comparative analysis")
    parser.add_argument(
        "--allow-full-exhaustive",
        action="store_true",
        help="Allow exhaustive deployment without a rerun manifest. Use only for deliberate full-pair exhaustive campaigns.",
    )
    parser.add_argument("--skip-completed", action="store_true", help="Skip rows with existing pose output")
    parser.add_argument("--gnina-image", help="Apptainer/Singularity image for GNINA")
    parser.add_argument("--gnina-binary", help="GNINA binary path inside or outside the container")
    parser.add_argument("--gnina-device", help="GNINA device id")
    parser.add_argument("--gnina-cpu", type=int, help="CPU cores for GNINA CPU mode")
    parser.add_argument("--gnina-cnn-scoring", default="rescore", help="GNINA CNN scoring mode")
    parser.add_argument("--vina-conda-env", help="Conda environment name for AutoDock Vina")
    parser.add_argument("--vina-binary", help="Vina binary name/path")
    parser.add_argument("--vina-cpu", type=int, help="CPU cores for Vina")
    parser.add_argument("--smina-conda-env", help="Conda environment name for Smina")
    parser.add_argument("--smina-binary", help="Smina binary name/path")
    parser.add_argument("--smina-cpu", type=int, help="CPU cores for Smina")
    parser.add_argument("--smina-scoring", help="Optional Smina scoring function")
    parser.add_argument("--autodock4-binary", help="AutoDock4 binary name/path")
    parser.add_argument("--autogrid4-binary", help="AutoGrid4 binary name/path")
    parser.add_argument("--autodock4-parameter-file", default="AD4.1_bound.dat", help="AutoDock4 parameter file")
    parser.add_argument("--autodock4-spacing", type=float, default=0.375, help="AutoGrid spacing (A)")
    parser.add_argument("--autodock4-ga-pop-size", type=int, default=150, help="AutoDock4 GA population size")
    parser.add_argument("--autodock4-ga-num-evals", type=int, default=25000000, help="AutoDock4 GA maximum evaluations")
    parser.add_argument("--autodock4-ga-num-generations", type=int, default=27000, help="AutoDock4 GA maximum generations")
    parser.add_argument("--autodock4-ga-run", type=int, default=100, help="AutoDock4 GA run count")
    parser.add_argument("--autodock4-ls-search-freq", type=float, default=0.06, help="AutoDock4 local-search frequency")
    parser.add_argument("--autodock4-torsdof", type=int, default=0, help="Override AutoDock4 torsional degrees of freedom (0 keeps default)")
    parser.add_argument("--autodocktools-python", help="Python executable used to run AutoDockTools prepare_gpf4.py/prepare_dpf4.py")
    parser.add_argument("--autodocktools-prepare-gpf4", help="Path to AutoDockTools prepare_gpf4.py script")
    parser.add_argument("--autodocktools-prepare-dpf4", help="Path to AutoDockTools prepare_dpf4.py script")
    parser.add_argument("--exhaustiveness", type=int, default=16, help="Shared exhaustiveness default")
    parser.add_argument("--num-modes", type=int, default=20, help="Shared num_modes default")
    parser.add_argument("--seed", type=int, help="Shared deterministic seed")
    parser.add_argument("--box-scale", type=float, default=1.0)
    parser.add_argument("--box-padding", type=float, default=0.0)
    parser.add_argument("--parameter-mode", choices=["basic", "advanced"], default="basic")
    parser.add_argument("--parameter-preset", choices=["screening_fast", "balanced", "exhaustive"], default="balanced")
    parser.add_argument("--slurm-time", help="SBATCH --time value")
    parser.add_argument("--slurm-mem", help="SBATCH --mem value")
    parser.add_argument("--slurm-cpus-per-task", type=int, help="SBATCH --cpus-per-task value")
    parser.add_argument("--slurm-gpus", type=int, help="SBATCH --gpus value override")
    parser.add_argument("--slurm-partition", help="SBATCH --partition value")
    parser.add_argument("--slurm-account", help="SBATCH --account value")
    parser.add_argument("--slurm-job-name-prefix", help="Prefix for Slurm job names")
    parser.add_argument(
        "--slurm-submit-mode",
        choices=["slurm_array", "single_job"],
        help="Slurm submission strategy: arrays per pair or one engine job that loops through the full manifest",
    )
    parser.add_argument("--slurm-array-parallelism", type=int, help="Optional cap for concurrent tasks in the generated Slurm array")
    parser.add_argument("--slurm-extra-args", help="Extra SBATCH directives separated by ';'")
    parser.add_argument("--condor-cpus", type=int, help="HTCondor request_cpus override")
    parser.add_argument("--condor-memory", help="HTCondor request_memory override (e.g. 32GB)")
    parser.add_argument("--condor-disk", help="HTCondor request_disk override (e.g. 20GB)")
    parser.add_argument("--condor-gpus", type=int, help="HTCondor request_gpus override")
    parser.add_argument("--condor-submit-mode", choices=["condor_array", "single_job"], help="HTCondor submission strategy")
    parser.add_argument("--condor-parallelism", type=int, help="max_materialize cap for condor_array (0 = unlimited)")
    parser.add_argument("--condor-requirements", help="HTCondor requirements expression")
    parser.add_argument("--condor-job-name-prefix", help="Prefix for HTCondor job names")
    parser.add_argument("--hpc-profile", help="Builtin or local HPC deployment profile name")
    parser.add_argument("--hpc-profile-file", help="Explicit YAML/JSON file with local HPC deployment settings")
    parser.add_argument("--remote-project-dir", help="Execution path of the project on the HPC so generated scripts are portable after transfer")
    parser.add_argument("--no-ligand-qc-gate", action="store_true", help="Disable ligand structural QC gate before deployment generation")
    parser.add_argument("--no-receptor-qc-gate", action="store_true", help="Disable receptor structural QC gate before deployment generation")
    parser.add_argument("--no-ligand-admet-filters", action="store_true", help="Disable ADMET-based ligand filters inside ligand QC gate")
    parser.add_argument("--ligand-qc-max-lipinski-violations", type=int, default=LIGAND_QC_MAX_LIPINSKI_DEFAULT, help="Maximum allowed Lipinski violations per ligand")
    parser.add_argument("--ligand-qc-max-molecular-weight", type=float, default=LIGAND_QC_MAX_MW_DEFAULT, help="Maximum allowed molecular weight (Da)")
    parser.add_argument("--ligand-qc-max-logp", type=float, default=LIGAND_QC_MAX_LOGP_DEFAULT, help="Maximum allowed LogP")
    parser.add_argument("--ligand-qc-max-tpsa", type=float, default=LIGAND_QC_MAX_TPSA_DEFAULT, help="Maximum allowed topological polar surface area")
    parser.add_argument("--ligand-qc-max-rotatable-bonds", type=int, default=LIGAND_QC_MAX_ROTATABLE_DEFAULT, help="Maximum allowed rotatable bonds")
    parser.add_argument("--ligand-qc-max-formal-charge-abs", type=int, default=LIGAND_QC_MAX_FORMAL_CHARGE_DEFAULT, help="Maximum allowed absolute formal charge")
    parser.add_argument("--ligand-qc-min-heavy-atom-count", type=int, default=LIGAND_QC_MIN_HEAVY_ATOM_DEFAULT, help="Minimum required heavy-atom count")
    parser.add_argument("--ligand-qc-allow-pains", action="store_true", help="Do not block ligands with PAINS alerts")
    parser.add_argument("--ligand-qc-allow-brenk", action="store_true", help="Do not block ligands with Brenk alerts")
    parser.add_argument("--ligand-qc-allow-reactive", action="store_true", help="Do not block ligands with reactive SMARTS alerts")
    parser.add_argument("--receptor-qc-min-atom-count", type=int, default=RECEPTOR_QC_MIN_ATOM_DEFAULT, help="Minimum receptor atom count for QC gate")
    parser.add_argument("--receptor-qc-min-heavy-atom-count", type=int, default=RECEPTOR_QC_MIN_HEAVY_ATOM_DEFAULT, help="Minimum receptor heavy-atom count for QC gate")
    parser.add_argument("--receptor-qc-min-chain-count", type=int, default=RECEPTOR_QC_MIN_CHAIN_DEFAULT, help="Minimum receptor chain count for QC gate")
    parser.add_argument("--receptor-qc-max-coordinate-span", type=float, default=RECEPTOR_QC_MAX_SPAN_DEFAULT, help="Maximum receptor coordinate span (A) for QC gate")

    args = parser.parse_args(argv)
    project_dir = Path(args.project_dir).expanduser().resolve()
    manifest = load_manifest(project_dir)
    _restore_saved_protocol(args, manifest, argv)
    rerun_manifest_file = Path(args.from_rerun_manifest).expanduser().resolve() if args.from_rerun_manifest else None
    if rerun_manifest_file and not rerun_manifest_file.exists():
        parser.error(f"Rerun manifest not found: {rerun_manifest_file}")
    _validate_exhaustive_deploy(parser, args, rerun_manifest_file)
    if rerun_manifest_file:
        rows, rerun_engines = _load_rerun_manifest_rows(rerun_manifest_file)
    else:
        rows = _load_pairlist_rows(project_dir)
        rerun_engines = []
    if not bool(args.no_ligand_qc_gate):
        ligand_audit = _audit_project_ligands_for_rows(
            project_dir,
            rows,
            enable_admet_filters=not bool(args.no_ligand_admet_filters),
            max_lipinski_violations=int(args.ligand_qc_max_lipinski_violations),
            max_molecular_weight=float(args.ligand_qc_max_molecular_weight),
            max_logp=float(args.ligand_qc_max_logp),
            max_tpsa=float(args.ligand_qc_max_tpsa),
            max_rotatable_bonds=int(args.ligand_qc_max_rotatable_bonds),
            max_formal_charge_abs=int(args.ligand_qc_max_formal_charge_abs),
            min_heavy_atom_count=int(args.ligand_qc_min_heavy_atom_count),
            block_pains=not bool(args.ligand_qc_allow_pains),
            block_brenk=not bool(args.ligand_qc_allow_brenk),
            block_reactive=not bool(args.ligand_qc_allow_reactive),
        )
        if ligand_audit.get("issue_count"):
            _print_ligand_validation_failure(ligand_audit)
            return 1
    if not bool(args.no_receptor_qc_gate):
        receptor_audit = _audit_project_receptors_for_rows(
            project_dir,
            rows,
            min_atom_count=int(args.receptor_qc_min_atom_count),
            min_heavy_atom_count=int(args.receptor_qc_min_heavy_atom_count),
            min_chain_count=int(args.receptor_qc_min_chain_count),
            max_coordinate_span=float(args.receptor_qc_max_coordinate_span),
        )
        if receptor_audit.get("issue_count"):
            _print_receptor_validation_failure(receptor_audit)
            return 1

    engines = normalize_engines([args.engines]) if args.engines else normalize_engines(rerun_engines or manifest.get("engines", []))
    round_id = args.round or str(manifest.get("latest_pair_round") or "round_001")
    try:
        hpc_profile = load_hpc_profile(
            project_dir,
            profile_name=args.hpc_profile or "",
            profile_file=args.hpc_profile_file or "",
        )
    except ValueError as exc:
        parser.error(str(exc))
    runtime_by_engine = merge_runtime_with_profile(
        _build_runtime_by_engine(args, docking_mode=args.mode, round_id=round_id),
        hpc_profile,
    )
    parameter_schema = resolve_parameter_schema(
        mode=args.parameter_mode, preset=args.parameter_preset,
        exhaustiveness=args.exhaustiveness, num_modes=args.num_modes,
        seed=args.seed, box_scale=args.box_scale, box_padding=args.box_padding,
        runtime_by_engine=runtime_by_engine,
    )
    schema_errors, schema_warnings = validate_parameter_schema(parameter_schema, engines)
    rows = transform_pairlist_rows(rows, parameter_schema)
    schema_errors.extend(validate_pairlist_geometry(rows))
    if schema_errors:
        for issue in schema_errors:
            print(f"❌ {issue}")
        return 1
    runtime_by_engine = apply_schema_to_runtime(parameter_schema, runtime_by_engine)
    manifest["parameter_schema"] = parameter_schema.to_dict()
    slurm_options = {
        "time": args.slurm_time,
        "mem": args.slurm_mem,
        "cpus_per_task": args.slurm_cpus_per_task,
        "gpus": args.slurm_gpus,
        "partition": args.slurm_partition,
        "account": args.slurm_account,
        "job_name_prefix": args.slurm_job_name_prefix,
        "submit_mode": args.slurm_submit_mode,
        "array_parallelism": args.slurm_array_parallelism,
        "extra_args": args.slurm_extra_args,
    }
    pair_source = str(rerun_manifest_file) if rerun_manifest_file else str(manifest.get("pairlist_file", ""))
    scheduler = str((hpc_profile or {}).get("scheduler") or "slurm").strip().lower()

    if scheduler == "condor":
        condor_options = {
            "cpus": args.condor_cpus,
            "memory": args.condor_memory,
            "disk": args.condor_disk,
            "gpus": args.condor_gpus,
            "submit_mode": args.condor_submit_mode,
            "parallelism": args.condor_parallelism,
            "requirements": args.condor_requirements,
            "job_name_prefix": args.condor_job_name_prefix,
        }
        try:
            summary = generate_condor_deployment(
                project_root=project_dir,
                engines=engines,
                pairlist_rows=rows,
                runtime_by_engine=runtime_by_engine,
                round_id=round_id,
                docking_mode=args.mode,
                condor_options=condor_options,
                hpc_profile=hpc_profile,
                remote_project_root=args.remote_project_dir or "",
                skip_completed=args.skip_completed,
                pair_source_file=pair_source,
                rerun_manifest_file=str(rerun_manifest_file) if rerun_manifest_file else "",
            )
        except ValueError as exc:
            print(f"❌ {exc}")
            return 1

        manifest.setdefault("engine_settings", {}).update(runtime_by_engine)
        manifest["deployment_root"] = str(summary.get("generation_stage_root") or summary["stage_root"])
        manifest["latest_pair_round"] = round_id
        if hpc_profile:
            manifest["hpc_profile"] = {
                "name": hpc_profile.get("name", ""),
                "source": hpc_profile.get("source", ""),
            }
        if rerun_manifest_file:
            manifest["latest_rerun_manifest_file"] = str(rerun_manifest_file)
        save_manifest(project_dir, manifest)

        print("✅ HTCondor deployment assets generated")
        print(f"   Stage root: {summary.get('generation_stage_root', summary['stage_root'])}")
        if args.remote_project_dir:
            print(f"   Execution project root: {args.remote_project_dir}")
            print(f"   Execution stage root: {summary['stage_root']}")
        print(f"   Round: {round_id}")
        print(f"   Mode: {args.mode}")
        print(f"   Pair source: {pair_source}")
        if hpc_profile:
            print(f"   HPC profile: {hpc_profile.get('name', '')} ({hpc_profile.get('source', '')})")
        for engine in engines:
            engine_summary = summary["engine_jobs"].get(engine, {})
            print(
                f"   {engine}: planned={engine_summary.get('planned', 0)} "
                f"skipped={engine_summary.get('skipped', 0)} "
                f"submit_mode={engine_summary.get('submit_mode', 'condor_array')}"
            )
            if engine_summary.get("submit_mode") == "condor_array":
                parallelism = engine_summary.get("array_parallelism", 0)
                cap = str(parallelism) if parallelism else "unlimited"
                print(f"      max_materialize={cap}")
        return 0

    try:
        summary = generate_slurm_deployment(
            project_root=project_dir,
            engines=engines,
            pairlist_rows=rows,
            runtime_by_engine=runtime_by_engine,
            round_id=round_id,
            docking_mode=args.mode,
            slurm_options=slurm_options,
            hpc_profile=hpc_profile,
            remote_project_root=args.remote_project_dir or "",
            skip_completed=args.skip_completed,
            pair_source_file=pair_source,
            rerun_manifest_file=str(rerun_manifest_file) if rerun_manifest_file else "",
        )
    except ValueError as exc:
        print(f"❌ {exc}")
        return 1

    manifest.setdefault("engine_settings", {}).update(runtime_by_engine)
    manifest["deployment_root"] = str(summary.get("generation_stage_root") or summary["stage_root"])
    manifest["latest_pair_round"] = round_id
    if hpc_profile:
        manifest["hpc_profile"] = {
            "name": hpc_profile.get("name", ""),
            "source": hpc_profile.get("source", ""),
        }
    if rerun_manifest_file:
        manifest["latest_rerun_manifest_file"] = str(rerun_manifest_file)
    save_manifest(project_dir, manifest)

    print("✅ Slurm deployment assets generated")
    print(f"   Stage root: {summary.get('generation_stage_root', summary['stage_root'])}")
    if args.remote_project_dir:
        print(f"   Execution project root: {args.remote_project_dir}")
        print(f"   Execution stage root: {summary['stage_root']}")
    print(f"   Round: {round_id}")
    print(f"   Mode: {args.mode}")
    print(f"   Pair source: {pair_source}")
    if hpc_profile:
        print(f"   HPC profile: {hpc_profile.get('name', '')} ({hpc_profile.get('source', '')})")
    for engine in engines:
        engine_summary = summary["engine_jobs"].get(engine, {})
        print(
            f"   {engine}: planned={engine_summary.get('planned', 0)} "
            f"skipped={engine_summary.get('skipped', 0)} "
            f"submit_mode={engine_summary.get('submit_mode', 'slurm_array')}"
        )
        if engine_summary.get("submit_mode") == "slurm_array":
            print(f"      array_parallelism={engine_summary.get('array_parallelism', 0)}")
    return 0


def sync_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Sync a prepared docking project from this machine to the HPC")
    parser.add_argument("--project-dir", required=True, help="Local docking project root")
    parser.add_argument("--hpc-profile", help="Builtin or local HPC deployment profile name")
    parser.add_argument("--hpc-profile-file", help="Explicit YAML/JSON file with local HPC deployment settings")
    parser.add_argument("--ssh-target", help="SSH target in the form user@host")
    parser.add_argument("--remote-project-dir", help="Project root on the HPC")
    parser.add_argument("--delete", action="store_true", help="Delete remote files that no longer exist locally")
    parser.add_argument("--dry-run", action="store_true", help="Print the sync commands without executing them")
    args = parser.parse_args(argv)

    project_dir = Path(args.project_dir).expanduser().resolve()
    try:
        hpc_profile = load_hpc_profile(
            project_dir,
            profile_name=args.hpc_profile or "",
            profile_file=args.hpc_profile_file or "",
        )
        ssh_target, remote_project_dir = resolve_remote_target(
            project_dir,
            hpc_profile,
            ssh_target=args.ssh_target or "",
            remote_project_dir=args.remote_project_dir or "",
        )
    except ValueError as exc:
        parser.error(str(exc))

    payload = sync_project_to_remote(
        local_project_root=project_dir,
        ssh_target=ssh_target,
        remote_project_dir=remote_project_dir,
        profile=hpc_profile,
        delete=args.delete,
        dry_run=args.dry_run,
    )
    if payload.get("status") == "failed":
        print("❌ Project sync failed")
        print(f"   SSH target: {ssh_target}")
        print(f"   Remote project root: {remote_project_dir}")
        if payload.get("failed_command"):
            print(f"   Failed step: {payload['failed_command']}")
        if payload.get("stderr"):
            print(f"   Error: {payload['stderr']}")
        return 1

    print("✅ Project sync ready" if args.dry_run else "✅ Project synced to HPC")
    print(f"   SSH target: {ssh_target}")
    print(f"   Remote project root: {remote_project_dir}")
    print(f"   mkdir command: {payload['commands']['mkdir']}")
    print(f"   rsync command: {payload['commands']['rsync']}")
    return 0


def submit_main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Submit a previously generated remote docking deployment on the HPC")
    parser.add_argument("--project-dir", required=True, help="Local docking project root")
    parser.add_argument("--engines", help="Comma-separated engine list (default: project manifest engines)")
    parser.add_argument("--round", help="Deployment round identifier (default: infer from manifest deployment_root)")
    parser.add_argument("--mode", choices=["screen", "exhaustive"], help="Deployment mode (default: infer from manifest deployment_root)")
    parser.add_argument("--hpc-profile", help="Builtin or local HPC deployment profile name")
    parser.add_argument("--hpc-profile-file", help="Explicit YAML/JSON file with local HPC deployment settings")
    parser.add_argument("--ssh-target", help="SSH target in the form user@host")
    parser.add_argument("--remote-project-dir", help="Project root on the HPC")
    parser.add_argument(
        "--allow-full-exhaustive",
        action="store_true",
        help="Allow exhaustive submission even when the local deployment manifest does not record a rerun manifest.",
    )
    parser.add_argument("--no-ligand-qc-gate", action="store_true", help="Disable ligand structural QC gate before remote submission")
    parser.add_argument("--no-receptor-qc-gate", action="store_true", help="Disable receptor structural QC gate before remote submission")
    parser.add_argument("--no-ligand-admet-filters", action="store_true", help="Disable ADMET-based ligand filters inside ligand QC gate")
    parser.add_argument("--ligand-qc-max-lipinski-violations", type=int, default=LIGAND_QC_MAX_LIPINSKI_DEFAULT, help="Maximum allowed Lipinski violations per ligand")
    parser.add_argument("--ligand-qc-max-molecular-weight", type=float, default=LIGAND_QC_MAX_MW_DEFAULT, help="Maximum allowed molecular weight (Da)")
    parser.add_argument("--ligand-qc-max-logp", type=float, default=LIGAND_QC_MAX_LOGP_DEFAULT, help="Maximum allowed LogP")
    parser.add_argument("--ligand-qc-max-tpsa", type=float, default=LIGAND_QC_MAX_TPSA_DEFAULT, help="Maximum allowed topological polar surface area")
    parser.add_argument("--ligand-qc-max-rotatable-bonds", type=int, default=LIGAND_QC_MAX_ROTATABLE_DEFAULT, help="Maximum allowed rotatable bonds")
    parser.add_argument("--ligand-qc-max-formal-charge-abs", type=int, default=LIGAND_QC_MAX_FORMAL_CHARGE_DEFAULT, help="Maximum allowed absolute formal charge")
    parser.add_argument("--ligand-qc-min-heavy-atom-count", type=int, default=LIGAND_QC_MIN_HEAVY_ATOM_DEFAULT, help="Minimum required heavy-atom count")
    parser.add_argument("--ligand-qc-allow-pains", action="store_true", help="Do not block ligands with PAINS alerts")
    parser.add_argument("--ligand-qc-allow-brenk", action="store_true", help="Do not block ligands with Brenk alerts")
    parser.add_argument("--ligand-qc-allow-reactive", action="store_true", help="Do not block ligands with reactive SMARTS alerts")
    parser.add_argument("--receptor-qc-min-atom-count", type=int, default=RECEPTOR_QC_MIN_ATOM_DEFAULT, help="Minimum receptor atom count for QC gate")
    parser.add_argument("--receptor-qc-min-heavy-atom-count", type=int, default=RECEPTOR_QC_MIN_HEAVY_ATOM_DEFAULT, help="Minimum receptor heavy-atom count for QC gate")
    parser.add_argument("--receptor-qc-min-chain-count", type=int, default=RECEPTOR_QC_MIN_CHAIN_DEFAULT, help="Minimum receptor chain count for QC gate")
    parser.add_argument("--receptor-qc-max-coordinate-span", type=float, default=RECEPTOR_QC_MAX_SPAN_DEFAULT, help="Maximum receptor coordinate span (A) for QC gate")
    parser.add_argument("--dry-run", action="store_true", help="Print the remote submit commands without executing them")
    args = parser.parse_args(argv)

    project_dir = Path(args.project_dir).expanduser().resolve()
    manifest = load_manifest(project_dir)
    try:
        hpc_profile = load_hpc_profile(
            project_dir,
            profile_name=args.hpc_profile or "",
            profile_file=args.hpc_profile_file or "",
        )
        ssh_target, remote_project_dir = resolve_remote_target(
            project_dir,
            hpc_profile,
            ssh_target=args.ssh_target or "",
            remote_project_dir=args.remote_project_dir or "",
        )
    except ValueError as exc:
        parser.error(str(exc))

    inferred_round, inferred_mode = infer_round_mode_from_deployment_root(str(manifest.get("deployment_root", "") or ""))
    round_id = args.round or inferred_round
    mode = args.mode or inferred_mode or "screen"
    if not round_id:
        parser.error("No deployment round could be inferred. Provide --round or generate a deployment bundle first.")
    _validate_exhaustive_submit(parser, args, project_dir, round_id, mode)
    rerun_manifest = Path(str(manifest.get("latest_rerun_manifest_file", "") or "")).expanduser().resolve() if mode == "exhaustive" and str(manifest.get("latest_rerun_manifest_file", "") or "").strip() else None
    if rerun_manifest and rerun_manifest.exists():
        rows, _ = _load_rerun_manifest_rows(rerun_manifest)
    else:
        rows = _load_pairlist_rows(project_dir)
    if not bool(args.no_ligand_qc_gate):
        ligand_audit = _audit_project_ligands_for_rows(
            project_dir,
            rows,
            enable_admet_filters=not bool(args.no_ligand_admet_filters),
            max_lipinski_violations=int(args.ligand_qc_max_lipinski_violations),
            max_molecular_weight=float(args.ligand_qc_max_molecular_weight),
            max_logp=float(args.ligand_qc_max_logp),
            max_tpsa=float(args.ligand_qc_max_tpsa),
            max_rotatable_bonds=int(args.ligand_qc_max_rotatable_bonds),
            max_formal_charge_abs=int(args.ligand_qc_max_formal_charge_abs),
            min_heavy_atom_count=int(args.ligand_qc_min_heavy_atom_count),
            block_pains=not bool(args.ligand_qc_allow_pains),
            block_brenk=not bool(args.ligand_qc_allow_brenk),
            block_reactive=not bool(args.ligand_qc_allow_reactive),
        )
        if ligand_audit.get("issue_count"):
            _print_ligand_validation_failure(ligand_audit)
            return 1
    if not bool(args.no_receptor_qc_gate):
        receptor_audit = _audit_project_receptors_for_rows(
            project_dir,
            rows,
            min_atom_count=int(args.receptor_qc_min_atom_count),
            min_heavy_atom_count=int(args.receptor_qc_min_heavy_atom_count),
            min_chain_count=int(args.receptor_qc_min_chain_count),
            max_coordinate_span=float(args.receptor_qc_max_coordinate_span),
        )
        if receptor_audit.get("issue_count"):
            _print_receptor_validation_failure(receptor_audit)
            return 1
    engines = normalize_engines([args.engines]) if args.engines else normalize_engines(manifest.get("engines", []))

    payload = submit_remote_deployment(
        ssh_target=ssh_target,
        remote_project_dir=remote_project_dir,
        round_id=round_id,
        mode=mode,
        engines=engines,
        dry_run=args.dry_run,
    )
    if payload.get("status") == "failed":
        print("❌ Remote submit failed")
        print(f"   SSH target: {ssh_target}")
        print(f"   Remote project root: {remote_project_dir}")
        for result in payload.get("results", []):
            print(
                f"   {result['engine']}: status={result['status']} returncode={result['returncode']}"
            )
            if result.get("job_id"):
                print(f"      job_id: {result['job_id']}")
            elif result.get("stdout"):
                print(f"      stdout: {result['stdout'].strip()}")
            if result.get("stderr"):
                print(f"      stderr: {result['stderr']}")
        return 1

    print("✅ Remote submit ready" if args.dry_run else "✅ Remote submit complete")
    print(f"   SSH target: {ssh_target}")
    print(f"   Remote project root: {remote_project_dir}")
    print(f"   Round: {round_id}")
    print(f"   Mode: {mode}")
    if args.dry_run:
        for command in payload.get("commands", []):
            print(f"   submit command: {command}")
        return 0
    for result in payload.get("results", []):
        line = f"   {result['engine']}: status={result['status']}"
        if result.get("job_id"):
            line += f" job_id={result['job_id']}"
        elif result.get("job_ids"):
            line += f" job_ids={','.join(result['job_ids'])}"
        print(line)
    return 0
