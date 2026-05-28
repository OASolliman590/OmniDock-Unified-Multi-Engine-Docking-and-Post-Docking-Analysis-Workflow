from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict, List

from .models import PairlistRow
from .project_layout import shared_ligands_dir, shared_receptors_dir
from .preparation.ligand_quality import audit_project_ligands
from .preparation.receptor_quality import audit_project_receptors


def _exists_nonempty_dir(path: Path) -> bool:
    if not path.exists() or not path.is_dir():
        return False
    return any(path.iterdir())


def _binary_resolved(binary: str) -> bool:
    candidate = str(binary or "").strip()
    if not candidate:
        return False
    as_path = Path(candidate).expanduser()
    if as_path.is_absolute():
        return as_path.exists()
    return shutil.which(candidate) is not None


def run_docking_preflight(
    *,
    project_dir: Path,
    engines: List[str],
    pairlist_rows: List[PairlistRow],
    runtime_by_engine: Dict[str, Dict[str, object]],
    execution_workdir: str = "",
    prerequisites_dir: str = "",
    required_files_dir: str = "",
    enable_ligand_qc: bool = False,
    enable_receptor_qc: bool = False,
    ligand_qc_report_path: str = "",
    receptor_qc_report_path: str = "",
    ligand_admet_filters_enabled: bool = True,
    ligand_max_lipinski_violations: int = 1,
    ligand_max_molecular_weight: float = 650.0,
    ligand_max_logp: float = 6.0,
    ligand_max_tpsa: float = 180.0,
    ligand_max_rotatable_bonds: int = 15,
    ligand_max_formal_charge_abs: int = 2,
    ligand_min_heavy_atom_count: int = 6,
    ligand_block_pains: bool = True,
    ligand_block_brenk: bool = True,
    ligand_block_reactive: bool = True,
    receptor_min_atom_count: int = 100,
    receptor_min_heavy_atom_count: int = 60,
    receptor_min_chain_count: int = 1,
    receptor_max_coordinate_span: float = 500.0,
    dry_run: bool = False,
) -> Dict[str, object]:
    root = Path(project_dir).expanduser().resolve()
    shared_receptors = shared_receptors_dir(root)
    shared_ligands = shared_ligands_dir(root)

    errors: List[str] = []
    warnings: List[str] = []
    details: Dict[str, object] = {}

    if not pairlist_rows:
        errors.append("No pairlist rows were found; docking cannot start.")
    details["pair_count"] = len(pairlist_rows)

    if execution_workdir:
        workdir = Path(execution_workdir).expanduser().resolve()
        if not workdir.exists() or not workdir.is_dir():
            errors.append(f"Execution workdir does not exist or is not a directory: {workdir}")
        else:
            details["execution_workdir"] = str(workdir)

    if prerequisites_dir:
        prereq = Path(prerequisites_dir).expanduser().resolve()
        if not prereq.exists() or not prereq.is_dir():
            errors.append(f"Prerequisites directory does not exist or is not a directory: {prereq}")
        else:
            details["prerequisites_dir"] = str(prereq)

    if required_files_dir:
        required = Path(required_files_dir).expanduser().resolve()
        if not required.exists() or not required.is_dir():
            errors.append(f"Required-files directory does not exist or is not a directory: {required}")
        elif not _exists_nonempty_dir(required):
            errors.append(f"Required-files directory is empty: {required}")
        else:
            details["required_files_dir"] = str(required)

    missing_receptors: List[str] = []
    missing_ligands: List[str] = []
    ad4_receptor_dir_candidates = [
        root / "4-Docking" / "receptors_engine_specific" / "autodock4",
        root / "receptors_engine_specific" / "autodock4",
    ]
    ad4_ligand_dir_candidates = [
        root / "4-Docking" / "ligands_engine_specific" / "autodock4",
        root / "ligands_engine_specific" / "autodock4",
    ]

    for row in pairlist_rows:
        receptor_name = str(row.receptor)
        ligand_name = str(row.ligand)

        receptor_ok = (shared_receptors / receptor_name).exists()
        ligand_ok = (shared_ligands / ligand_name).exists()

        if "autodock4" in engines:
            receptor_ok = receptor_ok or any((candidate / receptor_name).exists() for candidate in ad4_receptor_dir_candidates)
            ligand_ok = ligand_ok or any((candidate / ligand_name).exists() for candidate in ad4_ligand_dir_candidates)

        if not receptor_ok:
            missing_receptors.append(receptor_name)
        if not ligand_ok:
            missing_ligands.append(ligand_name)

    if missing_receptors:
        unique_missing = sorted(set(missing_receptors))
        preview = ", ".join(unique_missing[:10])
        errors.append(f"Missing receptor assets ({len(unique_missing)}): {preview}")
    if missing_ligands:
        unique_missing = sorted(set(missing_ligands))
        preview = ", ".join(unique_missing[:10])
        errors.append(f"Missing ligand assets ({len(unique_missing)}): {preview}")

    if enable_ligand_qc and pairlist_rows:
        try:
            admet_thresholds = {
                "enable_admet_filters": bool(ligand_admet_filters_enabled),
                "max_lipinski_violations": int(ligand_max_lipinski_violations),
                "max_molecular_weight": float(ligand_max_molecular_weight),
                "max_logp": float(ligand_max_logp),
                "max_tpsa": float(ligand_max_tpsa),
                "max_rotatable_bonds": int(ligand_max_rotatable_bonds),
                "max_formal_charge_abs": int(ligand_max_formal_charge_abs),
                "min_heavy_atom_count": int(ligand_min_heavy_atom_count),
                "block_pains": bool(ligand_block_pains),
                "block_brenk": bool(ligand_block_brenk),
                "block_reactive": bool(ligand_block_reactive),
            }
            ligand_report = audit_project_ligands(
                project_root=root,
                rows=pairlist_rows,
                enable_admet_filters=bool(ligand_admet_filters_enabled),
                admet_thresholds=admet_thresholds,
                report_path=Path(ligand_qc_report_path).expanduser() if ligand_qc_report_path else None,
            )
            details["ligand_qc"] = ligand_report
            warning_count = int(ligand_report.get("warning_count", 0) or 0)
            if warning_count > 0:
                warnings.append(
                    f"Ligand ADMET QC produced {warning_count} warning(s); see ligand QC report for details."
                )
            if int(ligand_report.get("issue_count", 0) or 0) > 0:
                report_hint = str(ligand_report.get("report_path", "") or "")
                hint = f" (report: {report_hint})" if report_hint else ""
                errors.append(
                    f"Ligand QC gate failed with {int(ligand_report.get('issue_count', 0) or 0)} issue(s){hint}"
                )
        except Exception as exc:
            warnings.append(f"Ligand QC gate could not complete: {exc}")

    if enable_receptor_qc and pairlist_rows:
        try:
            receptor_report = audit_project_receptors(
                project_root=root,
                rows=pairlist_rows,
                min_atom_count=int(receptor_min_atom_count),
                min_heavy_atom_count=int(receptor_min_heavy_atom_count),
                min_chain_count=int(receptor_min_chain_count),
                max_coordinate_span=float(receptor_max_coordinate_span),
                report_path=Path(receptor_qc_report_path).expanduser() if receptor_qc_report_path else None,
            )
            details["receptor_qc"] = receptor_report
            if int(receptor_report.get("issue_count", 0) or 0) > 0:
                report_hint = str(receptor_report.get("report_path", "") or "")
                hint = f" (report: {report_hint})" if report_hint else ""
                errors.append(
                    f"Receptor QC gate failed with {int(receptor_report.get('issue_count', 0) or 0)} issue(s){hint}"
                )
        except Exception as exc:
            warnings.append(f"Receptor QC gate could not complete: {exc}")

    if not dry_run:
        for engine in engines:
            runtime = runtime_by_engine.get(engine, {})
            if engine == "gnina":
                image = str(runtime.get("image") or "").strip()
                if image:
                    if not _binary_resolved("apptainer"):
                        warnings.append("GNINA image configured but `apptainer` is not available in PATH.")
                else:
                    binary = str(runtime.get("binary") or "gnina")
                    if not _binary_resolved(binary):
                        errors.append(f"GNINA binary not found: {binary}")
            elif engine == "vina":
                binary = str(runtime.get("binary") or "vina")
                conda_env = str(runtime.get("conda_env") or "").strip()
                if not conda_env and not _binary_resolved(binary):
                    errors.append(f"Vina binary not found: {binary}")
            elif engine == "smina":
                binary = str(runtime.get("binary") or "smina")
                conda_env = str(runtime.get("conda_env") or "").strip()
                if not conda_env and not _binary_resolved(binary):
                    errors.append(f"Smina binary not found: {binary}")
            elif engine == "autodock4":
                autodock_binary = str(runtime.get("binary") or runtime.get("autodock_binary") or "autodock4")
                autogrid_binary = str(runtime.get("autogrid_binary") or "autogrid4")
                if not _binary_resolved(autodock_binary):
                    errors.append(f"AutoDock4 binary not found: {autodock_binary}")
                if not _binary_resolved(autogrid_binary):
                    errors.append(f"AutoGrid4 binary not found: {autogrid_binary}")

    return {
        "ok": not errors,
        "errors": errors,
        "warnings": warnings,
        "details": details,
        "engines": list(engines),
        "pair_count": len(pairlist_rows),
        "dry_run": bool(dry_run),
    }
