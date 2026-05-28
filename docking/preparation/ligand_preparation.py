from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from docking.models import (
    DEFAULT_PREPARATION_PH,
    LIGAND_PREPARATION_PROFILES,
    MAX_PREPARATION_PH,
    MIN_PREPARATION_PH,
    normalize_engine_names,
    normalize_ligand_preparation_profile,
    resolve_effective_ligand_preparation_profile,
    validate_preparation_ph,
    validate_ligand_preparation_profile,
)


def _resolve_profile_alias(raw: str) -> str:
    normalized = normalize_ligand_preparation_profile(raw)
    return normalized if normalized in LIGAND_PREPARATION_PROFILES else "engine_aware_full"


def _run_command(
    command: list[str],
    env: Optional[Dict[str, str]] = None,
    cwd: Optional[Path] = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        env=env,
        cwd=str(cwd) if cwd else None,
    )


def _selected_profile() -> str:
    raw = str(
        os.environ.get("PDBWIZARD_LIGAND_PREP_PROFILE")
        or os.environ.get("PDBWIZARD_LIGAND_PREP_BACKEND")
        or "engine_aware_full"
    )
    return _resolve_profile_alias(raw)


def _selected_ph() -> float:
    raw = str(os.environ.get("PDBWIZARD_LIGAND_PREP_PH", str(DEFAULT_PREPARATION_PH)) or "").strip()
    ok, normalized, _ = validate_preparation_ph(raw)
    if not ok:
        return float(DEFAULT_PREPARATION_PH)
    return float(normalized)


def _selected_engines() -> Tuple[str, ...]:
    raw = str(os.environ.get("PDBWIZARD_SELECTED_ENGINES", "") or "").strip()
    if not raw:
        return ()
    engines = normalize_engine_names(raw.split(","))
    return tuple(engines)


def _resolve_prepare_ligand4_script() -> Optional[str]:
    repo_root = Path(__file__).resolve().parents[2]
    candidates = [
        os.environ.get("AUTODOCKTOOLS_PREPARE_LIGAND4"),
        os.environ.get("ADT_PREPARE_LIGAND4"),
        shutil.which("prepare_ligand4.py"),
        str(repo_root / ".workflow" / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"),
        str(repo_root / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"),
    ]
    for candidate in candidates:
        if not candidate:
            continue
        script = Path(candidate).expanduser()
        if script.exists():
            return str(script.resolve())
        if candidate.endswith(".py"):
            continue
        possible_py = Path(f"{candidate}.py").expanduser()
        if possible_py.exists():
            return str(possible_py.resolve())
    return None


def _is_zero_molecule_obabel_result(result: Optional[subprocess.CompletedProcess[str]]) -> bool:
    if result is None:
        return False
    stderr = str(result.stderr or "").lower()
    stdout = str(result.stdout or "").lower()
    return "0 molecules converted" in stderr or "0 molecules converted" in stdout


def _file_has_meaningful_content(path: Path) -> bool:
    file_path = Path(path)
    if not file_path.exists() or file_path.stat().st_size == 0:
        return False
    try:
        content = file_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return file_path.stat().st_size > 0
    return bool(content.strip())


def _mol2_has_atoms(path: Path) -> bool:
    file_path = Path(path)
    if not _file_has_meaningful_content(file_path):
        return False
    try:
        content = file_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return False
    return "@<TRIPOS>ATOM" in content


def _prepare_with_autodocktools(input_sdf: Path, output_pdbqt: Path) -> Dict[str, object]:
    script = _resolve_prepare_ligand4_script()
    if not script:
        raise FileNotFoundError(
            "AutoDockTools prepare_ligand4.py was not found. "
            "Set AUTODOCKTOOLS_PREPARE_LIGAND4 (or ADT_PREPARE_LIGAND4) to enable this fallback."
        )

    python_exe = os.environ.get("AUTODOCKTOOLS_PYTHON") or os.environ.get("ADT_PYTHON") or sys.executable
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    script_path = Path(script).expanduser().resolve()
    adt_repo_root = script_path.parent.parent.parent
    env = dict(os.environ)
    existing_pythonpath = str(env.get("PYTHONPATH", "") or "").strip()
    env["PYTHONPATH"] = str(adt_repo_root) if not existing_pythonpath else f"{adt_repo_root}:{existing_pythonpath}"
    input_file = Path(input_sdf).expanduser().resolve()
    output_file = Path(output_pdbqt).expanduser().resolve()
    work_dir = input_file.parent
    mol2_input = work_dir / f"{input_file.stem}_adt_input.mol2"
    local_output = work_dir / f"{input_file.stem}_adt.pdbqt"

    ligand_arg = input_file.name
    try:
        mol2_result = _run_command(["obabel", str(input_file), "-O", str(mol2_input)])
        if not _is_zero_molecule_obabel_result(mol2_result) and _mol2_has_atoms(mol2_input):
            ligand_arg = mol2_input.name
    except Exception:
        # Fall back to the original SDF when conversion is unavailable.
        ligand_arg = input_file.name

    command = [
        python_exe,
        str(script_path),
        "-l",
        ligand_arg,
        "-o",
        local_output.name,
        "-p",
        "Pd",
    ]
    result = _run_command(command, env=env, cwd=work_dir)
    if not local_output.exists():
        raise FileNotFoundError(f"AutoDockTools did not produce expected output file: {local_output}")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if local_output != output_file:
        shutil.copy2(local_output, output_file)
    return {
        "preparation_method": "autodocktools",
        "autodocktools_script": script,
        "autodocktools_python": python_exe,
        "autodocktools_command": command,
        "autodocktools_workdir": str(work_dir),
        "autodocktools_stdout": (result.stdout or "").strip(),
        "autodocktools_stderr": (result.stderr or "").strip(),
    }


def _load_rdkit() -> Tuple[Optional[Any], Optional[Any]]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return None, None
    return Chem, AllChem


def _load_rdkit_molecule(input_path: Path) -> Any:
    Chem, _ = _load_rdkit()
    if Chem is None:
        raise ImportError("RDKit is not available")

    suffix = input_path.suffix.lower()
    if suffix == ".sdf":
        supplier = Chem.SDMolSupplier(str(input_path), removeHs=False)
        for mol in supplier:
            if mol is not None:
                return mol
        return None
    if suffix == ".mol":
        return Chem.MolFromMolFile(str(input_path), removeHs=False)
    if suffix == ".mol2":
        return Chem.MolFromMol2File(str(input_path), removeHs=False)
    if suffix == ".pdb":
        return Chem.MolFromPDBFile(str(input_path), removeHs=False)
    raise ValueError(f"Unsupported ligand format for RDKit normalization: {input_path.suffix}")


def _has_3d_coordinates(mol: Any) -> bool:
    if mol is None or mol.GetNumConformers() == 0:
        return False
    conformer = mol.GetConformer()
    if conformer.Is3D():
        return True
    positions = conformer.GetPositions()
    if len(positions) == 0:
        return False
    return max(abs(float(position[2])) for position in positions) > 1e-3


def _embed_and_optimize_3d(mol: Any, all_chem: Any) -> str:
    params = all_chem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    status = all_chem.EmbedMolecule(mol, params)
    if status != 0:
        params.useRandomCoords = True
        status = all_chem.EmbedMolecule(mol, params)
    if status != 0:
        raise ValueError("RDKit could not embed a 3D conformer")

    try:
        if all_chem.MMFFHasAllMoleculeParams(mol):
            all_chem.MMFFOptimizeMolecule(mol)
            return "MMFF94"
        all_chem.UFFOptimizeMolecule(mol)
        return "UFF"
    except Exception:
        return ""


def _normalize_with_rdkit(input_path: Path, output_sdf: Path) -> Dict[str, object]:
    Chem, all_chem = _load_rdkit()
    if Chem is None or all_chem is None:
        raise ImportError("RDKit is not available")

    molecule = _load_rdkit_molecule(input_path)
    if molecule is None:
        raise ValueError(f"RDKit could not parse ligand file: {input_path}")

    Chem.SanitizeMol(molecule)
    had_3d = _has_3d_coordinates(molecule)
    molecule = Chem.AddHs(molecule, addCoords=had_3d)

    optimized_forcefield = ""
    if not had_3d:
        optimized_forcefield = _embed_and_optimize_3d(molecule, all_chem)

    molecule.SetProp("_Name", molecule.GetProp("_Name") if molecule.HasProp("_Name") else input_path.stem)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_sdf))
    try:
        writer.write(molecule)
    finally:
        writer.close()

    return {
        "normalization_backend": "rdkit",
        "source_file": str(input_path),
        "normalized_sdf": str(output_sdf),
        "had_3d_input": had_3d,
        "generated_3d": not had_3d,
        "optimized_forcefield": optimized_forcefield,
    }


def _normalize_with_openbabel(
    input_path: Path,
    output_sdf: Path,
    *,
    protonation_ph: float = 7.4,
    partial_charge_model: str = "gasteiger",
) -> Dict[str, object]:
    if not shutil.which("obabel"):
        raise FileNotFoundError("Open Babel (obabel) is required for ligand normalization fallback")

    with tempfile.TemporaryDirectory(prefix=f"ligand_normalize_{input_path.stem}_") as tmp_dir:
        temp_dir = Path(tmp_dir)
        temp_source = temp_dir / f"{input_path.stem}_source.sdf"
        temp_3d = temp_dir / f"{input_path.stem}_3d.sdf"
        temp_h = temp_dir / f"{input_path.stem}_h.sdf"
        temp_charged = temp_dir / f"{input_path.stem}_charged.sdf"
        temp_min = temp_dir / f"{input_path.stem}_min.sdf"

        source_result = _run_command(["obabel", str(input_path), "-O", str(temp_source)])
        if _is_zero_molecule_obabel_result(source_result) or not _file_has_meaningful_content(temp_source):
            raise ValueError(f"OpenBabel could not parse ligand source into SDF: {input_path}")
        gen3d_result = _run_command(["obabel", str(temp_source), "-O", str(temp_3d), "--gen3d"])
        if _is_zero_molecule_obabel_result(gen3d_result) or not _file_has_meaningful_content(temp_3d):
            raise ValueError(f"OpenBabel failed to generate 3D conformer for: {input_path}")
        protonation_applied = False
        try:
            addh_result = _run_command(
                [
                    "obabel",
                    str(temp_3d),
                    "-O",
                    str(temp_h),
                    "-h",
                    "-p",
                    f"{protonation_ph:.2f}",
                ]
            )
            protonation_applied = True
        except subprocess.CalledProcessError:
            addh_result = _run_command(["obabel", str(temp_3d), "-O", str(temp_h), "-h"])
        if _is_zero_molecule_obabel_result(addh_result) or not _file_has_meaningful_content(temp_h):
            raise ValueError(f"OpenBabel failed to add hydrogens for: {input_path}")

        final_source = temp_h
        partial_charge_applied = False
        if partial_charge_model:
            try:
                charge_result = _run_command(
                    [
                        "obabel",
                        str(temp_h),
                        "-O",
                        str(temp_charged),
                        "--partialcharge",
                        partial_charge_model,
                    ]
                )
                if not _is_zero_molecule_obabel_result(charge_result) and _file_has_meaningful_content(temp_charged):
                    final_source = temp_charged
                    partial_charge_applied = True
            except subprocess.CalledProcessError:
                partial_charge_applied = False

        optimized_forcefield = ""
        for forcefield in ("MMFF94", "UFF"):
            try:
                min_result = _run_command(
                    [
                        "obabel",
                        str(final_source),
                        "-O",
                        str(temp_min),
                        "--minimize",
                        "--steps",
                        "250",
                        "--ff",
                        forcefield,
                    ]
                )
                if not _is_zero_molecule_obabel_result(min_result) and _file_has_meaningful_content(temp_min):
                    final_source = temp_min
                    optimized_forcefield = forcefield
                    break
            except subprocess.CalledProcessError:
                continue

        output_sdf.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(final_source, output_sdf)
        if not _file_has_meaningful_content(output_sdf):
            raise ValueError(f"OpenBabel normalization produced an empty SDF for: {input_path}")

    return {
        "normalization_backend": "openbabel",
        "source_file": str(input_path),
        "normalized_sdf": str(output_sdf),
        "had_3d_input": False,
        "generated_3d": True,
        "protonation_ph": protonation_ph,
        "protonation_applied": protonation_applied,
        "partial_charge_model": partial_charge_model,
        "partial_charge_applied": partial_charge_applied,
        "optimized_forcefield": optimized_forcefield,
    }


def normalize_ligand_to_sdf(
    input_path: Path,
    output_sdf: Path,
    *,
    protonation_ph: float = 7.4,
) -> Dict[str, object]:
    source = Path(input_path).expanduser().resolve()
    destination = Path(output_sdf).expanduser().resolve()
    errors: list[str] = []

    # Enforce a single OpenBabel-first 3D normalization path across engines.
    try:
        return _normalize_with_openbabel(source, destination, protonation_ph=protonation_ph)
    except Exception as exc:
        errors.append(f"openbabel:{exc}")

    payload = _normalize_with_rdkit(source, destination)
    payload["fallback_reasons"] = errors
    payload["protonation_ph"] = protonation_ph
    return payload


def _prepare_with_meeko(input_sdf: Path, output_pdbqt: Path) -> None:
    _run_command(["mk_prepare_ligand.py", "-i", str(input_sdf), "-o", str(output_pdbqt)])


def _prepare_with_openbabel_pdbqt(input_sdf: Path, output_pdbqt: Path) -> None:
    _run_command(
        [
            "obabel",
            str(input_sdf),
            "-O",
            str(output_pdbqt),
            "-h",
            "--partialcharge",
            "gasteiger",
        ]
    )


def _effective_profile(profile: str, selected_engines: Tuple[str, ...]) -> str:
    return resolve_effective_ligand_preparation_profile(profile, list(selected_engines))


def prepare_ligand_for_vina_family(input_path: Path, output_pdbqt: Path) -> Dict[str, object]:
    source = Path(input_path).expanduser().resolve()
    destination = Path(output_pdbqt).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)

    selected_profile = _selected_profile()
    selected_engines = _selected_engines()
    compatibility = validate_ligand_preparation_profile(selected_profile, list(selected_engines))
    if not compatibility.is_valid:
        raise RuntimeError("; ".join(compatibility.errors))
    effective_profile = _effective_profile(selected_profile, selected_engines)
    protonation_ph = _selected_ph()

    with tempfile.TemporaryDirectory(prefix=f"ligand_prepare_{source.stem}_") as tmp_dir:
        normalized_sdf = Path(tmp_dir) / f"{source.stem}_normalized.sdf"
        meeko_error = ""
        obabel_error = ""
        autodocktools_error = ""
        prep_method = ""
        autodocktools_summary: Optional[Dict[str, object]] = None
        normalization: Dict[str, object] = {}

        if effective_profile == "meeko_only":
            try:
                _prepare_with_meeko(source, destination)
                prep_method = "meeko_direct"
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                meeko_error = str(exc)
                raise RuntimeError("Ligand preparation failed in Meeko-only mode.") from exc
        elif effective_profile == "autodocktools_only":
            try:
                autodocktools_summary = _prepare_with_autodocktools(source, destination)
                prep_method = "autodocktools_direct"
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                autodocktools_error = str(exc)
                raise RuntimeError("Ligand preparation failed in AutoDockTools-only mode.") from exc
        else:
            normalization = normalize_ligand_to_sdf(source, normalized_sdf, protonation_ph=protonation_ph)
            if effective_profile == "openbabel_only":
                try:
                    _prepare_with_openbabel_pdbqt(normalized_sdf, destination)
                    prep_method = "openbabel_pdbqt"
                except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                    obabel_error = str(exc)
                    raise RuntimeError("Ligand preparation failed in OpenBabel-only mode.") from exc
            elif effective_profile == "openbabel_meeko":
                try:
                    _prepare_with_meeko(normalized_sdf, destination)
                    prep_method = "openbabel_then_meeko"
                except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                    meeko_error = str(exc)
                    raise RuntimeError("Ligand preparation failed in OpenBabel->Meeko mode.") from exc
            elif effective_profile == "openbabel_autodocktools":
                try:
                    autodocktools_summary = _prepare_with_autodocktools(normalized_sdf, destination)
                    prep_method = "openbabel_then_autodocktools"
                except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                    autodocktools_error = str(exc)
                    raise RuntimeError("Ligand preparation failed in OpenBabel->AutoDockTools mode.") from exc
            else:
                # openbabel_meeko_autodock
                try:
                    _prepare_with_meeko(normalized_sdf, destination)
                    prep_method = "openbabel_then_meeko"
                except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                    meeko_error = str(exc)
                    try:
                        autodocktools_summary = _prepare_with_autodocktools(normalized_sdf, destination)
                        prep_method = "openbabel_then_meeko_fallback_autodocktools"
                    except (subprocess.CalledProcessError, FileNotFoundError) as adt_exc:
                        autodocktools_error = str(adt_exc)
                        raise RuntimeError(
                            "Ligand preparation failed in OpenBabel->Meeko->AutoDockTools mode."
                        ) from exc

    summary = {
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_file": str(source),
        "output_file": str(destination),
        "preparation_method": prep_method,
        "requested_backend": selected_profile,
        "requested_profile": selected_profile,
        "effective_profile": effective_profile,
        "selected_engines": list(selected_engines),
        "protonation_ph": protonation_ph,
        "normalization": normalization,
        "compatibility": compatibility.to_dict(),
    }
    if meeko_error:
        summary["meeko_error"] = meeko_error
    if obabel_error:
        summary["obabel_error"] = obabel_error
    if autodocktools_error:
        summary["autodocktools_error"] = autodocktools_error
    if autodocktools_summary:
        summary["autodocktools"] = autodocktools_summary
    report_dir = str(os.environ.get("PDBWIZARD_LIGAND_PREP_REPORT_DIR", "") or "").strip()
    if report_dir:
        report_path = write_ligand_preparation_step_report(summary, report_dir=Path(report_dir))
        if report_path:
            summary["step_report_file"] = str(report_path)
    validation = validate_ligand_preparation_output_contract(summary)
    summary["output_contract_validation"] = validation
    return summary


def _parse_float(raw: object) -> Optional[float]:
    try:
        return float(raw)
    except Exception:
        return None


def _count_line_token(lines: list[str], token: str) -> int:
    return sum(1 for line in lines if line.strip() == token)


def _validate_prepared_pdbqt_contract(output_file: Path) -> Dict[str, object]:
    errors: list[str] = []
    warnings: list[str] = []
    if not output_file.exists():
        errors.append("prepared output file is missing")
        return {"is_valid": False, "errors": errors, "warnings": warnings}
    text = output_file.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()
    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    root_count = _count_line_token(lines, "ROOT")
    torsdof_count = sum(1 for line in lines if line.startswith("TORSDOF"))
    if not atom_lines:
        errors.append("prepared PDBQT has no ATOM/HETATM records")
    if root_count != 1:
        errors.append(f"prepared PDBQT expected exactly one ROOT block, found {root_count}")
    if torsdof_count != 1:
        errors.append(f"prepared PDBQT expected exactly one TORSDOF record, found {torsdof_count}")
    return {"is_valid": not errors, "errors": errors, "warnings": warnings}


def validate_ligand_preparation_output_contract(summary: Dict[str, object]) -> Dict[str, object]:
    """
    Validate mode-specific output contracts for one prepared ligand summary payload.
    """
    requested_profile = normalize_ligand_preparation_profile(str(summary.get("requested_profile") or summary.get("requested_backend") or "engine_aware_full"))
    selected_engines = [str(engine).strip().lower() for engine in (summary.get("selected_engines") or []) if str(engine).strip()]
    compatibility = validate_ligand_preparation_profile(requested_profile, selected_engines)
    effective_profile = str(summary.get("effective_profile") or compatibility.effective_profile or "")
    output_file = Path(str(summary.get("output_file") or "")).expanduser()
    errors: list[str] = []
    warnings: list[str] = []

    if not compatibility.is_valid:
        errors.extend(compatibility.errors)
    warnings.extend(compatibility.warnings)

    expected_methods = {
        "meeko_only": {"meeko_direct"},
        "autodocktools_only": {"autodocktools_direct"},
        "openbabel_only": {"openbabel_pdbqt"},
        "openbabel_meeko": {"openbabel_then_meeko"},
        "openbabel_autodocktools": {"openbabel_then_autodocktools"},
        "openbabel_meeko_autodock": {"openbabel_then_meeko", "openbabel_then_meeko_fallback_autodocktools"},
    }
    preparation_method = str(summary.get("preparation_method") or "")
    allowed_methods = expected_methods.get(effective_profile)
    if allowed_methods and preparation_method not in allowed_methods:
        errors.append(
            f"preparation_method='{preparation_method}' does not match effective_profile='{effective_profile}' "
            f"(expected one of: {', '.join(sorted(allowed_methods))})"
        )

    ph_value = _parse_float(summary.get("protonation_ph"))
    ph_ok, normalized_ph, ph_error = validate_preparation_ph(ph_value if ph_value is not None else DEFAULT_PREPARATION_PH)
    if not ph_ok:
        errors.append(ph_error)
    elif ph_value is not None and abs(ph_value - normalized_ph) > 1e-6:
        warnings.append("protonation_ph was normalized to a valid range")

    normalization = summary.get("normalization") or {}
    if effective_profile.startswith("openbabel_") or effective_profile == "openbabel_only":
        if not isinstance(normalization, dict) or not normalization:
            errors.append("openbabel-based profiles require normalization details in summary")
        else:
            backend = str(normalization.get("normalization_backend") or "").strip().lower()
            if backend not in {"openbabel", "rdkit"}:
                errors.append(f"unexpected normalization backend: {backend or 'missing'}")
            if backend == "rdkit":
                warnings.append("normalization used RDKit fallback instead of Open Babel")

    output_contract = _validate_prepared_pdbqt_contract(output_file)
    if not output_contract.get("is_valid", False):
        errors.extend([str(issue) for issue in output_contract.get("errors", [])])
    warnings.extend([str(issue) for issue in output_contract.get("warnings", [])])

    return {
        "requested_profile": requested_profile,
        "effective_profile": effective_profile,
        "selected_engines": selected_engines,
        "output_file": str(output_file),
        "protonation_ph": float(ph_value) if ph_value is not None else float(DEFAULT_PREPARATION_PH),
        "ph_bounds": {
            "min": float(MIN_PREPARATION_PH),
            "max": float(MAX_PREPARATION_PH),
            "default": float(DEFAULT_PREPARATION_PH),
        },
        "is_valid": not errors,
        "errors": errors,
        "warnings": warnings,
    }


def write_ligand_preparation_step_report(
    summary: Dict[str, object],
    *,
    report_dir: Optional[Path] = None,
    report_file: Optional[Path] = None,
) -> Optional[Path]:
    """
    Persist a per-ligand preparation step report as JSON.
    """
    target_path: Optional[Path] = None
    if report_file is not None:
        target_path = Path(report_file).expanduser().resolve()
    elif report_dir is not None:
        directory = Path(report_dir).expanduser().resolve()
        input_file = Path(str(summary.get("input_file") or "ligand")).name
        stem = Path(input_file).stem or "ligand"
        safe_stem = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in stem).strip("_") or "ligand"
        target_path = directory / f"{safe_stem}.json"
    if target_path is None:
        return None
    target_path.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(summary)
    payload["output_contract_validation"] = payload.get("output_contract_validation") or validate_ligand_preparation_output_contract(payload)
    target_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return target_path


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize ligands to explicit-H 3D SDF and prepare AutoDock/Vina-family PDBQT output."
    )
    parser.add_argument("--input", required=True, help="Input ligand file (.sdf, .mol, .mol2, .pdb)")
    parser.add_argument("--output", required=True, help="Output PDBQT path")
    parser.add_argument(
        "--backend-profile",
        default="engine_aware_full",
        help=(
            "Ligand preparation profile: "
            "openbabel_only, meeko_only, autodocktools_only, openbabel_meeko, "
            "openbabel_meeko_autodock, openbabel_autodocktools, engine_aware_full"
        ),
    )
    parser.add_argument(
        "--ph",
        type=float,
        default=DEFAULT_PREPARATION_PH,
        help=f"Protonation pH for Open Babel normalization stages ({MIN_PREPARATION_PH:.1f}-{MAX_PREPARATION_PH:.1f})",
    )
    parser.add_argument("--engines", default="", help="Comma-separated selected engines for engine-aware profile decisions")
    parser.add_argument("--summary-file", help="Optional JSON summary output path")
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)
    ph_ok, normalized_ph, ph_error = validate_preparation_ph(args.ph)
    if not ph_ok:
        print(f"❌ Invalid pH value: {ph_error}", file=sys.stderr)
        return 2
    os.environ["PDBWIZARD_LIGAND_PREP_PROFILE"] = _resolve_profile_alias(args.backend_profile)
    os.environ["PDBWIZARD_LIGAND_PREP_PH"] = str(normalized_ph)
    if args.engines:
        os.environ["PDBWIZARD_SELECTED_ENGINES"] = args.engines
    summary = prepare_ligand_for_vina_family(Path(args.input), Path(args.output))
    if args.summary_file:
        write_ligand_preparation_step_report(summary, report_file=Path(args.summary_file))
    else:
        print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
