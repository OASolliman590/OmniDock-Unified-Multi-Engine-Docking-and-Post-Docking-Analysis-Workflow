"""Strict receptor preparation for Vina-family PDBQT receptors (Spec 033 R2).

Ported from main@7fc130a and extended for Spec 033:

- No silent backend substitution. A missing PDB2PQR, Open Babel or Meeko
  executable is a failure with reason ``backend_unavailable``.
- PDB2PQR protonation is on by default at an explicit pH and force field.
  Both are recorded in provenance. Disabling it requires an explicit JSON
  ``false`` and is recorded as ``template_selected_not_ph_titrated``.
- Heavy-atom coordinate conservation is checked before the output is accepted.
- Polar-hydrogen gate: the PDBQT must contain at least one ``HD`` atom, and at
  least one ``HD`` must lie within 1.1 A of a nitrogen donor (``N``/``NA``).
- Provenance JSON (backend, tool versions, commands, hashes, pH, protonation
  status) and a preparation log are written next to the output.

Scientific limitations (disclosed, not solved here): PDB2PQR removes HETATM
records, and its hydrogen-bond optimisation can flip His/Asn/Gln side chains.
It also adds terminal atoms such as C-terminal OXT. Any of these changes a
heavy atom and fails the conservation gate, by design.

The CLI prints one JSON object on stdout. Success is ``{"status": "completed"}``
and exit code 0. Failure is ``{"status": "failed", "reason": ..., "details": ...}``
and exit code 1, so the shell wrapper can record the reason.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .receptor_quality import validate_receptor_file
from .structure_contract import coordinates

PDB2PQR_EXECUTABLE = "pdb2pqr30"
OPEN_BABEL_EXECUTABLE = "obabel"
MEEKO_EXECUTABLE = "mk_prepare_receptor.py"
PDB2PQR_FORCE_FIELDS = ("AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON")
DEFAULT_PH = 7.4
DEFAULT_FORCE_FIELD = "AMBER"
HD_NITROGEN_MAX_DISTANCE = 1.1  # Angstrom; N-H bond length ceiling for HD-N association

PROTONATION_PDB2PQR = "pdb2pqr_applied"
PROTONATION_TEMPLATE = "template_selected_not_ph_titrated"
PROTONATION_AUTODOCKTOOLS = "input_state_not_ph_titrated"

_PDBQT_TYPE_TO_ELEMENT = {"A": "C", "NA": "N", "NS": "N", "OA": "O", "OS": "O", "SA": "S"}
_HYDROGEN_ELEMENTS = {"H", "HD", "HS", "D"}


class ReceptorPreparationError(RuntimeError):
    """Preparation refused. ``reason`` is the machine-readable failure reason."""

    def __init__(self, reason: str, details: str) -> None:
        super().__init__(f"{reason}: {details}")
        self.reason = reason
        self.details = details


def _distribution_version(executable: str, distribution: str) -> Optional[str]:
    if shutil.which(executable) is None:
        return None
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def receptor_backend_capabilities() -> Dict[str, Optional[str]]:
    """Versions of the receptor backends found on PATH; None when not runnable."""
    return {
        "pdb2pqr": _distribution_version(PDB2PQR_EXECUTABLE, "pdb2pqr"),
        "meeko": _distribution_version(MEEKO_EXECUTABLE, "meeko"),
    }


def _bool_option(preparation: Dict[str, Any], key: str, default: bool) -> bool:
    value = preparation.get(key, default)
    if not isinstance(value, bool):
        raise ReceptorPreparationError("invalid_configuration", f"preparation.{key} must be a JSON boolean, got {value!r}")
    return value


def _ph_option(value: Any) -> float:
    if isinstance(value, bool):
        raise ReceptorPreparationError("invalid_configuration", "preparation.ph must be a number")
    try:
        ph = float(value)
    except (TypeError, ValueError) as exc:
        raise ReceptorPreparationError("invalid_configuration", f"preparation.ph must be a number, got {value!r}") from exc
    if not math.isfinite(ph) or not 0.0 <= ph <= 14.0:
        raise ReceptorPreparationError("invalid_configuration", f"preparation.ph must be within [0, 14], got {ph}")
    return ph


def _atom_element(line: str, pdbqt: bool) -> str:
    element = line.split()[-1] if pdbqt else line[76:78].strip()
    if not element:
        element = line[12:16].strip().lstrip("0123456789")[:1]
    return element.upper()


def _heavy_atoms(path: Path, pdbqt: bool = False) -> List[Tuple[Tuple[str, float, float, float], str]]:
    """Heavy atoms as (signature, label); signature = (element, x, y, z) at 0.01 A."""
    atoms = []
    try:
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            element = _atom_element(line, pdbqt)
            if element in _HYDROGEN_ELEMENTS:
                continue
            if pdbqt:
                element = _PDBQT_TYPE_TO_ELEMENT.get(element, element)
            xyz = [round(value, 2) for value in coordinates(line)]
            label = f"{line[17:20].strip()} {line[22:26].strip()}:{line[12:16].strip()}"
            atoms.append(((element, *xyz), label))
    except ValueError as exc:
        raise ReceptorPreparationError("invalid_coordinates", f"{path.name}: {exc}") from exc
    return atoms


def _unmatched(source_atoms, target_atoms) -> List[str]:
    """Labels of source heavy atoms that have no identical counterpart in target."""
    remaining = Counter(signature for signature, _ in target_atoms)
    unmatched = []
    for signature, label in source_atoms:
        if remaining[signature] > 0:
            remaining[signature] -= 1
        else:
            unmatched.append(label)
    return unmatched


def polar_hydrogen_gate(pdbqt_path: Path) -> Dict[str, Any]:
    """Count HD atoms and HD atoms within 1.1 A of a nitrogen donor (N or NA)."""
    hydrogens: List[Tuple[float, float, float]] = []
    nitrogens: List[Tuple[float, float, float]] = []
    try:
        for line in Path(pdbqt_path).read_text(encoding="utf-8").splitlines():
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom_type = line.split()[-1].upper()
            if atom_type == "HD":
                hydrogens.append(coordinates(line))
            elif atom_type in {"N", "NA"}:
                nitrogens.append(coordinates(line))
    except ValueError as exc:
        raise ReceptorPreparationError("invalid_coordinates", f"{Path(pdbqt_path).name}: {exc}") from exc
    near_nitrogen = sum(
        1 for hydrogen in hydrogens
        if any(math.dist(hydrogen, nitrogen) <= HD_NITROGEN_MAX_DISTANCE for nitrogen in nitrogens)
    )
    return {
        "hd_count": len(hydrogens),
        "hd_near_nitrogen_count": near_nitrogen,
        "passed": len(hydrogens) >= 1 and near_nitrogen >= 1,
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _autodocktools_script(preparation: Dict[str, Any]) -> str:
    script = str(preparation.get("autodocktools_prepare_receptor4") or os.environ.get("AUTODOCKTOOLS_PREPARE_RECEPTOR4", "")).strip()
    if not script or not Path(script).is_file():
        raise ReceptorPreparationError(
            "backend_unavailable",
            "AutoDockTools prepare_receptor4.py is required by the autodocktools_only profile but was not found",
        )
    return script


def _autodocktools_python(preparation: Dict[str, Any]) -> str:
    python = str(preparation.get("autodocktools_python") or os.environ.get("AUTODOCKTOOLS_PYTHON", "") or "python3").strip()
    if shutil.which(python) is None and not Path(python).is_file():
        raise ReceptorPreparationError("backend_unavailable", f"AutoDockTools python interpreter not found: {python}")
    return python


def prepare_receptor(source: Path, destination: Path, config: Dict[str, Any]) -> Dict[str, Any]:
    """Prepare one receptor. Raises ReceptorPreparationError; never writes a partial output."""
    source, destination = Path(source).resolve(), Path(destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    log_path = destination.parent / (destination.name + ".preparation.log")
    # Always regenerate: a stale receptor or provenance must never survive a failed run.
    for stale in (destination, destination.parent / (destination.name + ".preparation.json")):
        if stale.exists():
            stale.unlink()
    log_lines: List[str] = []
    try:
        return _prepare(source, destination, dict(config.get("preparation") or {}), log_lines)
    finally:
        log_path.write_text("\n".join(log_lines), encoding="utf-8")


def _prepare(source: Path, destination: Path, preparation: Dict[str, Any], log_lines: List[str]) -> Dict[str, Any]:
    if _bool_option(preparation, "allow_bad_res", False):
        raise ReceptorPreparationError(
            "allow_bad_res_unsupported",
            "Permissive bad-residue deletion is unsupported: repair or explicitly remove individual residues first",
        )
    use_pdb2pqr = _bool_option(preparation, "receptor_use_pdb2pqr", True)
    ph = _ph_option(preparation.get("ph", DEFAULT_PH))
    force_field = str(preparation.get("force_field", DEFAULT_FORCE_FIELD)).strip().upper() or DEFAULT_FORCE_FIELD
    if use_pdb2pqr and force_field not in PDB2PQR_FORCE_FIELDS:
        raise ReceptorPreparationError("invalid_configuration", f"preparation.force_field {force_field!r} is not a PDB2PQR force field")
    profile = str(
        preparation.get("ligand_preparation_profile", preparation.get("ligand_preparation_backend", "engine_aware_full")) or "engine_aware_full"
    ).strip().lower()
    backend = "autodocktools" if profile == "autodocktools_only" else "meeko"

    lines = source.read_text(encoding="utf-8").splitlines()
    if sum(line.startswith("MODEL ") for line in lines) > 1:
        raise ReceptorPreparationError("multiple_models", "Select one receptor model before preparation")
    if any(line.startswith(("ATOM  ", "HETATM")) and line[16:17].strip() for line in lines):
        raise ReceptorPreparationError("altloc_ambiguous", "Select consistent alternate conformers before receptor preparation")
    original_atoms = _heavy_atoms(source)
    if not original_atoms:
        raise ReceptorPreparationError("no_heavy_atoms", "Receptor contains no finite heavy-atom coordinates")

    # Fail closed before running anything: every required executable must exist.
    required: List[str] = []
    if use_pdb2pqr:
        required += [PDB2PQR_EXECUTABLE, OPEN_BABEL_EXECUTABLE]
    if backend == "autodocktools":
        adt_script = _autodocktools_script(preparation)
        adt_python = _autodocktools_python(preparation)
    else:
        required.append(MEEKO_EXECUTABLE)
    missing = [name for name in required if shutil.which(name) is None]
    if missing:
        raise ReceptorPreparationError(
            "backend_unavailable",
            "Required receptor backend executable(s) not found on PATH: " + ", ".join(missing)
            + ". Install conda-forge pdb2pqr>=3.4.1 and meeko>=0.7.1 (see environment.yml).",
        )
    executables = {name: shutil.which(name) for name in required}

    commands: List[List[str]] = []

    def run(command: List[str], stage_reason: str, env: Optional[Dict[str, str]] = None) -> None:
        commands.append(list(command))
        log_lines.append("$ " + " ".join(command))
        try:
            process = subprocess.run(command, capture_output=True, text=True, env=env)
        except OSError as exc:
            raise ReceptorPreparationError("backend_unavailable", f"{command[0]} could not be executed: {exc}") from exc
        log_lines.append(process.stdout)
        log_lines.append(process.stderr)
        if process.returncode:
            tail = (process.stderr or process.stdout or "").strip().splitlines()[-1:] or [""]
            raise ReceptorPreparationError(stage_reason, f"{command[0]} exited with {process.returncode}: {tail[0]}")

    with tempfile.TemporaryDirectory(prefix="receptor_prepare_", dir=destination.parent) as temporary:
        folder = Path(temporary)
        prepared = folder / "receptor.pdbqt"
        prep_input = source
        protonation_status = PROTONATION_TEMPLATE
        if use_pdb2pqr:
            pqr, converted = folder / "protonated.pqr", folder / "protonated.pdb"
            run([PDB2PQR_EXECUTABLE, "--ff", force_field, "--with-ph", str(ph), str(source), str(pqr)], "pdb2pqr_failed")
            run([OPEN_BABEL_EXECUTABLE, str(pqr), "-O", str(converted)], "pqr_conversion_failed")
            prep_input = converted
            protonation_status = PROTONATION_PDB2PQR

        if backend == "autodocktools":
            env = dict(os.environ)
            try:
                adt_root = Path(adt_script).resolve().parents[2]
                env["PYTHONPATH"] = f"{adt_root}{os.pathsep}{env.get('PYTHONPATH', '')}".rstrip(os.pathsep)
            except IndexError:
                pass
            run([adt_python, adt_script, "-r", str(prep_input), "-o", str(prepared), "-A", "checkhydrogens", "-U", "nphs_lps"],
                "autodocktools_failed", env=env)
            if not use_pdb2pqr:
                protonation_status = PROTONATION_AUTODOCKTOOLS
        else:
            run([MEEKO_EXECUTABLE, "--read_pdb", str(prep_input), "-p", str(prepared), "-j", str(folder / "receptor.json")], "meeko_failed")

        if not prepared.exists():
            raise ReceptorPreparationError("backend_output_missing", f"{backend} backend did not create its PDBQT output")

        prepared_atoms = _heavy_atoms(prepared, pdbqt=True)
        removed = _unmatched(original_atoms, prepared_atoms)
        added = _unmatched(prepared_atoms, original_atoms)
        if removed or added:
            raise ReceptorPreparationError(
                "heavy_atom_conservation_failed",
                f"{len(removed)} source heavy atom(s) missing or moved and {len(added)} output heavy atom(s) added or moved. "
                f"Missing/moved: {removed[:20]}. Added/moved: {added[:20]}",
            )

        issues = validate_receptor_file(prepared, min_atom_count=1, min_heavy_atom_count=1)
        if issues:
            raise ReceptorPreparationError("receptor_contract_failed", "; ".join(issue.details for issue in issues))

        gate = polar_hydrogen_gate(prepared)
        if not gate["passed"]:
            raise ReceptorPreparationError(
                "polar_hydrogen_gate_failed",
                f"Vina-family receptor needs >=1 HD atom and >=1 HD within {HD_NITROGEN_MAX_DISTANCE} A of N/NA; "
                f"found hd_count={gate['hd_count']}, hd_near_nitrogen_count={gate['hd_near_nitrogen_count']}",
            )

        shutil.copy2(prepared, destination)

    prior_path = source.parent / (source.name + ".preparation.json")
    try:
        prior = json.loads(prior_path.read_text(encoding="utf-8")) if prior_path.exists() else {}
    except (OSError, ValueError):
        prior = {}
    capabilities = receptor_backend_capabilities()
    if backend == "meeko":
        tool_version = capabilities["meeko"] or "not_available"
    else:
        try:
            tool_version = importlib.metadata.version("AutoDockTools")
        except importlib.metadata.PackageNotFoundError:
            tool_version = "not_available"
    payload = {
        "input_file": str(source),
        "output_file": str(destination),
        "backend": backend,
        "tool_version": tool_version,
        "tool_versions": capabilities,
        "executables": executables,
        "input_sha256": _sha256(source),
        "output_sha256": _sha256(destination),
        "preparation_config": preparation,
        "commands": commands,
        "requested_ph": ph,
        "requested_force_field": force_field,
        "protonation_status": protonation_status,
        "receptor_frame_id": prior.get("receptor_frame_id", prior.get("reference_frame_id", "")),
        "heavy_atom_conservation": "passed",
        "heavy_atom_count": len(original_atoms),
        "hd_count": gate["hd_count"],
        "hd_near_nitrogen_count": gate["hd_near_nitrogen_count"],
        "polar_hydrogen_gate": "passed",
        "is_valid": True,
    }
    provenance_path = destination.parent / (destination.name + ".preparation.json")
    provenance_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Strict receptor preparation (Spec 033 R2).")
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("config", type=Path)
    args = parser.parse_args(argv)
    try:
        config = json.loads(args.config.read_text(encoding="utf-8"))
        payload = prepare_receptor(args.source, args.destination, config)
    except ReceptorPreparationError as exc:
        print(json.dumps({"status": "failed", "reason": exc.reason, "details": exc.details}))
        return 1
    except Exception as exc:  # fail closed with a machine-readable reason
        print(json.dumps({"status": "failed", "reason": "preparation_error", "details": f"{type(exc).__name__}: {exc}"}))
        return 1
    print(json.dumps({
        "status": "completed",
        "output_file": payload["output_file"],
        "backend": payload["backend"],
        "protonation_status": payload["protonation_status"],
        "hd_count": payload["hd_count"],
        "hd_near_nitrogen_count": payload["hd_near_nitrogen_count"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
