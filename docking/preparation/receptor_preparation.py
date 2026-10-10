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

Heavy-atom policy (Spec 033 decision, option 1, Scientific Lead 2026-10-08): with PDB2PQR
the run uses --noopt --nodebump, so no side chain is flipped. Removed or moved heavy atoms
always fail. An added heavy atom is accepted only as a C-terminal carboxylate O (name O or
OXT) on the last residue of a chain, at most one per terminus; it is recorded in provenance.
Any other addition fails. The explicit receptor_use_pdb2pqr=false path accepts no additions.
PDB2PQR removes HETATM records (ligands, ions, waters), so those fail conservation.

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
from typing import Any, Dict, List, NamedTuple, Optional, Sequence, Tuple

from .receptor_quality import validate_receptor_file
from .structure_contract import coordinates

PDB2PQR_EXECUTABLE = "pdb2pqr30"
OPEN_BABEL_EXECUTABLE = "obabel"
MEEKO_EXECUTABLE = "mk_prepare_receptor.py"
PDB2PQR_FORCE_FIELDS = ("AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON")
# Suggestions shown to the user only (interactive prompt). Preparation never falls back to them (Spec 034 R3a).
SUGGESTED_PH = 7.4
SUGGESTED_FORCE_FIELD = "AMBER"
# Where a pH or force field came from (Spec 034 R3b): an explicit CLI/prompt entry, or an explicit config value.
PARAMETER_SOURCES = ("user_entered", "config_file")
HD_NITROGEN_MAX_DISTANCE = 1.1  # Angstrom; N-H bond length ceiling for HD-N association
# Spec 033 decision (option 1): no side-chain flip or debump optimisation. --keep-chain only labels
# chains so that terminal additions can be attributed to a chain; it does not move atoms.
PDB2PQR_OPTIONS = ("--noopt", "--nodebump", "--keep-chain")
# Carboxylate oxygens that PDB2PQR may add at a C-terminus (names as in the PDB2PQR and Meeko outputs).
TERMINAL_CARBOXYLATE_ATOM_NAMES = ("O", "OXT")

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


def _required_parameter(preparation: Dict[str, Any], key: str, reason: str, guidance: str) -> Any:
    value = preparation.get(key)
    if value is None or (isinstance(value, str) and not value.strip()):
        raise ReceptorPreparationError(reason, f"preparation.{key} is required for PDB2PQR protonation. {guidance}")
    return value


def _parameter_source(preparation: Dict[str, Any], key: str) -> str:
    source = preparation.get(f"{key}_source")
    if source not in PARAMETER_SOURCES:
        raise ReceptorPreparationError(
            f"{key}_source_missing",
            f"preparation.{key}_source must be 'user_entered' (CLI or prompt) or 'config_file' "
            f"(an explicit configuration value), got {source!r}",
        )
    return str(source)


def _atom_element(line: str, pdbqt: bool) -> str:
    element = line.split()[-1] if pdbqt else line[76:78].strip()
    if not element:
        element = line[12:16].strip().lstrip("0123456789")[:1]
    return element.upper()


class _HeavyAtom(NamedTuple):
    signature: Tuple[str, float, float, float]  # (element, x, y, z) rounded to 0.01 A
    chain: str
    resname: str
    resseq: int
    atom: str

    @property
    def label(self) -> str:
        return f"{self.resname} {self.resseq}:{self.atom}"


def _heavy_atoms(path: Path, pdbqt: bool = False) -> List[_HeavyAtom]:
    atoms: List[_HeavyAtom] = []
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
            atoms.append(_HeavyAtom(
                signature=(element, *xyz),
                chain=line[21:22].strip(),
                resname=line[17:20].strip(),
                resseq=int(line[22:26]),
                atom=line[12:16].strip(),
            ))
    except ValueError as exc:
        raise ReceptorPreparationError("invalid_coordinates", f"{path.name}: {exc}") from exc
    return atoms


def _unmatched(source_atoms: List[_HeavyAtom], target_atoms: List[_HeavyAtom]) -> List[_HeavyAtom]:
    """Source heavy atoms that have no identical (element, position) counterpart in target."""
    remaining = Counter(atom.signature for atom in target_atoms)
    unmatched = []
    for atom in source_atoms:
        if remaining[atom.signature] > 0:
            remaining[atom.signature] -= 1
        else:
            unmatched.append(atom)
    return unmatched


def _chain_termini(source_path: Path) -> Dict[str, int]:
    """Highest residue number per chain in the input ATOM records (C-terminal residue).

    Internal chain breaks are not termini: only the highest residue number counts.
    """
    termini: Dict[str, int] = {}
    try:
        for line in Path(source_path).read_text(encoding="utf-8").splitlines():
            if line.startswith("ATOM  "):
                chain = line[21:22].strip()
                termini[chain] = max(termini.get(chain, int(line[22:26])), int(line[22:26]))
    except ValueError as exc:
        raise ReceptorPreparationError("invalid_coordinates", f"{Path(source_path).name}: {exc}") from exc
    return termini


def check_heavy_atom_conservation(
    source_atoms: List[_HeavyAtom],
    prepared_atoms: List[_HeavyAtom],
    termini: Dict[str, int],
    allow_terminal_additions: bool,
) -> List[Dict[str, Any]]:
    """Strict heavy-atom conservation (Spec 033 decision, option 1).

    Removed or moved source heavy atoms always fail. An added heavy atom is accepted only
    when allow_terminal_additions is set (PDB2PQR path), it is an oxygen named O or OXT on the
    last residue of its chain, and there is at most one such addition per chain terminus.
    Returns the accepted terminal additions. Raises ReceptorPreparationError otherwise.
    """
    removed = _unmatched(source_atoms, prepared_atoms)
    added = _unmatched(prepared_atoms, source_atoms)
    if removed:
        raise ReceptorPreparationError(
            "heavy_atom_conservation_failed",
            f"{len(removed)} source heavy atom(s) missing or moved. Missing/moved: {[a.label for a in removed[:20]]}. "
            f"Added/moved: {[a.label for a in added[:20]]}",
        )
    accepted: List[_HeavyAtom] = []
    rejected: List[_HeavyAtom] = []
    for atom in added:
        is_terminal_oxygen = (
            allow_terminal_additions
            and atom.signature[0] == "O"
            and atom.atom in TERMINAL_CARBOXYLATE_ATOM_NAMES
            and termini.get(atom.chain) == atom.resseq
        )
        (accepted if is_terminal_oxygen else rejected).append(atom)
    if rejected:
        raise ReceptorPreparationError(
            "heavy_atom_conservation_failed",
            f"{len(added)} output heavy atom(s) added; only a C-terminal carboxylate O on the last residue of a chain "
            f"is accepted. Rejected additions: {[a.label for a in rejected[:20]]}",
        )
    per_terminus = Counter((atom.chain, atom.resseq) for atom in accepted)
    overfilled = sorted(key for key, count in per_terminus.items() if count > 1)
    if overfilled:
        raise ReceptorPreparationError(
            "heavy_atom_conservation_failed",
            f"more than one added carboxylate O at terminus (chain, residue) {overfilled}",
        )
    return [
        {"chain": atom.chain, "residue": atom.resname, "resseq": atom.resseq, "atom": atom.atom, "element": "O"}
        for atom in accepted
    ]


def check_receptor_conservation(
    source_path: Path, prepared_path: Path, allow_terminal_additions: bool
) -> List[Dict[str, Any]]:
    """Read a source PDB and a prepared PDBQT, then apply check_heavy_atom_conservation."""
    return check_heavy_atom_conservation(
        _heavy_atoms(source_path),
        _heavy_atoms(prepared_path, pdbqt=True),
        _chain_termini(source_path),
        allow_terminal_additions,
    )

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
    # Spec 034 R3a: no silent pH or force field. PDB2PQR needs both as explicit values with a recorded source.
    ph: Optional[float] = None
    ph_source: Optional[str] = None
    force_field = ""
    force_field_source: Optional[str] = None
    if use_pdb2pqr:
        ph_guidance = "Enter it with --ph or in the configuration (preparation.ph and preparation.ph_source)."
        ph = _ph_option(_required_parameter(preparation, "ph", "ph_required", ph_guidance))
        ph_source = _parameter_source(preparation, "ph")
        force_field_guidance = "Choose one of: " + ", ".join(PDB2PQR_FORCE_FIELDS) + "."
        force_field = str(
            _required_parameter(preparation, "force_field", "force_field_required", force_field_guidance)
        ).strip().upper()
        force_field_source = _parameter_source(preparation, "force_field")
        if force_field not in PDB2PQR_FORCE_FIELDS:
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
            run([PDB2PQR_EXECUTABLE, "--ff", force_field, "--with-ph", str(ph), *PDB2PQR_OPTIONS, str(source), str(pqr)], "pdb2pqr_failed")
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
        accepted_additions = check_heavy_atom_conservation(
            original_atoms,
            prepared_atoms,
            _chain_termini(source),
            allow_terminal_additions=use_pdb2pqr,
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
        "ph_source": ph_source,
        "requested_force_field": force_field or None,
        "force_field_source": force_field_source,
        "protonation_status": protonation_status,
        "receptor_frame_id": prior.get("receptor_frame_id", prior.get("reference_frame_id", "")),
        "pdb2pqr_options": list(PDB2PQR_OPTIONS) if use_pdb2pqr else [],
        "accepted_terminal_additions": accepted_additions,
        "heavy_atom_conservation": "passed_with_terminal_additions" if accepted_additions else "passed",
        "heavy_atom_count": len(original_atoms),
        "output_heavy_atom_count": len(prepared_atoms),
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
