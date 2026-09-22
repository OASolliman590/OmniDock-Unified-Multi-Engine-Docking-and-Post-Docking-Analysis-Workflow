"""Strict receptor preparation: no silent backend substitution or atom deletion."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path

from ..models import validate_preparation_ph
from .structure_contract import coordinates


_PDBQT_TYPE_TO_ELEMENT = {
    "A": "C",
    "C": "C",
    "H": "H",
    "HD": "H",
    "HS": "H",
    "D": "H",
    "N": "N",
    "NA": "N",
    "NS": "N",
    "O": "O",
    "OA": "O",
    "OS": "O",
    "S": "S",
    "SA": "S",
    "P": "P",
    "F": "F",
    "CL": "CL",
    "BR": "BR",
    "I": "I",
}
_PDBQT_METAL_TYPES = {
    "Na": "NA",
    "Ca": "CA",
    "Mg": "MG",
    "Mn": "MN",
    "Zn": "ZN",
    "Fe": "FE",
    "Cu": "CU",
}


def _pdbqt_element(token: str) -> str:
    token = str(token or "").strip()
    if token in _PDBQT_METAL_TYPES:
        return _PDBQT_METAL_TYPES[token]
    return _PDBQT_TYPE_TO_ELEMENT.get(token.upper(), "".join(char for char in token if char.isalpha()).upper())


def _source_element(line: str) -> str:
    field = str(line[76:78]).strip() if len(line) >= 78 else ""
    if field:
        return "".join(char for char in field if char.isalpha()).upper()
    atom_name = str(line[12:16]).strip()
    letters = "".join(char for char in atom_name if char.isalpha()).upper()
    if len(letters) >= 2 and letters[:2] == "CA":
        # ``CA`` is the standard alpha-carbon atom name as well as a calcium
        # element symbol.  Without the fixed-column element field, inferring
        # calcium would silently change the receptor chemistry.
        return ""
    if len(letters) >= 2 and letters[:2] in {"CL", "BR", "ZN", "FE", "MG", "MN", "CU", "NA", "PD"}:
        return letters[:2]
    return letters[:1]


def _atom_records(path: Path, *, pdbqt: bool = False):
    records = []
    serials = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            serial = int(line[6:11])
        except ValueError as exc:
            raise ValueError("Receptor atom serial identifiers must be integers") from exc
        if serial in serials:
            raise ValueError(f"Receptor atom serial identifiers must be unique (duplicate {serial})")
        serials.add(serial)
        element = _pdbqt_element(line.split()[-1]) if pdbqt else _source_element(line)
        if not element:
            raise ValueError("Receptor atom is missing an unambiguous element symbol")
        if element in {"H", "D", "HD", "HS"}:
            continue
        xyz = coordinates(line)
        identity = (
            line[21:22].strip(),
            line[22:26].strip(),
            line[26:27].strip(),
            line[17:20].strip(),
            line[12:16].strip(),
        )
        records.append(
            {
                "serial": serial,
                "identity": identity,
                "element": element,
                "coordinates": xyz,
                "coordinate_key": tuple(round(value, 3) for value in xyz),
            }
        )
    return records


def _heavy_positions(path: Path, pdbqt=False):
    """Compatibility view of heavy atoms at the PDB/PDBQT 3-decimal precision."""
    return Counter((record["element"], *record["coordinate_key"]) for record in _atom_records(path, pdbqt=pdbqt))


def _format_identity(identity) -> str:
    chain, residue_number, insertion, residue_name, atom_name = identity
    return f"{chain or '_'}:{residue_number or '_'}{insertion or ''}:{residue_name or '_'}:{atom_name or '_'}"


def _conserve_heavy_atoms(source_records, prepared_records):
    source_by_identity = {}
    prepared_by_identity = {}
    for record in source_records:
        identity = record["identity"]
        if identity in source_by_identity:
            raise ValueError(f"Receptor heavy-atom/coordinate conservation failed: duplicate source identity {_format_identity(identity)}")
        source_by_identity[identity] = record
    for record in prepared_records:
        identity = record["identity"]
        if identity in prepared_by_identity:
            raise ValueError(f"Receptor heavy-atom/coordinate conservation failed: duplicate prepared identity {_format_identity(identity)}")
        prepared_by_identity[identity] = record

    missing = sorted(set(source_by_identity) - set(prepared_by_identity))
    added = sorted(set(prepared_by_identity) - set(source_by_identity))
    if missing or added:
        missing_text = ", ".join(_format_identity(identity) for identity in missing[:5]) or "none"
        added_text = ", ".join(_format_identity(identity) for identity in added[:5]) or "none"
        raise ValueError(
            "Receptor heavy-atom/coordinate conservation failed: identity mismatch; "
            f"missing={missing_text}; added={added_text}"
        )

    atom_map = []
    for identity in sorted(source_by_identity):
        source_record = source_by_identity[identity]
        prepared_record = prepared_by_identity[identity]
        if source_record["element"] != prepared_record["element"]:
            raise ValueError(
                "Receptor heavy-atom/coordinate conservation failed: element mismatch "
                f"for {_format_identity(identity)} ({source_record['element']} -> {prepared_record['element']})"
            )
        if source_record["coordinate_key"] != prepared_record["coordinate_key"]:
            raise ValueError(
                "Receptor heavy-atom/coordinate conservation failed: identity/coordinate mismatch "
                f"for {_format_identity(identity)} at represented 3-decimal precision"
            )
        atom_map.append(
            {
                "source_serial": source_record["serial"],
                "output_serial": prepared_record["serial"],
                "atom_name": identity[-1],
                "residue": list(identity[:4]),
                "element": source_record["element"],
                "source_coordinates": list(source_record["coordinates"]),
                "output_coordinates": list(prepared_record["coordinates"]),
            }
        )
    return atom_map


def prepare_receptor(source: Path, destination: Path, config: dict):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    preparation = config.get("preparation", {})
    requested_ph = preparation.get("ph", 7.4)
    ph_ok, normalized_ph, ph_error = validate_preparation_ph(requested_ph)
    if not ph_ok or not math.isfinite(normalized_ph):
        raise ValueError(ph_error or "pH must be finite")
    if preparation.get("allow_bad_res", False):
        raise ValueError("Permissive bad-residue deletion is unsupported: repair or explicitly remove individual residues first")
    lines = source.read_text(encoding="utf-8").splitlines()
    if sum(line.startswith("MODEL ") for line in lines) > 1:
        raise ValueError("Select one receptor model before preparation")
    if any(line.startswith(("ATOM  ", "HETATM")) and line[16:17].strip() for line in lines):
        raise ValueError("Select consistent alternate conformers before receptor preparation")
    source_records = _atom_records(source)
    if not source_records:
        raise ValueError("Receptor contains no finite heavy-atom coordinates")
    logs = []
    commands = []
    def run(command):
        commands.append(command)
        process = subprocess.run(command, capture_output=True, text=True)
        logs.append(process.stdout + "\n" + process.stderr)
        if process.returncode:
            raise RuntimeError(f"Required receptor preparation stage failed: {command[0]}\n{process.stderr}")
        return process
    try:
        with tempfile.TemporaryDirectory(prefix="receptor_prepare_", dir=destination.parent) as temporary:
            folder = Path(temporary)
            prepared = folder / "receptor.pdbqt"
            prep_input = source
            ph_status = "template_selected_not_ph_titrated"
            upstream_titration_status = "not_requested"
            if preparation.get("receptor_use_pdb2pqr", False):
                pqr, converted = folder / "protonated.pqr", folder / "protonated.pdb"
                run([
                    "pdb2pqr30",
                    "--ff",
                    str(preparation.get("force_field", "AMBER")),
                    "--titration-state-method",
                    "propka",
                    "--with-ph",
                    str(normalized_ph),
                    "--keep-chain",
                    str(source),
                    str(pqr),
                ])
                run(["obabel", str(pqr), "-O", str(converted)])
                prep_input = converted
                upstream_titration_status = "pdb2pqr_propka_applied"
                ph_status = "backend_state_not_ph_validated"
            profile = preparation.get("ligand_preparation_profile", preparation.get("ligand_preparation_backend", "engine_aware_full"))
            if profile == "autodocktools_only":
                script = preparation.get("autodocktools_prepare_receptor4")
                if not script:
                    raise ValueError("AutoDockTools receptor script is required for the requested profile")
                run([preparation.get("autodocktools_python") or "python3", str(script), "-r", str(prep_input), "-o", str(prepared), "-A", "checkhydrogens", "-U", "nphs_lps"])
                backend = "autodocktools"
                if upstream_titration_status == "not_requested":
                    ph_status = "input_state_not_ph_titrated"
            else:
                run(["mk_prepare_receptor.py", "--read_pdb", str(prep_input), "-p", str(prepared), "-j", str(folder / "receptor.json")])
                backend = "meeko"
            if not prepared.exists():
                raise ValueError("Receptor backend did not create its requested output")
            prepared_records = _atom_records(prepared, pdbqt=True)
            atom_map = _conserve_heavy_atoms(source_records, prepared_records)
            from .receptor_quality import validate_receptor_file
            issues = validate_receptor_file(prepared, min_atom_count=1, min_heavy_atom_count=1)
            if issues:
                raise ValueError("Receptor contract failed: " + "; ".join(issue.details for issue in issues))
            shutil.copy2(prepared, destination)
            prior = source.with_suffix(source.suffix + ".preparation.json")
            provenance = json.loads(prior.read_text(encoding="utf-8")) if prior.exists() else {}
            try:
                tool_version = importlib.metadata.version("meeko" if backend == "meeko" else "AutoDockTools")
            except importlib.metadata.PackageNotFoundError:
                tool_version = "not_available"
            payload = {"input_file": str(source), "output_file": str(destination), "backend": backend,
                       "tool_version": tool_version,
                       "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                       "output_sha256": hashlib.sha256(destination.read_bytes()).hexdigest(),
                       "preparation_config": preparation, "commands": commands,
                       "requested_ph": normalized_ph, "protonation_status": ph_status,
                       "upstream_titration_status": upstream_titration_status,
                       "final_backend_state_status": ph_status,
                       "final_pH_validated": False,
                       "receptor_frame_id": provenance.get("receptor_frame_id", provenance.get("reference_frame_id", "")),
                       "heavy_atom_conservation": "passed", "atom_identity_conservation": "passed",
                       "coordinate_precision_decimals": 3, "atom_map": atom_map, "is_valid": True}
            destination.with_suffix(destination.suffix + ".preparation.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
            return payload
    finally:
        destination.with_suffix(destination.suffix + ".preparation.log").write_text("\n".join(logs), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    prepare_receptor(args.source, args.destination, json.loads(args.config.read_text(encoding="utf-8")))


if __name__ == "__main__":
    main()
