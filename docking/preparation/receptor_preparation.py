"""Strict receptor preparation: no silent backend substitution or atom deletion."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

from .structure_contract import coordinates


def _heavy_positions(path: Path, pdbqt=False):
    positions = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line.split()[-1] if pdbqt else line[76:78].strip()
        if not element:
            element = line[12:16].strip().lstrip("0123456789")[:1]
        if element.upper() in {"H", "HD", "HS", "D"}:
            continue
        atom_elements = {"A": "C", "NA": "N", "NS": "N", "OA": "O", "OS": "O", "SA": "S"}
        if pdbqt:
            element = atom_elements.get(element, element)
        positions.append((element.upper(), *[round(value, 2) for value in coordinates(line)]))
    return Counter(positions)


def prepare_receptor(source: Path, destination: Path, config: dict):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    preparation = config.get("preparation", {})
    if preparation.get("allow_bad_res", False):
        raise ValueError("Permissive bad-residue deletion is unsupported: repair or explicitly remove individual residues first")
    lines = source.read_text(encoding="utf-8").splitlines()
    if sum(line.startswith("MODEL ") for line in lines) > 1:
        raise ValueError("Select one receptor model before preparation")
    if any(line.startswith(("ATOM  ", "HETATM")) and line[16:17].strip() for line in lines):
        raise ValueError("Select consistent alternate conformers before receptor preparation")
    original_atoms = _heavy_positions(source)
    if not original_atoms:
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
            if preparation.get("receptor_use_pdb2pqr", False):
                pqr, converted = folder / "protonated.pqr", folder / "protonated.pdb"
                run(["pdb2pqr30", "--ff", str(preparation.get("force_field", "AMBER")), "--with-ph", str(preparation.get("ph", 7.4)), str(source), str(pqr)])
                run(["obabel", str(pqr), "-O", str(converted)])
                prep_input = converted
                ph_status = "pdb2pqr_applied"
            profile = preparation.get("ligand_preparation_profile", preparation.get("ligand_preparation_backend", "engine_aware_full"))
            if profile == "autodocktools_only":
                script = preparation.get("autodocktools_prepare_receptor4")
                if not script:
                    raise ValueError("AutoDockTools receptor script is required for the requested profile")
                run([preparation.get("autodocktools_python") or "python3", str(script), "-r", str(prep_input), "-o", str(prepared), "-A", "checkhydrogens", "-U", "nphs_lps"])
                backend = "autodocktools"
                if ph_status != "pdb2pqr_applied":
                    ph_status = "input_state_not_ph_titrated"
            else:
                run(["mk_prepare_receptor.py", "--read_pdb", str(prep_input), "-p", str(prepared), "-j", str(folder / "receptor.json")])
                backend = "meeko"
            if not prepared.exists():
                raise ValueError("Receptor backend did not create its requested output")
            prepared_atoms = _heavy_positions(prepared, pdbqt=True)
            removed, added = original_atoms - prepared_atoms, prepared_atoms - original_atoms
            if removed or added:
                raise ValueError(f"Receptor heavy-atom/coordinate conservation failed: {sum(removed.values())} missing, {sum(added.values())} added or moved; inspect the preparation log")
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
                       "requested_ph": preparation.get("ph", 7.4), "protonation_status": ph_status,
                       "receptor_frame_id": provenance.get("receptor_frame_id", provenance.get("reference_frame_id", "")),
                       "heavy_atom_conservation": "passed", "is_valid": True}
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
