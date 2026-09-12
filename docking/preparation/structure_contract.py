"""Lossless PDB selection and format-level contracts, without chemistry inference."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def coordinates(line: str) -> tuple[float, float, float]:
    try:
        xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
    except ValueError as exc:
        raise ValueError("Invalid fixed-column atom coordinates") from exc
    if not all(math.isfinite(value) for value in xyz):
        raise ValueError("Atom coordinates must be finite")
    return xyz


def residue_key(line: str) -> tuple[str, str, str, str]:
    return line[21:22].strip(), line[22:26].strip(), line[26:27].strip(), line[17:20].strip()


def select_pdb_lines(lines, *, residue=None, remove_instances=(), preferred_altloc="A"):
    """Select a single model and coherent residue conformers; preserve atom metadata.

    Ambiguous multi-model input must be explicitly split before this operation.
    Residue selectors are (chain, residue number, insertion code, residue name).
    """
    if sum(line.startswith("MODEL ") for line in lines) > 1:
        raise ValueError("Multiple models: select one explicit model before preparation")
    groups = {}
    removal = {tuple(str(v).strip() for v in key) for key in remove_instances}
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = residue_key(line)
        if residue is not None and key != tuple(str(v).strip() for v in residue):
            continue
        if key in removal:
            continue
        coordinates(line)
        groups.setdefault(key, []).append(line)
    selected = []
    selections = {}
    for key, atoms in groups.items():
        locations = {line[16:17].strip() for line in atoms} - {""}
        chosen = preferred_altloc if preferred_altloc in locations else None
        if locations and chosen is None:
            if len(locations) != 1:
                raise ValueError(f"Ambiguous alternate conformer for residue {key}: {sorted(locations)}")
            chosen = next(iter(locations))
        picked = [line for line in atoms if not line[16:17].strip() or line[16:17].strip() == chosen]
        names = [line[12:16].strip() for line in picked]
        all_names = {line[12:16].strip() for line in atoms}
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate atom names in selected residue {key}")
        if set(names) != all_names:
            raise ValueError(f"Incomplete alternate conformer {chosen} for residue {key}")
        selected.extend(line[:16] + " " + line[17:] for line in picked)
        selections[":".join(key)] = chosen or "shared"
    if not selected:
        raise ValueError("Selection contains no atoms")
    serials = {int(line[6:11]) for line in selected}
    if len(serials) != len(selected):
        raise ValueError("Duplicate atom serial identifiers")
    conect = []
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        values = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
        if values and values[0] in serials:
            # Repetitions can encode bond order and are intentionally preserved.
            neighbors = [value for value in values[1:] if value in serials]
            if neighbors:
                conect.append("CONECT" + "".join(f"{value:5d}" for value in [values[0], *neighbors]))
    return selected + conect + ["END"], selections


def write_selection(source: Path, destination: Path, *, residue=None, remove_instances=(), preferred_altloc="A"):
    source = Path(source)
    raw = source.read_bytes()
    original = raw.decode("utf-8", errors="strict").splitlines()
    output, selections = select_pdb_lines(original, residue=residue, remove_instances=remove_instances, preferred_altloc=preferred_altloc)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(output) + "\n", encoding="utf-8")
    atoms = [line for line in output if line.startswith(("ATOM  ", "HETATM"))]
    removed = {tuple(str(value).strip() for value in key) for key in remove_instances}
    eligible_count = sum(line.startswith(("ATOM  ", "HETATM")) and (residue is None or residue_key(line) == tuple(map(str, residue))) and residue_key(line) not in removed for line in original)
    metadata = {
        "source_file": str(source.resolve()), "source_sha256": hashlib.sha256(raw).hexdigest(),
        "destination_file": str(destination.resolve()), "atom_count": len(atoms),
        "removed_altloc_atoms": eligible_count - len(atoms),
        "removed_atom_count": sum(line.startswith(("ATOM  ", "HETATM")) for line in original) - len(atoms),
        "alternate_locations": selections,
        "chemical_identity_status": "not_evaluated_pdb_is_not_authoritative_chemistry",
        "reference_frame_id": "coordinates:" + hashlib.sha256("\n".join(line[30:54] for line in original if line.startswith(("ATOM  ", "HETATM"))).encode()).hexdigest(),
        "atom_map": [{"source_serial": int(line[6:11]), "output_serial": int(line[6:11]), "atom_name": line[12:16].strip(), "residue": residue_key(line)} for line in atoms],
    }
    destination.with_suffix(destination.suffix + ".preparation.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata
