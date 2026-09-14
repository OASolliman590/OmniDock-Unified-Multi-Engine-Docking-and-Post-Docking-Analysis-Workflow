"""Lossless PDB selection and format-level contracts, without chemistry inference."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def _atom_identity(line: str) -> tuple[tuple[str, str, str, str], str]:
    """Return the PDB residue/atom identity without using serial numbers."""
    return residue_key(line), line[12:16].strip()


def _link_endpoint(line: str, start: int) -> tuple[tuple[str, str, str, str], str, str]:
    """Read one fixed-column LINK endpoint.

    ``start`` is the first column of the four-character atom name (12 or 42
    in zero-based indexing for the PDB v3.3 LINK record).
    """
    atom_name = line[start:start + 4].strip()
    altloc = line[start + 4:start + 5].strip()
    resname = line[start + 5:start + 8].strip()
    chain = line[start + 9:start + 10].strip()
    resnum = line[start + 10:start + 14].strip()
    insertion = line[start + 14:start + 15].strip()
    if not atom_name or not resname or not resnum:
        raise ValueError("Malformed LINK record cannot be safely selected")
    return (chain, resnum, insertion, resname), atom_name, altloc


def _ssbond_endpoint(line: str, start: int) -> tuple[str, str, str, str]:
    """Read one fixed-column SSBOND residue endpoint."""
    resname = line[start:start + 3].strip()
    chain = line[start + 4:start + 5].strip()
    resnum = line[start + 6:start + 10].strip()
    insertion = line[start + 10:start + 11].strip()
    if not resname or not resnum:
        raise ValueError("Malformed SSBOND record cannot be safely selected")
    return chain, resnum, insertion, resname


def _has_non_identity_symmetry(line: str) -> bool:
    """Return whether a LINK/SSBOND endpoint uses a non-local symmetry mate."""
    symmetry_1 = line[59:65].strip()
    symmetry_2 = line[66:72].strip()
    return any(value not in {"", "1555"} for value in (symmetry_1, symmetry_2))


def _normalize_link_altlocs(line: str) -> str:
    """Match selected atom records whose alternate-location fields are blanked."""
    chars = list(line.ljust(80))
    chars[16] = " "
    chars[46] = " "
    return "".join(chars)


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
    all_serials = set()
    for index, line in enumerate(lines):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            serial = int(line[6:11])
        except ValueError as exc:
            raise ValueError("Invalid atom serial identifiers") from exc
        if serial in all_serials:
            raise ValueError("Duplicate atom serial identifiers")
        all_serials.add(serial)
        key = residue_key(line)
        if residue is not None and key != tuple(str(v).strip() for v in residue):
            continue
        if key in removal:
            continue
        coordinates(line)
        groups.setdefault(key, []).append((index, line))
    selected_entries = []
    selections = {}
    for key, atoms in groups.items():
        locations = {line[16:17].strip() for _, line in atoms} - {""}
        chosen = preferred_altloc if preferred_altloc in locations else None
        if locations and chosen is None:
            if len(locations) != 1:
                raise ValueError(f"Ambiguous alternate conformer for residue {key}: {sorted(locations)}")
            chosen = next(iter(locations))
        picked = [(index, line) for index, line in atoms if not line[16:17].strip() or line[16:17].strip() == chosen]
        names = [line[12:16].strip() for _, line in picked]
        all_names = {line[12:16].strip() for _, line in atoms}
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate atom names in selected residue {key}")
        if set(names) != all_names:
            raise ValueError(f"Incomplete alternate conformer {chosen} for residue {key}")
        selected_entries.extend((index, line[:16] + " " + line[17:], line) for index, line in picked)
        selections[":".join(key)] = chosen or "shared"
    if not selected_entries:
        raise ValueError("Selection contains no atoms")

    selected_entries.sort(key=lambda entry: entry[0])
    selected = [entry[1] for entry in selected_entries]
    selected_serials = {int(line[6:11]) for line in selected}
    selected_residues = {residue_key(line) for line in selected}
    selected_atom_variants = {}
    for _, _, source_line in selected_entries:
        identity = _atom_identity(source_line)
        selected_atom_variants.setdefault(identity, set()).add(source_line[16:17].strip())

    retained_ter = {}
    # TER is a polymer-chain delimiter.  A ligand-only extraction can share a
    # chain identifier with a protein and must not inherit that chain's TER.
    selected_polymer_chains = {
        residue_key(source_line)[0]
        for _, _, source_line in selected_entries
        if source_line.startswith("ATOM  ")
    }
    for index, line in enumerate(lines):
        if not line.startswith("TER"):
            continue
        key = residue_key(line)
        if key[0] not in selected_polymer_chains:
            continue
        if key not in selected_residues:
            raise ValueError(f"TER record crosses the selected residue boundary: {key}")
        retained_ter[index] = line

    retained_annotations = []
    for index, line in enumerate(lines):
        if line.startswith("LINK"):
            endpoint1 = _link_endpoint(line, 12)
            endpoint2 = _link_endpoint(line, 42)
            endpoints = (endpoint1, endpoint2)
            selected_flags = [endpoint[0] in selected_residues for endpoint in endpoints]
            if not any(selected_flags):
                continue
            if not all(selected_flags):
                raise ValueError(f"LINK record crosses the selected residue/atom selection boundary: {line[:80]!r}")
            if _has_non_identity_symmetry(line):
                raise ValueError(f"LINK record uses an unsupported non-identity symmetry operation: {line[:80]!r}")
            for residue_id, atom_name, altloc in endpoints:
                identity = (residue_id, atom_name)
                variants = selected_atom_variants.get(identity)
                if not variants:
                    raise ValueError(f"LINK record references an atom removed by selection: {line[:80]!r}")
                if altloc and altloc not in variants:
                    raise ValueError(f"LINK record references an alternate location removed by selection: {line[:80]!r}")
            retained_annotations.append((index, _normalize_link_altlocs(line)))
        elif line.startswith("SSBOND"):
            endpoint1 = _ssbond_endpoint(line, 11)
            endpoint2 = _ssbond_endpoint(line, 25)
            endpoints = (endpoint1, endpoint2)
            selected_flags = [endpoint in selected_residues for endpoint in endpoints]
            if not any(selected_flags):
                continue
            if not all(selected_flags):
                raise ValueError(f"SSBOND record crosses the selected residue boundary: {line[:80]!r}")
            if _has_non_identity_symmetry(line):
                raise ValueError(f"SSBOND record uses an unsupported non-identity symmetry operation: {line[:80]!r}")
            for residue_id in endpoints:
                if (residue_id, "SG") not in selected_atom_variants:
                    raise ValueError(f"SSBOND record references an SG atom removed by selection: {line[:80]!r}")
            retained_annotations.append((index, line))

    conect = []
    for index, line in enumerate(lines):
        if not line.startswith("CONECT"):
            continue
        values = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
        if not values:
            continue
        selected_flags = [value in selected_serials for value in values]
        if any(selected_flags) and not all(selected_flags):
            raise ValueError(f"CONECT record crosses the selected atom selection boundary: {line[:80]!r}")
        if all(selected_flags):
            # Repetitions can encode bond order and are intentionally preserved.
            conect.append((index, line))

    output = []
    selected_by_index = {index: normalized for index, normalized, _ in selected_entries}
    annotations_by_index = {index: line for index, line in retained_annotations}
    conect_by_index = {index: line for index, line in conect}
    for index, line in enumerate(lines):
        if index in selected_by_index:
            output.append(selected_by_index[index])
        if index in retained_ter:
            output.append(retained_ter[index])
        if index in annotations_by_index:
            output.append(annotations_by_index[index])
        if index in conect_by_index:
            output.append(conect_by_index[index])
    return output + ["END"], selections


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
    source_sha256 = hashlib.sha256(raw).hexdigest()
    atom_identity_payload = [
        {
            "serial": int(line[6:11]),
            "atom": line[12:16].strip(),
            "residue": residue_key(line),
            "altloc": line[16:17].strip(),
            "element": line[76:78].strip(),
            "coordinates": line[30:54],
        }
        for line in original
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    frame_atom_identity_digest = hashlib.sha256(
        json.dumps(atom_identity_payload, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()
    metadata = {
        "source_file": str(source.resolve()), "source_sha256": hashlib.sha256(raw).hexdigest(),
        "destination_file": str(destination.resolve()), "atom_count": len(atoms),
        "removed_altloc_atoms": eligible_count - len(atoms),
        "removed_atom_count": sum(line.startswith(("ATOM  ", "HETATM")) for line in original) - len(atoms),
        "alternate_locations": selections,
        "chemical_identity_status": "not_evaluated_pdb_is_not_authoritative_chemistry",
        "reference_frame_id": "source-pdb-v1:" + source_sha256,
        "frame_source_sha256": source_sha256,
        "frame_atom_identity_digest": frame_atom_identity_digest,
        "retained_connectivity_records": {
            "TER": sum(1 for line in output if line.startswith("TER")),
            "LINK": sum(1 for line in output if line.startswith("LINK")),
            "SSBOND": sum(1 for line in output if line.startswith("SSBOND")),
            "CONECT": sum(1 for line in output if line.startswith("CONECT")),
        },
        "atom_map": [{"source_serial": int(line[6:11]), "output_serial": int(line[6:11]), "atom_name": line[12:16].strip(), "residue": residue_key(line)} for line in atoms],
    }
    destination.with_suffix(destination.suffix + ".preparation.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata
