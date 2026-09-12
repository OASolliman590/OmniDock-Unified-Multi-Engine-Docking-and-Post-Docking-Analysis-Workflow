"""Retrieve instance chemistry without replacing the experimental coordinates.

RCSB's chemical-component-instance endpoint is documented at
https://www.rcsb.org/docs/programmatic-access/file-download-services .
Ideal CCD coordinates and coordinate-inferred bonds are never reference poses.
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode, quote
from urllib.request import urlopen


def _native_atoms(pdb_file: Path):
    atoms = []
    identity = None
    for line in pdb_file.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        current = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip(), line[17:20].strip())
        if identity is not None and current != identity:
            raise ValueError("Reference PDB contains multiple residue instances")
        identity = current
        element = line[76:78].strip().upper()
        if not element:
            raise ValueError("Reference atom is missing an explicit element")
        if element in {"H", "D"}:
            continue
        coords = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        atoms.append((element, coords, line[12:16].strip()))
    if not atoms or identity is None:
        raise ValueError("Reference has no heavy atoms")
    if identity[2]:
        raise ValueError("RCSB instance lookup with insertion codes requires an explicit verified SDF")
    return identity, atoms


def retrieve_reference_sdf(pdb_file: Path, pdb_id: str, output: Path, *, fetch_text=None):
    """Write authoritative instance SDF only after exact native-pose matching.

    Failure leaves no new SDF; the caller records the reason. RDKit is required
    for graph validation. ``fetch_text`` permits offline endpoint regression tests.
    """
    from rdkit import Chem
    import numpy as np

    pdb_file, output = Path(pdb_file), Path(output)
    provenance_path = pdb_file.with_suffix(pdb_file.suffix + ".preparation.json")
    provenance = json.loads(provenance_path.read_text(encoding="utf-8")) if provenance_path.exists() else {}
    source = str(provenance.get("reference_source", ""))
    if source not in {"cocrystal", f"experimental:{pdb_id.upper()}"} or not provenance.get("reference_pdb_id"):
        raise ValueError("Experimental reference provenance is missing")
    identity, atoms = _native_atoms(pdb_file)
    query = urlencode({"auth_asym_id": identity[0], "auth_seq_id": identity[1], "encoding": "sdf"})
    url = f"https://models.rcsb.org/v1/{quote(pdb_id.lower(), safe='')}/ligand?{query}"
    if fetch_text is None:
        with urlopen(url, timeout=30) as response:
            text = response.read(10 * 1024 * 1024 + 1)
        if len(text) > 10 * 1024 * 1024:
            raise ValueError("Unexpectedly large RCSB ligand response")
    else:
        text = fetch_text(url)
        if isinstance(text, str):
            text = text.encode("utf-8")
    mols = list(Chem.ForwardSDMolSupplier(io.BytesIO(text), removeHs=False, sanitize=True))
    if len(mols) != 1 or mols[0] is None:
        raise ValueError("RCSB instance did not provide one valid molecular graph")
    mol = mols[0]
    if not mol.GetNumConformers():
        raise ValueError("RCSB instance lacks experimental coordinates")
    heavy = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    if len(heavy) != len(atoms):
        raise ValueError("RCSB/native reference heavy-atom counts differ")
    native = np.array([row[1] for row in atoms], dtype=float)
    if not np.isfinite(native).all():
        raise ValueError("Nonfinite experimental coordinates")
    conf = mol.GetConformer()
    used = set()
    atom_map = []
    for atom in heavy:
        coord = np.array(conf.GetAtomPosition(atom.GetIdx()), dtype=float)
        matches = [i for i, row in enumerate(atoms) if row[0] == atom.GetSymbol().upper()
                   and np.linalg.norm(coord - native[i]) <= 0.02]
        if len(matches) != 1 or matches[0] in used:
            raise ValueError("RCSB instance does not uniquely match the selected native coordinates/altloc")
        i = matches[0]
        used.add(i)
        conf.SetAtomPosition(atom.GetIdx(), tuple(native[i]))
        atom.SetProp("_TriposAtomName", atoms[i][2])
        atom_map.append({"sdf_atom_index": atom.GetIdx(), "source_atom_name": atoms[i][2]})
    frame = str(provenance.get("reference_frame_id", ""))
    if not frame:
        raise ValueError("Native coordinate frame provenance is missing")
    mol.SetProp("reference_source", f"experimental:{pdb_id.upper()}")
    mol.SetProp("reference_frame_id", frame)
    mol.SetProp("reference_pdb_id", pdb_id.upper())
    mol.SetProp("chemistry_source", url)
    output.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(output)) as writer:
        writer.write(mol)
    provenance.update({"reference_source": f"experimental:{pdb_id.upper()}",
                       "reference_pdb_id": pdb_id.upper(), "reference_pose_file": str(output.resolve()),
                       "reference_ligand_file": str(output.resolve()), "chemistry_status": "verified_instance_graph",
                       "chemistry_source": url, "sdf_atom_map": atom_map,
                       "sdf_sha256": hashlib.sha256(output.read_bytes()).hexdigest()})
    output.with_suffix(output.suffix + ".preparation.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    return provenance
