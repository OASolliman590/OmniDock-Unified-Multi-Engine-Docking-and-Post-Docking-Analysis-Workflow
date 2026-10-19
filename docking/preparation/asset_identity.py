"""Exact prepared assets and coordinate/reference provenance shared by builders."""
from __future__ import annotations

import json
import re
from pathlib import Path

REFERENCE_FIELDS = ["reference_source", "reference_frame_id", "receptor_frame_id", "reference_pdb_id", "reference_pose_file", "reference_ligand_file"]


class AssetIndex(dict):
    ambiguous: set[str]

    def __init__(self):
        super().__init__()
        self.ambiguous = set()


def active_assets(directory: Path, suffixes: set[str]):
    paths = sorted(path for path in Path(directory).iterdir() if path.is_file() and path.suffix.lower() in suffixes)
    by_name = {path.name.lower(): path for path in paths}
    return [path for path in paths
            if f"{path.stem}_cleaned{path.suffix}".lower() not in by_name
            and (path.suffix.lower() == ".pdbqt" or f"{path.stem}.pdbqt".lower() not in by_name)]


def index_assets(directory: Path, suffixes: set[str], normalize):
    index = AssetIndex()
    for path in active_assets(directory, suffixes):
        tokens = {path.name.lower(), path.stem.lower(), normalize(path.stem)}
        if "_ligand_" in path.stem.lower():
            ligand_suffix = path.stem.lower().split("_ligand_", 1)[1]
            tokens.update({ligand_suffix, normalize(ligand_suffix)})
        else:
            match = re.search(r"([0-9][A-Za-z0-9]{3})", path.stem)
            if match:
                tokens.add(match.group(1).lower())
        for token in tokens:
            if token in index and index[token] != path:
                index.ambiguous.add(token)
            else:
                index[token] = path
    if not index:
        raise ValueError(f"No usable prepared assets in {directory}")
    return index


def sidecar(path: Path):
    metadata = path.with_suffix(path.suffix + ".preparation.json")
    if metadata.exists():
        return json.loads(metadata.read_text(encoding="utf-8"))
    return {}


def pair_reference_metadata(receptor: Path, ligand: Path):
    receptor_data, ligand_data = sidecar(receptor), sidecar(ligand)
    result = {key: str(ligand_data.get(key) or "") for key in REFERENCE_FIELDS}
    result["receptor_frame_id"] = str(receptor_data.get("receptor_frame_id") or "")
    # Cross-docking does not inherit another target's native reference pose.
    if result["reference_frame_id"] and result["receptor_frame_id"] and result["reference_frame_id"] != result["receptor_frame_id"]:
        result = {key: "" for key in REFERENCE_FIELDS}
        result["receptor_frame_id"] = str(receptor_data.get("receptor_frame_id") or "")
    return result
