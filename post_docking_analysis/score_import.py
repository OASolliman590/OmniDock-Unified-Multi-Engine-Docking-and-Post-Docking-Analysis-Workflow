"""Explicit user-supplied score tables are sources, never inferred cache files."""
import json
from pathlib import Path
import pandas as pd
import numpy as np
from post_docking_analysis.pose_geometry import content_hash


def read_explicit_score_import(scores_directory: Path):
    table = Path(scores_directory) / "normalized_scores.csv"
    marker = Path(scores_directory) / "normalized_scores.import.json"
    if not marker.is_file():
        return None
    payload = json.loads(marker.read_text(encoding="utf-8"))
    if payload.get("kind") != "imported_scores" or not table.is_file() or payload.get("sha256") != content_hash(table):
        raise ValueError("Explicit score import requires kind=imported_scores and matching CSV sha256")
    frame = pd.read_csv(table)
    required = {"engine", "tag", "protein", "ligand", "site_id", "pose", "affinity_kcal_mol"}
    if not required.issubset(frame.columns):
        raise ValueError("Imported score table is missing required identity or score columns")
    for column in ("engine", "tag", "protein", "ligand", "site_id"):
        values = frame[column].astype("string").str.strip()
        if values.isna().any() or values.eq("").any():
            raise ValueError(f"Imported score identity {column} must be nonempty")
        frame[column] = values
    frame["engine"] = frame["engine"].str.lower()
    for column in ("pose", "affinity_kcal_mol"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
        if not np.isfinite(frame[column].to_numpy(dtype=float)).all():
            raise ValueError(f"Imported score {column} must be finite")
    if ((frame["pose"] < 1) | (frame["pose"] % 1 != 0)).any():
        raise ValueError("Imported pose indices must be positive integers")
    if frame.duplicated(["engine", "tag", "pose"]).any():
        raise ValueError("Imported score table has duplicate pose identities")
    frame["score_provenance"] = "explicit_import_unverified_execution"
    return frame
