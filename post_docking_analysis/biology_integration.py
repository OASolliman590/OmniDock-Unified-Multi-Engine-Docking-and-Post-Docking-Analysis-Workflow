"""
External biology table ingestion and docking-correlation utilities.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


SUPPORTED_BIOLOGY_EXTENSIONS = {".csv", ".tsv", ".txt", ".json"}


def _normalize_text(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip()


def _key_tuples_to_rows(
    key_names: Iterable[str],
    values: Iterable[Tuple[str, ...]],
    *,
    max_rows: int = 250,
) -> List[Dict[str, str]]:
    keys = [str(value) for value in key_names]
    rows: List[Dict[str, str]] = []
    for tuple_values in sorted(set(values)):
        row = {key: str(value) for key, value in zip(keys, tuple_values)}
        row["composite_key"] = "|".join(str(value) for value in tuple_values)
        rows.append(row)
        if len(rows) >= max_rows:
            break
    return rows


def load_biology_table(path: str) -> Tuple[pd.DataFrame, Dict[str, object]]:
    bio_path = Path(str(path or "")).expanduser().resolve()
    if not bio_path.exists():
        raise FileNotFoundError(f"Biology file does not exist: {bio_path}")
    if bio_path.suffix.lower() not in SUPPORTED_BIOLOGY_EXTENSIONS:
        raise ValueError(
            f"Unsupported biology file extension: {bio_path.suffix}. "
            f"Supported: {', '.join(sorted(SUPPORTED_BIOLOGY_EXTENSIONS))}"
        )

    suffix = bio_path.suffix.lower()
    if suffix == ".json":
        payload = json.loads(bio_path.read_text(encoding="utf-8"))
        if isinstance(payload, list):
            frame = pd.DataFrame(payload)
        elif isinstance(payload, dict):
            if "records" in payload and isinstance(payload["records"], list):
                frame = pd.DataFrame(payload["records"])
            else:
                frame = pd.DataFrame([payload])
        else:
            raise ValueError("JSON biology file must contain an object or a list of objects")
    elif suffix == ".csv":
        frame = pd.read_csv(bio_path)
    else:
        # Treat .tsv/.txt as tab-separated.
        frame = pd.read_csv(bio_path, sep="\t")

    frame.columns = [str(column).strip() for column in frame.columns]
    info = {
        "path": str(bio_path),
        "rows": int(len(frame)),
        "columns": [str(column) for column in frame.columns],
    }
    return frame, info


def _resolve_mapping_keys(
    docking_df: pd.DataFrame,
    biology_df: pd.DataFrame,
    preferred_mode: str = "auto",
) -> Tuple[List[str], List[str], str]:
    mode = str(preferred_mode or "auto").strip().lower()
    dock_cols = {str(column).strip().lower(): str(column) for column in docking_df.columns}
    bio_cols = {str(column).strip().lower(): str(column) for column in biology_df.columns}

    def _has(*keys: str) -> bool:
        return all(key in dock_cols and key in bio_cols for key in keys)

    if mode in {"tag", "pair_tag"} and _has("tag"):
        return [dock_cols["tag"]], [bio_cols["tag"]], "tag"
    if mode in {"protein_ligand", "protein+ligand"} and _has("protein", "ligand"):
        return [dock_cols["protein"], dock_cols["ligand"]], [bio_cols["protein"], bio_cols["ligand"]], "protein_ligand"
    if mode in {"ligand"} and _has("ligand"):
        return [dock_cols["ligand"]], [bio_cols["ligand"]], "ligand"
    if mode in {"protein"} and _has("protein"):
        return [dock_cols["protein"]], [bio_cols["protein"]], "protein"

    if _has("tag"):
        return [dock_cols["tag"]], [bio_cols["tag"]], "tag"
    if _has("protein", "ligand"):
        return [dock_cols["protein"], dock_cols["ligand"]], [bio_cols["protein"], bio_cols["ligand"]], "protein_ligand"
    if _has("ligand"):
        return [dock_cols["ligand"]], [bio_cols["ligand"]], "ligand"
    if _has("protein"):
        return [dock_cols["protein"]], [bio_cols["protein"]], "protein"
    return [], [], "none"


def attach_biology_annotations(
    docking_df: pd.DataFrame,
    biology_df: pd.DataFrame,
    mapping_mode: str = "auto",
    biology_prefix: str = "bio_",
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    if docking_df is None or docking_df.empty:
        return docking_df.copy(), {"mapping_mode": "none", "matched_rows": 0, "docking_rows": 0, "biology_rows": 0}
    if biology_df is None or biology_df.empty:
        return docking_df.copy(), {"mapping_mode": "none", "matched_rows": 0, "docking_rows": int(len(docking_df)), "biology_rows": 0}

    dock_keys, bio_keys, resolved_mode = _resolve_mapping_keys(docking_df, biology_df, preferred_mode=mapping_mode)
    if not dock_keys or not bio_keys:
        return docking_df.copy(), {
            "mapping_mode": "none",
            "matched_rows": 0,
            "docking_rows": int(len(docking_df)),
            "biology_rows": int(len(biology_df)),
            "unresolved_reason": "No shared mapping columns found (expected tag and/or protein/ligand).",
        }

    dock = docking_df.copy()
    bio = biology_df.copy()
    for column in dock_keys:
        dock[column] = _normalize_text(dock[column])
    for column in bio_keys:
        bio[column] = _normalize_text(bio[column])

    rename_map = {}
    for column in bio.columns:
        if column in bio_keys:
            continue
        rename_map[column] = f"{biology_prefix}{column}"
    bio = bio.rename(columns=rename_map)

    matched = dock.merge(
        bio,
        left_on=dock_keys,
        right_on=bio_keys,
        how="left",
        suffixes=("", "_bio"),
    )
    bio_payload_cols = [column for column in matched.columns if str(column).startswith(biology_prefix)]
    matched_rows = int((matched[bio_payload_cols].notna().any(axis=1)).sum()) if bio_payload_cols else 0

    unmatched_biology_rows = 0
    unmatched_docking_rows = 0
    unmatched_biology_examples: List[Dict[str, str]] = []
    unmatched_docking_examples: List[Dict[str, str]] = []
    if bio_keys:
        dock_key_tuples = set(tuple(values) for values in dock[dock_keys].drop_duplicates().itertuples(index=False, name=None))
        bio_key_tuples = set(tuple(values) for values in bio[bio_keys].drop_duplicates().itertuples(index=False, name=None))
        unmatched_biology = bio_key_tuples - dock_key_tuples
        unmatched_docking = dock_key_tuples - bio_key_tuples
        unmatched_biology_rows = int(len(unmatched_biology))
        unmatched_docking_rows = int(len(unmatched_docking))
        unmatched_biology_examples = _key_tuples_to_rows(bio_keys, unmatched_biology)
        unmatched_docking_examples = _key_tuples_to_rows(dock_keys, unmatched_docking)

    report = {
        "mapping_mode": resolved_mode,
        "mapped_on_docking_columns": dock_keys,
        "mapped_on_biology_columns": bio_keys,
        "docking_rows": int(len(docking_df)),
        "biology_rows": int(len(biology_df)),
        "matched_rows": matched_rows,
        "match_rate": float(matched_rows / len(docking_df)) if len(docking_df) else 0.0,
        "unmatched_biology_keys": unmatched_biology_rows,
        "unmatched_docking_keys": unmatched_docking_rows,
        "unmatched_biology_examples": unmatched_biology_examples,
        "unmatched_docking_examples": unmatched_docking_examples,
        "biology_columns_prefixed": bio_payload_cols,
    }
    return matched, report


def compute_biology_correlations(
    annotated_df: pd.DataFrame,
    score_columns: Optional[List[str]] = None,
    biology_prefix: str = "bio_",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if annotated_df is None or annotated_df.empty:
        empty_cols = ["score_column", "biology_column", "paired_rows", "spearman_corr", "pearson_corr"]
        return pd.DataFrame(columns=empty_cols), pd.DataFrame(columns=["protein"] + empty_cols)

    scores = [column for column in (score_columns or ["consensus_score", "best_affinity_kcal_mol"]) if column in annotated_df.columns]
    biology_numeric = []
    for column in annotated_df.columns:
        if not str(column).startswith(biology_prefix):
            continue
        numeric = pd.to_numeric(annotated_df[column], errors="coerce")
        if numeric.notna().sum() >= 2:
            biology_numeric.append(column)
    if not scores or not biology_numeric:
        empty_cols = ["score_column", "biology_column", "paired_rows", "spearman_corr", "pearson_corr"]
        return pd.DataFrame(columns=empty_cols), pd.DataFrame(columns=["protein"] + empty_cols)

    def _pair_corr_rows(frame: pd.DataFrame, include_protein: Optional[str] = None) -> List[Dict[str, object]]:
        rows: List[Dict[str, object]] = []
        for score_column in scores:
            score = pd.to_numeric(frame[score_column], errors="coerce")
            for bio_column in biology_numeric:
                bio = pd.to_numeric(frame[bio_column], errors="coerce")
                mask = score.notna() & bio.notna()
                paired_rows = int(mask.sum())
                if paired_rows < 3:
                    spearman = np.nan
                    pearson = np.nan
                else:
                    spearman = score[mask].corr(bio[mask], method="spearman")
                    pearson = score[mask].corr(bio[mask], method="pearson")
                row = {
                    "score_column": score_column,
                    "biology_column": bio_column,
                    "paired_rows": paired_rows,
                    "spearman_corr": float(spearman) if pd.notna(spearman) else np.nan,
                    "pearson_corr": float(pearson) if pd.notna(pearson) else np.nan,
                }
                if include_protein is not None:
                    row["protein"] = include_protein
                rows.append(row)
        return rows

    global_rows = _pair_corr_rows(annotated_df)
    per_protein_rows: List[Dict[str, object]] = []
    if "protein" in annotated_df.columns:
        for protein, group in annotated_df.groupby("protein", dropna=False):
            per_protein_rows.extend(_pair_corr_rows(group, include_protein=str(protein)))

    global_df = pd.DataFrame(
        global_rows,
        columns=["score_column", "biology_column", "paired_rows", "spearman_corr", "pearson_corr"],
    )
    per_protein_df = pd.DataFrame(
        per_protein_rows,
        columns=["protein", "score_column", "biology_column", "paired_rows", "spearman_corr", "pearson_corr"],
    )
    return global_df, per_protein_df
