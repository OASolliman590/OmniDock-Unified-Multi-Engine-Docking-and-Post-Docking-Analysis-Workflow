"""
Helpers for optional protein-ligand complex filtering in post-docking analysis.

Query syntax
------------
- Clauses are separated by ``;``.
- Values inside a clause are separated by ``,``
- Supported include fields: ``protein``, ``ligand``, ``site``, ``tag``
- Supported exclude fields: ``exclude_protein``, ``exclude_ligand``,
  ``exclude_site``, ``exclude_tag``

Examples
--------
- ``protein=2FVD;ligand=Sorafenib,MS3``
- ``protein=5GRN;exclude_ligand=Acetazolamide``
- ``tag=1M17_cleaned.pdbqt_site_1_ML1.pdbqt``
- ``0`` or ``skip`` to disable filtering
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple
import re

import pandas as pd


FIELD_ALIASES: Dict[str, Tuple[str, ...]] = {
    "protein": ("protein", "receptor", "protein_label", "protein_display_name", "receptor_name"),
    "ligand": ("ligand", "ligand_name", "ligand_display_name", "ligand_label"),
    "site": ("site_id",),
    "tag": ("tag", "complex_name"),
}

QUERY_HELP_TEXT = (
    "Optional complex filter query. Examples: "
    "'protein=2FVD;ligand=Sorafenib,MS3' or "
    "'protein=5GRN;exclude_ligand=Acetazolamide'. "
    "Supported fields: protein, ligand, site, tag, and exclude_* variants. "
    "Use 0/skip/all to disable filtering."
)

DISABLED_QUERY_TOKENS = {"0", "all", "none", "skip", "no", "n", "false", "f", "empty", "-"}


def parse_complex_query(query: str) -> List[Tuple[bool, str, List[str]]]:
    clauses: List[Tuple[bool, str, List[str]]] = []
    text = str(query or "").strip()
    if not text:
        return clauses
    if text.lower() in DISABLED_QUERY_TOKENS:
        return clauses

    for raw_clause in re.split(r"[;\n]+", text):
        clause = str(raw_clause).strip()
        if not clause:
            continue
        if "=" not in clause:
            raise ValueError(
                f"Invalid complex-query clause '{clause}'. Expected field=value syntax."
            )
        raw_field, raw_values = clause.split("=", 1)
        field = str(raw_field).strip().lower()
        exclude = False
        if field.startswith("exclude_"):
            exclude = True
            field = field[len("exclude_") :]
        elif field.startswith("not_"):
            exclude = True
            field = field[len("not_") :]
        if field not in FIELD_ALIASES:
            raise ValueError(
                f"Unsupported complex-query field '{raw_field}'. "
                f"Supported fields: {', '.join(sorted(FIELD_ALIASES))} and exclude_* variants."
            )
        values = [value.strip() for value in str(raw_values).split(",") if value.strip()]
        if not values:
            raise ValueError(f"Complex-query clause '{clause}' does not contain any values.")
        clauses.append((exclude, field, values))
    return clauses


def build_tag_series(frame: pd.DataFrame) -> pd.Series:
    if frame is None or frame.empty:
        return pd.Series(dtype="object")
    if "tag" in frame.columns:
        return frame["tag"].astype(str)
    if {"receptor", "site_id", "ligand"}.issubset(frame.columns):
        return frame.apply(
            lambda row: f"{row['receptor']}_{row['site_id']}_{row['ligand']}",
            axis=1,
        ).astype(str)
    if "complex_name" in frame.columns:
        return frame["complex_name"].astype(str)
    return pd.Series([""] * len(frame), index=frame.index, dtype="object")


def tags_from_frame(frame: pd.DataFrame) -> Set[str]:
    tag_series = build_tag_series(frame)
    if tag_series.empty:
        return set()
    return {str(value).strip() for value in tag_series.dropna().tolist() if str(value).strip()}


def filter_frame_by_complex_query(frame: pd.DataFrame, query: str) -> pd.DataFrame:
    if frame is None or frame.empty:
        return frame.copy() if isinstance(frame, pd.DataFrame) else pd.DataFrame()

    clauses = parse_complex_query(query)
    if not clauses:
        return frame.copy()

    filtered = frame.copy()
    working = filtered.copy()
    working["_complex_query_tag"] = build_tag_series(filtered)

    mask = pd.Series(True, index=working.index)
    for exclude, field, values in clauses:
        candidate_columns = [column for column in FIELD_ALIASES[field] if column in working.columns]
        if field == "tag":
            candidate_columns = list(dict.fromkeys(candidate_columns + ["_complex_query_tag"]))
        if not candidate_columns:
            raise ValueError(
                f"Complex-query field '{field}' is not available in this analysis context."
            )

        clause_mask = pd.Series(False, index=working.index)
        for column in candidate_columns:
            series = working[column].astype(str)
            for value in values:
                clause_mask = clause_mask | series.str.contains(re.escape(value), case=False, na=False)

        mask = mask & (~clause_mask if exclude else clause_mask)

    filtered = filtered.loc[mask].copy()
    return filtered


def normalize_result_tag(path_or_name: object) -> str:
    stem = Path(str(path_or_name)).stem
    for suffix in ("_poses", "_top", "_out", "_log"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem


def filter_paths_by_tags(paths: Sequence[Path], tags: Iterable[str]) -> List[Path]:
    allowed = {str(tag).strip() for tag in tags if str(tag).strip()}
    if not allowed:
        return list(paths)
    return [path for path in paths if normalize_result_tag(path) in allowed]
