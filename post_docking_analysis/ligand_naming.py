"""
Ligand naming utilities for reports and visualization labeling.

Provides light-weight heuristics and override support to convert ligand file
names/identifiers into stable display names suitable for figures.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Optional
import re

import pandas as pd

LIGAND_CODE_PATTERN = re.compile(r"_ligand_([A-Za-z0-9]{1,20})(?:_[A-Za-z]_\d+)?", re.IGNORECASE)

NOISE_TOKENS = {
    "ligand",
    "prepared",
    "prep",
    "pdb",
    "pdbqt",
    "sdf",
    "poses",
    "pose",
    "top",
    "out",
    "cleaned",
}


def infer_ligand_display_name(identifier: str) -> str:
    """Infer a compact ligand label from a free-form ligand identifier."""
    text = str(identifier or "").strip()
    if not text:
        return "Unknown Ligand"

    stem = Path(text).stem
    match = LIGAND_CODE_PATTERN.search(stem)
    if match:
        return match.group(1)

    tokens = [tok for tok in re.split(r"[^A-Za-z0-9]+", stem) if tok]
    kept = [tok for tok in tokens if tok.lower() not in NOISE_TOKENS]
    if not kept:
        return stem or "Unknown Ligand"

    return kept[0] if len(kept) == 1 else " ".join(kept)


def build_ligand_name_mapping(
    identifiers: Iterable[str],
    pairlist_df: Optional[pd.DataFrame] = None,
    overrides_file: Optional[Path] = None,
) -> Dict[str, str]:
    """
    Build a mapping from ligand identifiers/stems to display names.

    The mapping stores both the original identifier and its stem for robust
    lookup across file names, tags, and already-trimmed ligand strings.
    """
    mapping: Dict[str, str] = {}

    def _register(identifier: str, display_name: Optional[str] = None) -> None:
        text = str(identifier or "").strip()
        if not text:
            return
        display = str(display_name or infer_ligand_display_name(text)).strip()
        if not display:
            return
        stem = Path(text).stem
        mapping[text] = display
        mapping[stem] = display

    for identifier in identifiers:
        _register(str(identifier or ""))

    if pairlist_df is not None and not pairlist_df.empty:
        if {"ligand", "ligand_display_name"}.issubset(pairlist_df.columns):
            subset = pairlist_df[["ligand", "ligand_display_name"]].dropna()
            for _, row in subset.iterrows():
                ligand = str(row.get("ligand", "")).strip()
                display = str(row.get("ligand_display_name", "")).strip()
                if ligand and display:
                    _register(ligand, display)
        if {"cocrystal_ligand_name", "cocrystal_ligand_display_name"}.issubset(pairlist_df.columns):
            subset = pairlist_df[["cocrystal_ligand_name", "cocrystal_ligand_display_name"]].dropna()
            for _, row in subset.iterrows():
                ligand = str(row.get("cocrystal_ligand_name", "")).strip()
                display = str(row.get("cocrystal_ligand_display_name", "")).strip()
                if ligand and display:
                    _register(ligand, display)

    if overrides_file and overrides_file.exists():
        try:
            override_df = pd.read_csv(overrides_file)
            if {"ligand", "display_name"}.issubset(override_df.columns):
                for _, row in override_df.iterrows():
                    ligand = str(row.get("ligand", "")).strip()
                    display = str(row.get("display_name", "")).strip()
                    if ligand and display:
                        _register(ligand, display)
        except Exception:
            pass

    return mapping


def resolve_ligand_display_name(identifier: str, mapping: Optional[Dict[str, str]] = None) -> str:
    """Resolve a ligand display name from a mapping with heuristic fallback."""
    text = str(identifier or "").strip()
    if not text:
        return "Unknown Ligand"

    if mapping:
        if text in mapping:
            return mapping[text]
        stem = Path(text).stem
        if stem in mapping:
            return mapping[stem]

    return infer_ligand_display_name(text)


def ligand_mapping_to_dataframe(mapping: Dict[str, str]) -> pd.DataFrame:
    """Flatten a ligand mapping into an audit-friendly dataframe."""
    rows = [
        {"ligand": key, "display_name": value}
        for key, value in sorted(mapping.items())
    ]
    return pd.DataFrame(rows, columns=["ligand", "display_name"])
