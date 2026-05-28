"""
Protein naming utilities for visualization labeling.

Provides light-weight heuristics to convert receptor/PDB-like identifiers
into stable display names suitable for figures and output manifests.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional
import re

import pandas as pd

PDB_CODE_PATTERN = re.compile(r"(?<![A-Za-z0-9])([0-9][A-Za-z0-9]{3})(?![A-Za-z0-9])")

# Tokens that are useful in file names but noisy in human labels.
NOISE_TOKENS = {
    "cleaned",
    "prep",
    "prepared",
    "receptor",
    "protein",
    "pdb",
    "pdbqt",
    "apo",
    "holo",
    "chain",
    "site",
    "series",
}


def extract_pdb_code(value: str) -> Optional[str]:
    """Extract a 4-character PDB code from a free-form identifier."""
    if not value:
        return None

    text = str(value).strip()
    match = PDB_CODE_PATTERN.search(text.upper())
    if match:
        return match.group(1).upper()

    # Fallback token scan for patterns like "1abc".
    for token in re.split(r"[^A-Za-z0-9]+", text):
        if len(token) == 4 and token[:1].isdigit() and token.isalnum():
            return token.upper()
    return None


def _to_tokens(value: str) -> List[str]:
    stem = Path(str(value)).stem
    return [t for t in re.split(r"[^A-Za-z0-9]+", stem) if t]


def infer_display_name(identifier: str) -> str:
    """
    Infer a readable protein display name from an identifier.

    Examples:
    - "Caspase3_3H0E_cleaned.pdbqt" -> "Caspase3"
    - "1M17_cleaned.pdbqt" -> "Protein 1M17"
    """
    code = extract_pdb_code(identifier)
    tokens = _to_tokens(identifier)

    kept: List[str] = []
    for token in tokens:
        lower = token.lower()
        if lower in NOISE_TOKENS:
            continue
        if code and lower == code.lower():
            continue
        kept.append(token)

    if kept:
        # Keep compact naming for files while making labels readable.
        return " ".join(t[:1].upper() + t[1:] for t in kept)
    if code:
        return f"Protein {code}"
    return Path(str(identifier)).stem


def _display_priority(name: str) -> int:
    """Prefer explicit names over generic 'Protein <PDB>' labels."""
    generic = name.lower().startswith("protein ")
    return (0 if generic else 1000) + len(name)


def build_protein_name_mapping(
    receptor_files: List[Path],
    pairlist_df: Optional[pd.DataFrame] = None,
    overrides_file: Optional[Path] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Build mapping dictionaries for receptor/protein display names.

    Returns a dict:
    - 'pdb_code_to_name': mapping from PDB code to display name
    - 'receptor_to_name': mapping from receptor stem/name to display name
    """
    pdb_code_to_name: Dict[str, str] = {}
    receptor_to_name: Dict[str, str] = {}

    def _consider(source_identifier: str) -> None:
        code = extract_pdb_code(source_identifier)
        display = infer_display_name(source_identifier)
        stem = Path(str(source_identifier)).stem

        receptor_to_name[stem] = display
        receptor_to_name[source_identifier] = display

        if code:
            current = pdb_code_to_name.get(code)
            if current is None or _display_priority(display) > _display_priority(current):
                pdb_code_to_name[code] = display

    for receptor_file in receptor_files:
        _consider(receptor_file.name)
        _consider(receptor_file.stem)

    if pairlist_df is not None and not pairlist_df.empty and "receptor" in pairlist_df.columns:
        for receptor_name in pairlist_df["receptor"].dropna().astype(str):
            _consider(receptor_name)
        if "protein_display_name" in pairlist_df.columns:
            subset = pairlist_df[["receptor", "protein_display_name"]].dropna()
            for _, row in subset.iterrows():
                receptor = str(row.get("receptor", "")).strip()
                display = str(row.get("protein_display_name", "")).strip()
                if not receptor or not display:
                    continue
                receptor_to_name[receptor] = display
                receptor_to_name[Path(receptor).stem] = display
                code = extract_pdb_code(receptor)
                if code:
                    pdb_code_to_name[code] = display
        if {"pdb_id", "protein_display_name"}.issubset(pairlist_df.columns):
            subset = pairlist_df[["pdb_id", "protein_display_name"]].dropna()
            for _, row in subset.iterrows():
                pdb_id = str(row.get("pdb_id", "")).strip().upper()
                display = str(row.get("protein_display_name", "")).strip()
                if pdb_id and display:
                    pdb_code_to_name[pdb_id] = display

    # Optional user overrides (columns: pdb_code, display_name) and/or
    # (columns: receptor, display_name).
    if overrides_file and overrides_file.exists():
        try:
            override_df = pd.read_csv(overrides_file)
            if {"pdb_code", "display_name"}.issubset(override_df.columns):
                for _, row in override_df.iterrows():
                    code = str(row.get("pdb_code", "")).strip().upper()
                    display = str(row.get("display_name", "")).strip()
                    if code and display:
                        pdb_code_to_name[code] = display
            if {"receptor", "display_name"}.issubset(override_df.columns):
                for _, row in override_df.iterrows():
                    receptor = str(row.get("receptor", "")).strip()
                    display = str(row.get("display_name", "")).strip()
                    if receptor and display:
                        receptor_to_name[receptor] = display
                        receptor_to_name[Path(receptor).stem] = display
                        code = extract_pdb_code(receptor)
                        if code:
                            pdb_code_to_name[code] = display
        except Exception:
            # Non-fatal: keep auto-generated mapping.
            pass

    return {
        "pdb_code_to_name": pdb_code_to_name,
        "receptor_to_name": receptor_to_name,
    }


def resolve_protein_display_name(identifier: str, mapping: Optional[Dict[str, Dict[str, str]]] = None) -> str:
    """
    Resolve display name from mapping with robust fallback.
    """
    if not identifier:
        return "Unknown Protein"

    if not mapping:
        return infer_display_name(identifier)

    receptor_to_name = mapping.get("receptor_to_name", {})
    pdb_code_to_name = mapping.get("pdb_code_to_name", {})

    text = str(identifier)
    stem = Path(text).stem
    if text in receptor_to_name:
        return receptor_to_name[text]
    if stem in receptor_to_name:
        return receptor_to_name[stem]

    code = extract_pdb_code(text)
    if code and code in pdb_code_to_name:
        return pdb_code_to_name[code]
    return infer_display_name(text)


def format_protein_label(identifier: str, mapping: Optional[Dict[str, Dict[str, str]]] = None) -> str:
    """
    Resolve a human label and append the detected PDB code as `(CODE)`.

    Examples
    --------
    - "1M17_cleaned.pdbqt" -> "Protein 1M17 (1M17)"
    - "VEGFR2_4AG8_cleaned.pdbqt" -> "VEGFR2 (4AG8)"
    """
    if not identifier:
        return "Unknown Protein"

    display = resolve_protein_display_name(identifier, mapping)
    code = extract_pdb_code(identifier)
    if not code:
        return display

    # Keep labels stable if code is already explicitly present.
    if re.search(rf"\({re.escape(code)}\)", str(display), flags=re.IGNORECASE):
        return str(display)
    return f"{display} ({code})"


def mapping_to_dataframe(mapping: Dict[str, Dict[str, str]]) -> pd.DataFrame:
    """Convert mapping into a flat table for reporting/auditing."""
    rows = []
    pdb_code_to_name = mapping.get("pdb_code_to_name", {})
    receptor_to_name = mapping.get("receptor_to_name", {})

    for code, display in sorted(pdb_code_to_name.items()):
        rows.append(
            {
                "record_type": "pdb_code",
                "key": code,
                "display_name": display,
            }
        )
    for receptor_key, display in sorted(receptor_to_name.items()):
        rows.append(
            {
                "record_type": "receptor",
                "key": receptor_key,
                "display_name": display,
            }
        )

    return pd.DataFrame(rows)
