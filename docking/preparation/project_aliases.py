from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd

from ..project_layout import detect_layout_profile, ensure_project_layout
from post_docking_analysis.ligand_naming import infer_ligand_display_name
from post_docking_analysis.protein_naming import extract_pdb_code, infer_display_name

PROTEIN_ALIAS_COLUMNS = ["pdb_id", "receptor", "display_name", "notes"]
LIGAND_ALIAS_COLUMNS = ["ligand", "display_name", "notes"]

RECEPTOR_SUFFIXES = {".pdbqt", ".pdb"}
LIGAND_SUFFIXES = {".pdbqt", ".pdb", ".sdf", ".mol2"}


def protein_aliases_path(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    layout = ensure_project_layout(root, profile)
    return layout["metadata"] / "protein_aliases.csv"


def ligand_aliases_path(project_root: Path, layout_profile: Optional[str] = None) -> Path:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    layout = ensure_project_layout(root, profile)
    return layout["metadata"] / "ligand_aliases.csv"


def _read_alias_frame(path: Path, columns: List[str]) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=columns)
    try:
        frame = pd.read_csv(path)
    except Exception:
        return pd.DataFrame(columns=columns)
    for column in columns:
        if column not in frame.columns:
            frame[column] = ""
    return frame[columns].copy()


def _list_assets(directory: Path, allowed_suffixes: set[str]) -> List[Path]:
    path = Path(directory).expanduser().resolve()
    if not path.exists():
        return []
    return sorted(
        child
        for child in path.iterdir()
        if child.is_file() and child.suffix.lower() in allowed_suffixes
    )


def _existing_value(frame: pd.DataFrame, key_column: str, key: str) -> str:
    if frame.empty or not key:
        return ""
    series = frame[key_column].fillna("").astype(str).str.strip()
    matches = frame.loc[series == key]
    if matches.empty:
        return ""
    return str(matches.iloc[0].get("display_name") or "").strip()


def _first_nonempty(values: Iterable[str]) -> str:
    for value in values:
        text = str(value or "").strip()
        if text:
            return text
    return ""


def _canonical_ligand_alias_key(ligand: str) -> str:
    stem = Path(str(ligand).strip()).stem
    lowered = stem.lower()
    if "_ligand_" in lowered:
        return lowered.split("_ligand_", 1)[1]
    return lowered


def ensure_project_alias_files(
    project_root: Path,
    prepared_proteins: Path,
    prepared_ligands: Path,
    site_catalog: Optional[pd.DataFrame] = None,
    layout_profile: Optional[str] = None,
) -> Dict[str, str]:
    protein_alias_file = protein_aliases_path(project_root, layout_profile)
    ligand_alias_file = ligand_aliases_path(project_root, layout_profile)

    existing_proteins = _read_alias_frame(protein_alias_file, PROTEIN_ALIAS_COLUMNS)
    existing_ligands = _read_alias_frame(ligand_alias_file, LIGAND_ALIAS_COLUMNS)

    protein_rows: List[Dict[str, str]] = []
    seen_proteins: set[tuple[str, str]] = set()
    seen_protein_pdb_ids: set[str] = set()
    protein_files = _list_assets(prepared_proteins, RECEPTOR_SUFFIXES)
    protein_catalog = site_catalog if isinstance(site_catalog, pd.DataFrame) else pd.DataFrame()

    for receptor_file in protein_files:
        receptor_name = receptor_file.name
        receptor_stem = receptor_file.stem
        pdb_id = extract_pdb_code(receptor_name) or extract_pdb_code(receptor_stem) or ""
        workbook_alias = ""
        if not protein_catalog.empty and pdb_id and "pdb_id" in protein_catalog.columns:
            matches = protein_catalog.loc[
                protein_catalog["pdb_id"].fillna("").astype(str).str.upper() == pdb_id.upper()
            ]
            if not matches.empty:
                workbook_alias = str(matches.iloc[0].get("protein_display_name") or "").strip()
        display_name = (
            _existing_value(existing_proteins, "receptor", receptor_name)
            or _existing_value(existing_proteins, "receptor", receptor_stem)
            or _existing_value(existing_proteins, "pdb_id", pdb_id)
            or workbook_alias
            or infer_display_name(receptor_name)
        )
        key = (pdb_id.upper(), receptor_name)
        if key in seen_proteins:
            continue
        seen_proteins.add(key)
        if pdb_id:
            seen_protein_pdb_ids.add(pdb_id.upper())
        protein_rows.append(
            {
                "pdb_id": pdb_id.upper(),
                "receptor": receptor_name,
                "display_name": display_name,
                "notes": "Edit display_name if needed; keep pdb_id/receptor stable.",
            }
        )

    if not protein_catalog.empty and "pdb_id" in protein_catalog.columns:
        for _, row in protein_catalog.iterrows():
            pdb_id = str(row.get("pdb_id") or "").strip().upper()
            if not pdb_id:
                continue
            if pdb_id in seen_protein_pdb_ids:
                continue
            key = (pdb_id, "")
            if key in seen_proteins:
                continue
            display_name = (
                _existing_value(existing_proteins, "pdb_id", pdb_id)
                or str(row.get("protein_display_name") or "").strip()
                or f"Protein {pdb_id}"
            )
            seen_proteins.add(key)
            protein_rows.append(
                {
                    "pdb_id": pdb_id,
                    "receptor": "",
                    "display_name": display_name,
                    "notes": "Edit display_name if needed; keep pdb_id/receptor stable.",
                }
            )

    ligand_rows: List[Dict[str, str]] = []
    seen_ligands: set[str] = set()
    ligand_files = _list_assets(prepared_ligands, LIGAND_SUFFIXES)
    for ligand_file in ligand_files:
        ligand_key = ligand_file.stem
        display_name = (
            _existing_value(existing_ligands, "ligand", ligand_key)
            or _existing_value(existing_ligands, "ligand", ligand_file.name)
            or infer_ligand_display_name(ligand_file.name)
        )
        if ligand_key in seen_ligands:
            continue
        seen_ligands.add(ligand_key)
        ligand_rows.append(
            {
                "ligand": ligand_key,
                "display_name": display_name,
                "notes": "Edit display_name if needed; keep ligand stable.",
            }
        )

    if not protein_catalog.empty and "selected_ligand" in protein_catalog.columns:
        for raw_value in protein_catalog["selected_ligand"].dropna().astype(str):
            ligand_key = raw_value.strip()
            if not ligand_key or ligand_key in seen_ligands:
                continue
            display_name = (
                _existing_value(existing_ligands, "ligand", ligand_key)
                or infer_ligand_display_name(ligand_key)
            )
            seen_ligands.add(ligand_key)
            ligand_rows.append(
                {
                    "ligand": ligand_key,
                    "display_name": display_name,
                    "notes": "Optional benchmark/raw ligand alias.",
                }
            )

    protein_df = pd.DataFrame(protein_rows, columns=PROTEIN_ALIAS_COLUMNS)
    ligand_df = pd.DataFrame(ligand_rows, columns=LIGAND_ALIAS_COLUMNS)
    protein_alias_file.parent.mkdir(parents=True, exist_ok=True)
    protein_df.to_csv(protein_alias_file, index=False)
    ligand_df.to_csv(ligand_alias_file, index=False)
    return {
        "protein_alias_file": str(protein_alias_file),
        "ligand_alias_file": str(ligand_alias_file),
    }


def load_protein_aliases(project_root: Path, layout_profile: Optional[str] = None) -> pd.DataFrame:
    return _read_alias_frame(protein_aliases_path(project_root, layout_profile), PROTEIN_ALIAS_COLUMNS)


def load_ligand_aliases(project_root: Path, layout_profile: Optional[str] = None) -> pd.DataFrame:
    return _read_alias_frame(ligand_aliases_path(project_root, layout_profile), LIGAND_ALIAS_COLUMNS)


def build_protein_alias_lookup(project_root: Path, layout_profile: Optional[str] = None) -> Dict[str, str]:
    frame = load_protein_aliases(project_root, layout_profile)
    mapping: Dict[str, str] = {}
    for _, row in frame.iterrows():
        display_name = str(row.get("display_name") or "").strip()
        receptor = str(row.get("receptor") or "").strip()
        pdb_id = str(row.get("pdb_id") or "").strip().upper()
        if display_name:
            if receptor:
                mapping[receptor] = display_name
                mapping[Path(receptor).stem] = display_name
            if pdb_id:
                mapping[pdb_id] = display_name
    return mapping


def build_ligand_alias_lookup(project_root: Path, layout_profile: Optional[str] = None) -> Dict[str, str]:
    frame = load_ligand_aliases(project_root, layout_profile)
    mapping: Dict[str, str] = {}
    for _, row in frame.iterrows():
        ligand = str(row.get("ligand") or "").strip()
        display_name = str(row.get("display_name") or "").strip()
        if not ligand or not display_name:
            continue
        mapping[ligand] = display_name
        mapping[Path(ligand).stem] = display_name
    return mapping


def prompt_project_aliases(
    project_root: Path,
    prompt_protein_aliases: bool = False,
    prompt_ligand_aliases: bool = False,
    layout_profile: Optional[str] = None,
) -> Dict[str, str]:
    protein_alias_file = protein_aliases_path(project_root, layout_profile)
    ligand_alias_file = ligand_aliases_path(project_root, layout_profile)

    if prompt_protein_aliases:
        protein_df = _read_alias_frame(protein_alias_file, PROTEIN_ALIAS_COLUMNS)
        if not protein_df.empty:
            groups: Dict[str, List[int]] = {}
            ordered_keys: List[str] = []
            for idx, row in protein_df.iterrows():
                pdb_id = str(row.get("pdb_id") or "").strip().upper()
                receptor = str(row.get("receptor") or "").strip()
                group_key = f"pdb:{pdb_id}" if pdb_id else f"receptor:{Path(receptor).stem.lower()}"
                if group_key not in groups:
                    groups[group_key] = []
                    ordered_keys.append(group_key)
                groups[group_key].append(idx)

            for group_key in ordered_keys:
                indices = groups[group_key]
                rows = [protein_df.loc[idx] for idx in indices]
                preferred_row = next(
                    (row for row in rows if str(row.get("receptor") or "").strip()),
                    rows[0],
                )
                pdb_id = str(preferred_row.get("pdb_id") or "").strip().upper()
                receptor = str(preferred_row.get("receptor") or "").strip()
                label = f"PDB {pdb_id}" if pdb_id else (receptor or "protein")
                current = _first_nonempty(row.get("display_name") for row in rows) or (
                    f"Protein {pdb_id}" if pdb_id else receptor
                )
                try:
                    updated = input(f"  Protein alias for {label} [{current}]: ").strip()
                except Exception:
                    updated = ""
                final_value = updated or current
                for idx in indices:
                    protein_df.at[idx, "display_name"] = final_value
            protein_df.to_csv(protein_alias_file, index=False)

    if prompt_ligand_aliases:
        ligand_df = _read_alias_frame(ligand_alias_file, LIGAND_ALIAS_COLUMNS)
        if not ligand_df.empty:
            groups: Dict[str, List[int]] = {}
            ordered_keys: List[str] = []
            for idx, row in ligand_df.iterrows():
                ligand = str(row.get("ligand") or "").strip()
                group_key = _canonical_ligand_alias_key(ligand) or ligand.lower()
                if group_key not in groups:
                    groups[group_key] = []
                    ordered_keys.append(group_key)
                groups[group_key].append(idx)

            for group_key in ordered_keys:
                indices = groups[group_key]
                rows = [ligand_df.loc[idx] for idx in indices]
                preferred_row = next(
                    (row for row in rows if "_ligand_" not in Path(str(row.get("ligand") or "")).stem.lower()),
                    rows[0],
                )
                ligand = str(preferred_row.get("ligand") or "").strip()
                label = Path(ligand).stem or ligand or "ligand"
                current = _first_nonempty(row.get("display_name") for row in rows) or infer_ligand_display_name(label)
                try:
                    updated = input(f"  Ligand alias for {label} [{current}]: ").strip()
                except Exception:
                    updated = ""
                final_value = updated or current
                for idx in indices:
                    ligand_df.at[idx, "display_name"] = final_value
            ligand_df.to_csv(ligand_alias_file, index=False)

    return {
        "protein_alias_file": str(protein_alias_file),
        "ligand_alias_file": str(ligand_alias_file),
    }
