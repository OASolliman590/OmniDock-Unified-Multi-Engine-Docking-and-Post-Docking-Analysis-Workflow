from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

from ..project_layout import (
    LAYOUT_DOCKING_LEGACY,
    detect_layout_profile,
    ensure_project_layout,
    pair_intent_path,
    pairlist_path,
)
from .excel_sites import load_site_catalog_from_summary
from .project_aliases import (
    build_ligand_alias_lookup,
    build_protein_alias_lookup,
    ensure_project_alias_files,
)


RECEPTOR_SUFFIXES = {".pdbqt", ".pdb"}
LIGAND_SUFFIXES = {".pdbqt", ".pdb", ".sdf", ".mol2"}

PAIR_INTENT_COLUMNS = [
    "pdb_id",
    "receptor",
    "protein_display_name",
    "ligand",
    "ligand_display_name",
    "site_id",
    "pair_source",
    "is_cocrystal_benchmark",
    "cocrystal_ligand_name",
    "cocrystal_ligand_display_name",
    "selection_mode",
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
]

PAIRLIST_COLUMNS = [
    "receptor",
    "site_id",
    "ligand",
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
    "protein_display_name",
    "ligand_display_name",
]


def _normalize_key(value: str) -> str:
    raw = Path(str(value).strip()).stem.lower()
    for suffix in ("_cleaned", "_prepared", "_prep"):
        if raw.endswith(suffix):
            raw = raw[: -len(suffix)]
    return re.sub(r"[^a-z0-9]+", "", raw)


def _extract_pdb_id(value: str) -> str:
    match = re.search(r"([0-9][A-Za-z0-9]{3})", Path(str(value)).stem)
    if not match:
        raise ValueError(f"Could not extract PDB ID from receptor name: {value}")
    return match.group(1).upper()


def _index_assets(directory: Path, allowed_suffixes: set[str]) -> Dict[str, Path]:
    directory = Path(directory)
    if not directory.exists():
        raise FileNotFoundError(f"Asset directory not found: {directory}")

    index: Dict[str, Path] = {}
    for path in sorted(directory.iterdir()):
        if not path.is_file() or path.suffix.lower() not in allowed_suffixes:
            continue
        for token in {path.name.lower(), path.stem.lower(), _normalize_key(path.name), _normalize_key(path.stem)}:
            index.setdefault(token, path)
        if "_ligand_" in path.stem.lower():
            ligand_suffix = path.stem.lower().split("_ligand_", 1)[1]
            index.setdefault(ligand_suffix, path)
            index.setdefault(_normalize_key(ligand_suffix), path)
        try:
            index.setdefault(_extract_pdb_id(path.name).lower(), path)
        except ValueError:
            pass
    if not index:
        raise ValueError(f"No usable assets found in {directory}")
    return index


def _list_unique_assets(directory: Path, allowed_suffixes: set[str]) -> List[Path]:
    return sorted(
        path
        for path in Path(directory).iterdir()
        if path.is_file() and path.suffix.lower() in allowed_suffixes
    )


def _resolve_asset(raw_name: str, asset_index: Dict[str, Path], asset_type: str) -> Path:
    candidates = [
        str(raw_name).strip().lower(),
        Path(str(raw_name).strip()).stem.lower(),
        _normalize_key(raw_name),
    ]
    if asset_type == "receptor":
        try:
            candidates.append(_extract_pdb_id(raw_name).lower())
        except ValueError:
            pass
    for candidate in candidates:
        if candidate in asset_index:
            return asset_index[candidate]
    raise ValueError(f"Could not match {asset_type} '{raw_name}' to a prepared file")


def list_prepared_asset_names(directory: Path, asset_type: str) -> List[str]:
    suffixes = RECEPTOR_SUFFIXES if asset_type == "receptor" else LIGAND_SUFFIXES
    return [path.name for path in _list_unique_assets(directory, suffixes)]


def _write_csv(path: Path, rows: Iterable[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_pairlists(
    project_root: Path,
    prepared_proteins: Path,
    prepared_ligands: Path,
    excel_path: Path,
    mode: str,
    default_site_id: str = "site_1",
    default_box_size: float = 20.0,
    curated_receptors: Optional[List[str]] = None,
    curated_ligands: Optional[List[str]] = None,
    curated_mapping: Optional[Dict[str, List[str]]] = None,
    layout_profile: str = LAYOUT_DOCKING_LEGACY,
) -> Dict[str, object]:
    summary = generate_pairlists(
        project_root=project_root,
        prepared_proteins=prepared_proteins,
        prepared_ligands=prepared_ligands,
        excel_path=excel_path,
        mode=mode,
        default_site_id=default_site_id,
        default_box_size=default_box_size,
        curated_receptors=curated_receptors,
        curated_ligands=curated_ligands,
        curated_mapping=curated_mapping,
        layout_profile=layout_profile,
    )
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    pairlist_file = pairlist_path(root, profile)
    intent_file = pair_intent_path(root, profile)
    _write_csv(intent_file, summary["pair_intent_rows"], PAIR_INTENT_COLUMNS)
    _write_csv(pairlist_file, summary["pairlist_rows"], PAIRLIST_COLUMNS)
    summary["pairlist_file"] = str(pairlist_file)
    summary["pair_intent_file"] = str(intent_file)
    return {key: value for key, value in summary.items() if key not in {"pairlist_rows", "pair_intent_rows"}}


def generate_pairlists(
    project_root: Path,
    prepared_proteins: Path,
    prepared_ligands: Path,
    excel_path: Path,
    mode: str,
    default_site_id: str = "site_1",
    default_box_size: float = 20.0,
    curated_receptors: Optional[List[str]] = None,
    curated_ligands: Optional[List[str]] = None,
    curated_mapping: Optional[Dict[str, List[str]]] = None,
    layout_profile: str = LAYOUT_DOCKING_LEGACY,
) -> Dict[str, object]:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    ensure_project_layout(root, profile)

    site_catalog, missing_coordinate_ids = load_site_catalog_from_summary(
        excel_path,
        allow_missing_coordinates=True,
        return_missing_coordinates=True,
    )
    alias_files = ensure_project_alias_files(
        project_root=root,
        prepared_proteins=prepared_proteins,
        prepared_ligands=prepared_ligands,
        site_catalog=site_catalog,
        layout_profile=profile,
    )
    receptor_index = _index_assets(prepared_proteins, RECEPTOR_SUFFIXES)
    ligand_index = _index_assets(prepared_ligands, LIGAND_SUFFIXES)
    all_ligands = _list_unique_assets(prepared_ligands, LIGAND_SUFFIXES)
    site_lookup = {str(row["pdb_id"]).upper(): row for _, row in site_catalog.iterrows()}
    protein_alias_lookup = build_protein_alias_lookup(root, profile)
    ligand_alias_lookup = build_ligand_alias_lookup(root, profile)

    warnings: List[str] = []
    if missing_coordinate_ids:
        warnings.append(
            "Skipped PDB entries with missing active-site coordinates: "
            + ", ".join(sorted(missing_coordinate_ids))
        )
    pair_intent_rows: List[Dict[str, object]] = []
    pairlist_rows: List[Dict[str, object]] = []
    seen: set[Tuple[str, str, str]] = set()

    def add_row(
        receptor_source: Path,
        ligand_source: Path,
        site_row: pd.Series,
        pair_source: str,
        is_cocrystal_benchmark: bool,
        selection_mode: str,
        cocrystal_ligand_name: str,
        cocrystal_ligand_display_name: str,
    ) -> None:
        key = (receptor_source.name, ligand_source.name, default_site_id)
        if key in seen:
            return
        seen.add(key)
        pdb_id = str(site_row.get("pdb_id") or _extract_pdb_id(receptor_source.name)).upper()
        protein_display_name = (
            protein_alias_lookup.get(receptor_source.name)
            or protein_alias_lookup.get(receptor_source.stem)
            or protein_alias_lookup.get(pdb_id)
            or str(site_row.get("protein_display_name") or "").strip()
            or pdb_id
        )
        ligand_display_name = (
            ligand_alias_lookup.get(ligand_source.stem)
            or ligand_alias_lookup.get(ligand_source.name)
            or ligand_source.stem
        )
        core_row = {
            "receptor": receptor_source.name,
            "ligand": ligand_source.name,
            "site_id": default_site_id,
            "center_x": float(site_row["center_x"]),
            "center_y": float(site_row["center_y"]),
            "center_z": float(site_row["center_z"]),
            "size_x": float(default_box_size),
            "size_y": float(default_box_size),
            "size_z": float(default_box_size),
            "protein_display_name": protein_display_name,
            "ligand_display_name": ligand_display_name,
        }
        pairlist_rows.append(core_row.copy())
        pair_intent_rows.append(
            {
                **core_row,
                "pdb_id": pdb_id,
                "protein_display_name": protein_display_name,
                "ligand_display_name": ligand_display_name,
                "pair_source": pair_source,
                "is_cocrystal_benchmark": bool(is_cocrystal_benchmark),
                "cocrystal_ligand_name": cocrystal_ligand_name,
                "cocrystal_ligand_display_name": cocrystal_ligand_display_name,
                "selection_mode": selection_mode,
            }
        )

    if mode in {"cocrystal_only", "cocrystal_plus_all", "cocrystal_plus_nonreference"}:
        reference_ligand_names: set[str] = set()
        if mode == "cocrystal_plus_nonreference":
            for _, candidate_site_row in site_catalog.iterrows():
                candidate_selected = str(candidate_site_row.get("selected_ligand") or "").strip()
                if not candidate_selected:
                    continue
                try:
                    reference_ligand_names.add(
                        _resolve_asset(candidate_selected, ligand_index, "ligand").name
                    )
                except ValueError:
                    continue
        for _, site_row in site_catalog.iterrows():
            pdb_id = str(site_row["pdb_id"]).upper()
            selected_ligand = str(site_row.get("selected_ligand") or "").strip()
            try:
                receptor_source = _resolve_asset(pdb_id, receptor_index, "receptor")
            except ValueError as exc:
                warnings.append(str(exc))
                continue

            resolved_cocrystal: Optional[Path] = None
            if not selected_ligand:
                warnings.append(f"No selected_ligand in Excel for PDB {pdb_id}; cocrystal benchmark row skipped.")
            else:
                try:
                    resolved_cocrystal = _resolve_asset(selected_ligand, ligand_index, "ligand")
                    add_row(
                        receptor_source,
                        resolved_cocrystal,
                        site_row,
                        pair_source="cocrystal",
                        is_cocrystal_benchmark=True,
                        selection_mode=mode,
                        cocrystal_ligand_name=selected_ligand,
                        cocrystal_ligand_display_name=(
                            ligand_alias_lookup.get(resolved_cocrystal.stem)
                            or ligand_alias_lookup.get(selected_ligand)
                            or resolved_cocrystal.stem
                        ),
                    )
                except ValueError as exc:
                    warnings.append(
                        f"{exc}; cocrystal benchmark row skipped for receptor {receptor_source.name}."
                    )
            if mode in {"cocrystal_plus_all", "cocrystal_plus_nonreference"}:
                ligand_pool = all_ligands
                pair_source = "expanded_all_ligands"
                if mode == "cocrystal_plus_nonreference":
                    ligand_pool = [
                        ligand_source
                        for ligand_source in all_ligands
                        if ligand_source.name not in reference_ligand_names
                    ]
                    pair_source = "expanded_nonreference_ligands"
                for ligand_source in ligand_pool:
                    add_row(
                        receptor_source,
                        ligand_source,
                        site_row,
                        pair_source=pair_source,
                        is_cocrystal_benchmark=bool(
                            resolved_cocrystal and ligand_source.name == resolved_cocrystal.name
                        ),
                        selection_mode=mode,
                        cocrystal_ligand_name=selected_ligand,
                        cocrystal_ligand_display_name=(
                            ligand_alias_lookup.get(resolved_cocrystal.stem)
                            if resolved_cocrystal
                            else ligand_alias_lookup.get(selected_ligand, "")
                        )
                        or "",
                    )

    elif mode == "curated_cartesian":
        if not curated_receptors or not curated_ligands:
            raise ValueError("curated_cartesian requires curated_receptors and curated_ligands")
        for receptor_name in curated_receptors:
            receptor_source = _resolve_asset(receptor_name, receptor_index, "receptor")
            pdb_id = _extract_pdb_id(receptor_source.name)
            site_row = site_lookup.get(pdb_id)
            if site_row is None:
                warnings.append(f"No site metadata found for receptor {receptor_source.name}")
                continue
            cocrystal_ligand_name = str(site_row.get("selected_ligand") or "")
            resolved_cocrystal = None
            if cocrystal_ligand_name:
                try:
                    resolved_cocrystal = _resolve_asset(cocrystal_ligand_name, ligand_index, "ligand")
                    add_row(
                        receptor_source,
                        resolved_cocrystal,
                        site_row,
                        pair_source="cocrystal",
                        is_cocrystal_benchmark=True,
                        selection_mode=mode,
                        cocrystal_ligand_name=cocrystal_ligand_name,
                        cocrystal_ligand_display_name=(
                            ligand_alias_lookup.get(resolved_cocrystal.stem)
                            or ligand_alias_lookup.get(cocrystal_ligand_name)
                            or resolved_cocrystal.stem
                        ),
                    )
                except ValueError as exc:
                    warnings.append(str(exc))
            for ligand_name in curated_ligands:
                ligand_source = _resolve_asset(ligand_name, ligand_index, "ligand")
                add_row(
                    receptor_source,
                    ligand_source,
                    site_row,
                    pair_source="curated_cartesian",
                    is_cocrystal_benchmark=bool(resolved_cocrystal and ligand_source.name == resolved_cocrystal.name),
                    selection_mode=mode,
                    cocrystal_ligand_name=cocrystal_ligand_name,
                    cocrystal_ligand_display_name=(
                        ligand_alias_lookup.get(resolved_cocrystal.stem)
                        if resolved_cocrystal
                        else ligand_alias_lookup.get(cocrystal_ligand_name, "")
                    )
                    or "",
                )

    elif mode == "curated_per_protein":
        if not curated_mapping:
            raise ValueError("curated_per_protein requires curated_mapping")
        for receptor_name, ligand_names in curated_mapping.items():
            receptor_source = _resolve_asset(receptor_name, receptor_index, "receptor")
            pdb_id = _extract_pdb_id(receptor_source.name)
            site_row = site_lookup.get(pdb_id)
            if site_row is None:
                warnings.append(f"No site metadata found for receptor {receptor_source.name}")
                continue
            cocrystal_ligand_name = str(site_row.get("selected_ligand") or "")
            resolved_cocrystal = None
            if cocrystal_ligand_name:
                try:
                    resolved_cocrystal = _resolve_asset(cocrystal_ligand_name, ligand_index, "ligand")
                    add_row(
                        receptor_source,
                        resolved_cocrystal,
                        site_row,
                        pair_source="cocrystal",
                        is_cocrystal_benchmark=True,
                        selection_mode=mode,
                        cocrystal_ligand_name=cocrystal_ligand_name,
                        cocrystal_ligand_display_name=(
                            ligand_alias_lookup.get(resolved_cocrystal.stem)
                            or ligand_alias_lookup.get(cocrystal_ligand_name)
                            or resolved_cocrystal.stem
                        ),
                    )
                except ValueError as exc:
                    warnings.append(str(exc))
                    resolved_cocrystal = None
            for ligand_name in ligand_names:
                ligand_source = _resolve_asset(ligand_name, ligand_index, "ligand")
                add_row(
                    receptor_source,
                    ligand_source,
                    site_row,
                    pair_source="curated_per_protein",
                    is_cocrystal_benchmark=bool(resolved_cocrystal and ligand_source.name == resolved_cocrystal.name),
                    selection_mode=mode,
                    cocrystal_ligand_name=cocrystal_ligand_name,
                    cocrystal_ligand_display_name=(
                        ligand_alias_lookup.get(resolved_cocrystal.stem)
                        if resolved_cocrystal
                        else ligand_alias_lookup.get(cocrystal_ligand_name, "")
                    )
                    or "",
                )
    else:
        raise ValueError(f"Unsupported pairlist mode: {mode}")

    if not pairlist_rows:
        raise ValueError("No pairlist rows were generated")

    return {
        "project_root": str(root),
        "layout_profile": profile,
        "pairlist_rows": pairlist_rows,
        "pair_intent_rows": pair_intent_rows,
        "pair_count": len(pairlist_rows),
        "receptor_count": len({row["receptor"] for row in pairlist_rows}),
        "ligand_count": len({row["ligand"] for row in pairlist_rows}),
        "warning_count": len(warnings),
        "warnings": warnings,
        "pair_mode": mode,
        "has_cocrystal_benchmark_rows": any(bool(row["is_cocrystal_benchmark"]) for row in pair_intent_rows),
        "protein_alias_file": alias_files["protein_alias_file"],
        "ligand_alias_file": alias_files["ligand_alias_file"],
    }
