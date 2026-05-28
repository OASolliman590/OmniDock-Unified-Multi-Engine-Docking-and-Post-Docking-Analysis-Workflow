from __future__ import annotations

import csv
import json
import os
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from ..models import PairlistRow, ProjectManifest
from ..project_layout import (
    LAYOUT_DOCKING_LEGACY,
    deployment_root,
    ensure_engine_layout,
    ensure_gnina_hpc_compat_layout,
    ensure_project_layout,
    pair_curation_state_path,
    pairlist_path,
    save_manifest,
)
from .excel_sites import load_site_catalog_from_summary
from .ligand_quality import validate_prepared_ligand_pdbqt


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


def _write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _copy_if_different(source: Path, destination: Path) -> None:
    src = Path(source).expanduser().resolve()
    dest = Path(destination).expanduser().resolve()
    if src == dest:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, destination)


@dataclass
class DockingPreparationConfig:
    prepared_proteins: Path
    prepared_ligands: Path
    pair_intent: Optional[Path]
    excel_path: Optional[Path]
    output_dir: Path
    pairlist_file: Optional[Path] = None
    raw_proteins: Optional[Path] = None
    raw_ligands: Optional[Path] = None
    default_site_id: str = "site_1"
    default_box_size: float = 20.0
    asset_mode: str = "symlink"
    engines: Optional[List[str]] = None
    project_name: str = ""
    pair_mode: str = "manual"
    dock_all_proteins_to_all_ligands: bool = False
    layout_profile: str = LAYOUT_DOCKING_LEGACY


class DockingProjectBuilder:
    RECEPTOR_SUFFIXES = {".pdbqt", ".pdb"}
    LIGAND_SUFFIXES = {".pdbqt", ".pdb", ".sdf", ".mol2"}

    def __init__(self, config: DockingPreparationConfig):
        self.config = config

    def build(self) -> Dict[str, object]:
        layout = ensure_project_layout(self.config.output_dir, self.config.layout_profile)
        site_catalog = pd.DataFrame(
            columns=["pdb_id", "selected_ligand", "protein_display_name", "center_x", "center_y", "center_z"]
        )
        missing_coordinate_ids: List[str] = []
        if self.config.excel_path:
            excel_path = Path(self.config.excel_path).expanduser().resolve()
            if excel_path.exists():
                site_catalog, missing_coordinate_ids = load_site_catalog_from_summary(
                    excel_path,
                    allow_missing_coordinates=True,
                    return_missing_coordinates=True,
                )
            elif not self.config.pairlist_file:
                raise FileNotFoundError(f"Site workbook not found: {excel_path}")
        elif not self.config.pairlist_file:
            raise ValueError("excel_path is required unless pairlist_file is provided")
        site_catalog_path = layout["metadata"] / "protein_site_catalog.csv"
        site_catalog.to_csv(site_catalog_path, index=False)
        warnings: List[str] = []
        if missing_coordinate_ids:
            warnings.append(
                "Skipped PDB entries with missing active-site coordinates in site catalog: "
                + ", ".join(sorted(missing_coordinate_ids))
            )
        if self.config.pairlist_file and site_catalog.empty:
            warnings.append(
                "No Excel site catalog provided; pairlist-file coordinates were used as the source of docking boxes."
            )

        receptors = self._index_assets(self.config.prepared_proteins, self.RECEPTOR_SUFFIXES)
        ligands = self._index_assets(self.config.prepared_ligands, self.LIGAND_SUFFIXES)
        if self.config.pairlist_file:
            rows, receptor_links, ligand_links, row_warnings = self._load_existing_pairlist_rows(receptors, ligands)
        else:
            rows, receptor_links, ligand_links, row_warnings = self._build_pairlist_rows(
                site_catalog,
                receptors,
                ligands,
            )
        warnings.extend(row_warnings)
        if not rows:
            raise ValueError("No valid pairlist rows were generated")

        for source, dest_dir in receptor_links:
            self._materialize_asset(source, layout["receptors"] / source.name)
        for source, dest_dir in ligand_links:
            self._materialize_asset(source, layout["ligands"] / source.name)

        engines = self.config.engines or ["gnina", "vina", "smina", "autodock4"]
        if "autodock4" in engines:
            self._materialize_engine_specific_ligands(layout["docking_root"], layout["ligands"], "autodock4")
            self._materialize_engine_specific_receptors(layout["docking_root"], layout["receptors"], "autodock4")

        ligand_validation_issues: List[Dict[str, object]] = []
        seen_ligands = set()
        for source, _dest in ligand_links:
            candidate = layout["ligands"] / source.name
            if candidate.name in seen_ligands:
                continue
            seen_ligands.add(candidate.name)
            issues = validate_prepared_ligand_pdbqt(candidate)
            if issues:
                ligand_validation_issues.extend(issue.to_dict() for issue in issues)

        if ligand_validation_issues:
            validation_report = layout["metadata"] / "ligand_materialization_validation_report.json"
            validation_report.write_text(
                json.dumps(
                    {
                        "project_root": str(self.config.output_dir),
                        "ligands_dir": str(layout["ligands"]),
                        "issue_count": len(ligand_validation_issues),
                        "issues": ligand_validation_issues,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            raise ValueError(
                "Invalid AutoDock ligand outputs were materialized into the active docking set. "
                f"See {validation_report}"
            )

        pairlist_file = pairlist_path(self.config.output_dir, self.config.layout_profile)
        _write_csv(
            pairlist_file,
            [row.to_dict() for row in rows],
            ["receptor", "site_id", "ligand", "center_x", "center_y", "center_z", "size_x", "size_y", "size_z"],
        )

        if self.config.excel_path:
            excel_path = Path(self.config.excel_path).expanduser().resolve()
            if excel_path.exists():
                excel_snapshot = layout["metadata"] / excel_path.name
                _copy_if_different(excel_path, excel_snapshot)
        if self.config.pair_intent:
            intent_snapshot = layout["metadata"] / self.config.pair_intent.name
            _copy_if_different(self.config.pair_intent, intent_snapshot)

        for engine in engines:
            ensure_engine_layout(self.config.output_dir, engine)
        if "gnina" in engines:
            ensure_gnina_hpc_compat_layout(self.config.output_dir)

        manifest = ProjectManifest(
            project_name=self.config.project_name or self.config.output_dir.name,
            project_root=str(self.config.output_dir),
            created_at=datetime.now(timezone.utc).isoformat(),
            asset_mode=self.config.asset_mode,
            engines=engines,
            pairlist_file=str(pairlist_file),
            receptors_dir=str(layout["receptors"]),
            ligands_dir=str(layout["ligands"]),
            metadata_dir=str(layout["metadata"]),
            layout_profile=self.config.layout_profile,
            docking_root=str(layout["docking_root"]),
            post_docking_root=str(layout["post_docking_root"]),
            raw_proteins_dir=str(layout["raw_proteins"]),
            raw_ligands_dir=str(layout["raw_ligands"]),
            raw_ligands_sdf_dir=str(layout["raw_ligands_sdf"]),
            prepared_proteins_dir=str(layout["prepared_proteins"]),
            prepared_ligands_dir=str(layout["prepared_ligands"]),
            enabled_panels=engines,
            pairlist_mode=self.config.pair_mode,
            has_cocrystal_benchmark_rows=bool(self.config.pairlist_file and self.config.pair_intent),
            source_paths={
                "prepared_proteins": str(self.config.prepared_proteins),
                "prepared_ligands": str(self.config.prepared_ligands),
                "pair_intent": str(self.config.pair_intent) if self.config.pair_intent else "",
                "pairlist_file": str(self.config.pairlist_file) if self.config.pairlist_file else "",
                "excel_path": str(self.config.excel_path) if self.config.excel_path else "",
                "raw_proteins": str(self.config.raw_proteins) if self.config.raw_proteins else "",
                "raw_ligands": str(self.config.raw_ligands) if self.config.raw_ligands else "",
            },
            engine_settings={engine: {} for engine in engines},
            pair_curation_state_file=str(pair_curation_state_path(self.config.output_dir, self.config.layout_profile)),
            deployment_root=str(deployment_root(self.config.output_dir, self.config.layout_profile)),
            notes=[
                f"pair_mode={self.config.pair_mode}",
                f"dock_all_proteins_to_all_ligands={self.config.dock_all_proteins_to_all_ligands}",
                f"default site_id={self.config.default_site_id}",
                f"default box size={self.config.default_box_size}",
            ],
        )
        if "gnina" in engines:
            manifest.notes.append("compatibility_profile=gnina_hpc")
        if "autodock4" in engines:
            manifest.notes.append("autodock4_assets=isolated")
        if warnings:
            warnings_file = layout["metadata"] / "preparation_warnings.txt"
            warnings_file.write_text("\n".join(warnings), encoding="utf-8")
            manifest.notes.append(f"warnings_file={warnings_file}")
        save_manifest(self.config.output_dir, manifest.to_dict())

        return {
            "project_root": str(self.config.output_dir),
            "docking_root": str(layout["docking_root"]),
            "pair_count": len(rows),
            "receptor_count": len({row.receptor for row in rows}),
            "ligand_count": len({row.ligand for row in rows}),
            "engines": engines,
            "pairlist_file": str(pairlist_file),
            "site_catalog_file": str(site_catalog_path),
            "warning_count": len(warnings),
            "pair_mode": self.config.pair_mode,
            "layout_profile": self.config.layout_profile,
            "autodock4_isolated_assets": "autodock4" in engines,
        }

    def _index_assets(self, directory: Path, allowed_suffixes: set[str]) -> Dict[str, Path]:
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
            suffix_list = ", ".join(sorted(allowed_suffixes))
            raise ValueError(
                f"No usable assets found in {directory}. "
                f"Expected at least one file with suffix: {suffix_list}."
            )
        return index

    def _list_unique_assets(self, directory: Path, allowed_suffixes: set[str]) -> List[Path]:
        directory = Path(directory)
        return sorted(
            path for path in directory.iterdir()
            if path.is_file() and path.suffix.lower() in allowed_suffixes
        )

    def _resolve_asset(self, raw_name: str, asset_index: Dict[str, Path], asset_type: str) -> Path:
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

    def _build_pairlist_rows(
        self,
        site_catalog: pd.DataFrame,
        receptor_index: Dict[str, Path],
        ligand_index: Dict[str, Path],
    ) -> Tuple[List[PairlistRow], List[Tuple[Path, Path]], List[Tuple[Path, Path]], List[str]]:
        if self.config.pair_mode == "protein-based":
            return self._build_protein_based_rows(site_catalog, receptor_index, ligand_index)
        if self.config.pair_mode != "manual":
            raise ValueError(f"Unsupported pair_mode: {self.config.pair_mode}")
        if not self.config.pair_intent:
            raise ValueError("--pair-intent is required when pair_mode=manual")

        intent_df = pd.read_csv(self.config.pair_intent)
        if not {"receptor", "ligand"}.issubset(intent_df.columns):
            raise ValueError("Pair intent CSV must contain at least 'receptor' and 'ligand' columns")

        site_lookup = {
            str(row["pdb_id"]).upper(): row
            for _, row in site_catalog.iterrows()
        }
        rows: List[PairlistRow] = []
        receptor_links: List[Tuple[Path, Path]] = []
        ligand_links: List[Tuple[Path, Path]] = []
        warnings: List[str] = []

        for _, raw_row in intent_df.iterrows():
            try:
                receptor_source = self._resolve_asset(str(raw_row["receptor"]), receptor_index, "receptor")
                ligand_source = self._resolve_asset(str(raw_row["ligand"]), ligand_index, "ligand")
            except ValueError as exc:
                warnings.append(str(exc))
                continue
            pdb_id = _extract_pdb_id(receptor_source.name)
            raw_has_centers = all(
                pd.notna(raw_row.get(axis))
                for axis in ("center_x", "center_y", "center_z")
            )
            site_row = site_lookup.get(pdb_id)
            if site_row is None and not raw_has_centers:
                raise ValueError(
                    f"No site metadata found in Excel for receptor {receptor_source.name} ({pdb_id}) and pair intent row has no explicit center coordinates"
                )

            site_id = str(raw_row.get("site_id") or self.config.default_site_id)
            size_value = float(raw_row.get("size_x") or self.config.default_box_size)
            row = PairlistRow(
                receptor=receptor_source.name,
                site_id=site_id,
                ligand=ligand_source.name,
                center_x=float(raw_row.get("center_x") if pd.notna(raw_row.get("center_x")) else site_row["center_x"]),
                center_y=float(raw_row.get("center_y") if pd.notna(raw_row.get("center_y")) else site_row["center_y"]),
                center_z=float(raw_row.get("center_z") if pd.notna(raw_row.get("center_z")) else site_row["center_z"]),
                size_x=float(raw_row.get("size_x") if pd.notna(raw_row.get("size_x")) else size_value),
                size_y=float(raw_row.get("size_y") if pd.notna(raw_row.get("size_y")) else size_value),
                size_z=float(raw_row.get("size_z") if pd.notna(raw_row.get("size_z")) else size_value),
            )
            rows.append(row)
            receptor_links.append((receptor_source, Path("receptors")))
            ligand_links.append((ligand_source, Path("ligands")))

        return rows, receptor_links, ligand_links, warnings

    def _load_existing_pairlist_rows(
        self,
        receptor_index: Dict[str, Path],
        ligand_index: Dict[str, Path],
    ) -> Tuple[List[PairlistRow], List[Tuple[Path, Path]], List[Tuple[Path, Path]], List[str]]:
        pairlist_df = pd.read_csv(self.config.pairlist_file)
        missing = [column for column in ["receptor", "site_id", "ligand", "center_x", "center_y", "center_z", "size_x", "size_y", "size_z"] if column not in pairlist_df.columns]
        if missing:
            raise ValueError(f"Pairlist file is missing required columns: {', '.join(missing)}")

        rows: List[PairlistRow] = []
        receptor_links: List[Tuple[Path, Path]] = []
        ligand_links: List[Tuple[Path, Path]] = []
        warnings: List[str] = []
        for _, raw_row in pairlist_df.iterrows():
            try:
                receptor_source = self._resolve_asset(str(raw_row["receptor"]), receptor_index, "receptor")
                ligand_source = self._resolve_asset(str(raw_row["ligand"]), ligand_index, "ligand")
            except ValueError as exc:
                warnings.append(str(exc))
                continue

            rows.append(
                PairlistRow(
                    receptor=receptor_source.name,
                    site_id=str(raw_row["site_id"]),
                    ligand=ligand_source.name,
                    center_x=float(raw_row["center_x"]),
                    center_y=float(raw_row["center_y"]),
                    center_z=float(raw_row["center_z"]),
                    size_x=float(raw_row["size_x"]),
                    size_y=float(raw_row["size_y"]),
                    size_z=float(raw_row["size_z"]),
                )
            )
            receptor_links.append((receptor_source, Path("receptors")))
            ligand_links.append((ligand_source, Path("ligands")))
        return rows, receptor_links, ligand_links, warnings

    def _build_protein_based_rows(
        self,
        site_catalog: pd.DataFrame,
        receptor_index: Dict[str, Path],
        ligand_index: Dict[str, Path],
    ) -> Tuple[List[PairlistRow], List[Tuple[Path, Path]], List[Tuple[Path, Path]], List[str]]:
        rows: List[PairlistRow] = []
        receptor_links: List[Tuple[Path, Path]] = []
        ligand_links: List[Tuple[Path, Path]] = []
        warnings: List[str] = []
        all_ligands = self._list_unique_assets(self.config.prepared_ligands, self.LIGAND_SUFFIXES)

        for _, site_row in site_catalog.iterrows():
            pdb_id = str(site_row["pdb_id"]).upper()
            try:
                receptor_source = self._resolve_asset(pdb_id, receptor_index, "receptor")
            except ValueError as exc:
                warnings.append(str(exc))
                continue

            ligand_sources: List[Path] = []
            if self.config.dock_all_proteins_to_all_ligands:
                ligand_sources = all_ligands
            else:
                selected_ligand = str(site_row.get("selected_ligand") or "").strip()
                if not selected_ligand:
                    warnings.append(f"No selected_ligand in Excel for receptor {receptor_source.name} ({pdb_id})")
                    continue
                try:
                    ligand_sources = [self._resolve_asset(selected_ligand, ligand_index, "ligand")]
                except ValueError as exc:
                    warnings.append(str(exc))
                    continue

            for ligand_source in ligand_sources:
                row = PairlistRow(
                    receptor=receptor_source.name,
                    site_id=self.config.default_site_id,
                    ligand=ligand_source.name,
                    center_x=float(site_row["center_x"]),
                    center_y=float(site_row["center_y"]),
                    center_z=float(site_row["center_z"]),
                    size_x=float(self.config.default_box_size),
                    size_y=float(self.config.default_box_size),
                    size_z=float(self.config.default_box_size),
                )
                rows.append(row)
                receptor_links.append((receptor_source, Path("receptors")))
                ligand_links.append((ligand_source, Path("ligands")))

        return rows, receptor_links, ligand_links, warnings

    def _materialize_asset(self, source: Path, destination: Path) -> None:
        if destination.exists():
            return
        destination.parent.mkdir(parents=True, exist_ok=True)
        if self.config.asset_mode == "copy":
            shutil.copy2(source, destination)
            return
        try:
            os.symlink(source, destination)
        except OSError:
            shutil.copy2(source, destination)

    def _materialize_engine_specific_ligands(self, docking_root: Path, shared_ligands_dir: Path, engine: str) -> None:
        target_dir = docking_root / "ligands_engine_specific" / engine
        target_dir.mkdir(parents=True, exist_ok=True)
        for ligand_file in sorted(shared_ligands_dir.iterdir()):
            if not ligand_file.is_file():
                continue
            destination = target_dir / ligand_file.name
            if destination.exists():
                continue
            if self.config.asset_mode == "copy":
                shutil.copy2(ligand_file, destination)
                continue
            try:
                os.symlink(ligand_file, destination)
            except OSError:
                shutil.copy2(ligand_file, destination)

    def _materialize_engine_specific_receptors(self, docking_root: Path, shared_receptors_dir: Path, engine: str) -> None:
        target_dir = docking_root / "receptors_engine_specific" / engine
        target_dir.mkdir(parents=True, exist_ok=True)
        for receptor_file in sorted(shared_receptors_dir.iterdir()):
            if not receptor_file.is_file():
                continue
            destination = target_dir / receptor_file.name
            if destination.exists():
                continue
            if self.config.asset_mode == "copy":
                shutil.copy2(receptor_file, destination)
                continue
            try:
                os.symlink(receptor_file, destination)
            except OSError:
                shutil.copy2(receptor_file, destination)
