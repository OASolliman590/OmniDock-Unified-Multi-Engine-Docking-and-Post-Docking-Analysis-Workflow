"""Spec 033 R1: every workflow step resolves the layout profile from the project manifest."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

import docking.cli as docking_cli
import docking.preparation.excel_sites as excel_sites
import docking.preparation.pairlist_builder as pairlist_builder
from docking.project_layout import (
    detect_layout_profile,
    ensure_project_layout,
    pairlist_path,
    resolve_layout_profile,
)
from workflow import execution
from workflow.interactive import _project_layout


SITE_ROW = {
    "pdb_id": "1IEP",
    "selected_ligand": "1IEP_ligand_STI_A_201.pdbqt",
    "center_x": 1.0,
    "center_y": 2.0,
    "center_z": 3.0,
}


@pytest.fixture
def offline_site_catalog(monkeypatch):
    """Replace the Excel site catalog with an in-memory one (no workbook parsing)."""
    catalog = pd.DataFrame([SITE_ROW])

    def _fake_catalog(*_args, **kwargs):
        if kwargs.get("return_missing_coordinates"):
            return catalog, []
        return catalog

    monkeypatch.setattr(excel_sites, "load_site_catalog_from_summary", _fake_catalog)
    monkeypatch.setattr(pairlist_builder, "load_site_catalog_from_summary", _fake_catalog)


def _prepare_assets(root: Path, layout_profile: str) -> None:
    """Fake prepared receptor/ligand files and the site workbook, in the profile's own folders."""
    layout = ensure_project_layout(root, layout_profile)
    (layout["prepared_proteins"] / "1IEP_receptor.pdbqt").write_text("ATOM\n", encoding="utf-8")
    (layout["prepared_ligands"] / "1IEP_ligand_STI_A_201.pdbqt").write_text("ATOM\n", encoding="utf-8")
    (layout["raw_proteins"] / "multi_pdb_analysis.xlsx").write_bytes(b"")  # faked catalog, content unused


def _dry_run_vina(root: Path, layout_profile: str) -> int:
    # Shared receptor/ligand assets are normally materialized by `prep project`.
    layout = ensure_project_layout(root, layout_profile)
    (layout["receptors"] / "1IEP_receptor.pdbqt").write_text("ATOM\n", encoding="utf-8")
    (layout["ligands"] / "1IEP_ligand_STI_A_201.pdbqt").write_text("ATOM\n", encoding="utf-8")
    return docking_cli.dock_main(
        [
            "--project-dir",
            str(root),
            "--engines",
            "vina",
            "--dry-run",
            "--no-ligand-qc-gate",
            "--no-receptor-qc-gate",
        ]
    )


def test_canonical_prep_pairlist_then_dry_run_uses_canonical_pairlist(tmp_path, offline_site_catalog):
    root = tmp_path / "canonical_project"
    execution.run_workflow_init(str(root), layout_profile="canonical")
    _prepare_assets(root, "canonical")

    result = execution.run_prepare_pairlist(str(root), mode="cocrystal_only")

    assert result.status == "completed"
    assert not (root / "4-Docking").exists(), "canonical prep must not create a legacy 4-Docking folder"
    assert not (root / "2-Raw_Protien").exists()
    assert detect_layout_profile(root) == "canonical"

    canonical_pairlist = root / "pairlist.csv"
    assert canonical_pairlist.is_file()
    assert pairlist_path(root) == canonical_pairlist.resolve()
    assert pairlist_path(root, "canonical") == canonical_pairlist.resolve()
    assert len(pd.read_csv(canonical_pairlist)) == 1

    # The dry run reads the same file that prep wrote; it must not raise FileNotFoundError.
    assert _dry_run_vina(root, "canonical") == 0
    assert not (root / "4-Docking").exists()


def test_legacy_prep_pairlist_writes_into_docking_folder(tmp_path, offline_site_catalog):
    root = tmp_path / "legacy_project"
    execution.run_workflow_init(str(root), layout_profile="docking_legacy")
    _prepare_assets(root, "docking_legacy")

    result = execution.run_prepare_pairlist(str(root), mode="cocrystal_only")

    assert result.status == "completed"
    legacy_pairlist = root / "4-Docking" / "pairlist.csv"
    assert legacy_pairlist.is_file()
    assert not (root / "pairlist.csv").exists()
    assert detect_layout_profile(root) == "docking_legacy"
    assert pairlist_path(root) == legacy_pairlist.resolve()
    assert len(pd.read_csv(legacy_pairlist)) == 1
    assert _dry_run_vina(root, "docking_legacy") == 0


def test_interactive_project_layout_follows_manifest_profile(tmp_path):
    canonical = tmp_path / "canonical"
    execution.run_workflow_init(str(canonical), layout_profile="canonical")
    layout = _project_layout(canonical)
    assert layout["layout_profile"] == "canonical"
    assert layout["prepared_proteins"] == (canonical / "prepared_proteins").resolve()
    assert not (canonical / "4-Docking").exists()

    legacy = tmp_path / "legacy"
    execution.run_workflow_init(str(legacy), layout_profile="docking_legacy")
    assert _project_layout(legacy)["layout_profile"] == "docking_legacy"


def _write_manifest(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def test_resolver_manifest_profile_wins_over_folder_presence(tmp_path):
    # A canonical project where a legacy folder exists: folder presence alone says legacy.
    root = tmp_path / "project"
    _write_manifest(root / "project_manifest.json", {"layout_profile": "canonical"})
    (root / "4-Docking").mkdir()

    assert detect_layout_profile(root) == "canonical"
    assert resolve_layout_profile(root) == "canonical"
    assert pairlist_path(root) == (root / "pairlist.csv").resolve()

    # A legacy manifest is honoured as recorded.
    legacy_root = tmp_path / "legacy_project"
    _write_manifest(legacy_root / "4-Docking" / "project_manifest.json", {"layout_profile": "docking_legacy"})
    assert detect_layout_profile(legacy_root) == "docking_legacy"
    assert resolve_layout_profile(legacy_root) == "docking_legacy"


def test_resolver_prefers_manifest_at_its_recorded_location(tmp_path):
    # A stray canonical manifest left in 4-Docking/ by an older run must not override
    # the canonical manifest at the project root.
    root = tmp_path / "project"
    _write_manifest(root / "project_manifest.json", {"layout_profile": "canonical"})
    _write_manifest(root / "4-Docking" / "project_manifest.json", {"layout_profile": "canonical"})
    assert resolve_layout_profile(root) == "canonical"
    assert detect_layout_profile(root) == "canonical"

    # Damaged project: the only manifest is the stray canonical one inside 4-Docking/.
    damaged = tmp_path / "damaged"
    _write_manifest(damaged / "4-Docking" / "project_manifest.json", {"layout_profile": "canonical"})
    assert resolve_layout_profile(damaged) == "canonical"
    assert detect_layout_profile(damaged) == "canonical"


def test_resolver_defaults_and_fallbacks(tmp_path):
    # No manifest: steps keep the legacy default; the low-level heuristic is unchanged.
    empty = tmp_path / "empty"
    empty.mkdir()
    assert resolve_layout_profile(empty) == "docking_legacy"
    assert detect_layout_profile(empty) == "canonical"

    # Corrupt manifest or a manifest without a recorded profile falls back to the heuristic.
    corrupt = tmp_path / "corrupt"
    (corrupt / "project_manifest.json").parent.mkdir(parents=True)
    (corrupt / "project_manifest.json").write_text("{not json", encoding="utf-8")
    assert resolve_layout_profile(corrupt) == "docking_legacy"
    assert detect_layout_profile(corrupt) == "canonical"

    unrecorded = tmp_path / "unrecorded"
    _write_manifest(unrecorded / "4-Docking" / "project_manifest.json", {"engines": ["vina"]})
    assert resolve_layout_profile(unrecorded) == "docking_legacy"
    assert detect_layout_profile(unrecorded) == "docking_legacy"

    # An explicit request is still honoured by the low-level helper.
    assert detect_layout_profile(empty, "docking_legacy") == "docking_legacy"
    with pytest.raises(ValueError):
        detect_layout_profile(empty, "unknown_profile")
