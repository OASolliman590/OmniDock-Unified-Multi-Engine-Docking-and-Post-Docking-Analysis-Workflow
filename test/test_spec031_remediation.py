from __future__ import annotations

import hashlib
import inspect
import json
from pathlib import Path

import numpy as np
import pandas as pd

from docking.preparation import pairlist_builder
from core_pipeline import MolecularDockingPipeline
from post_docking_analysis.docking_parser import parse_vina_pdbqt
from post_docking_analysis.report_generator import _safe_relative
from post_docking_analysis.value_normalization import normalize_boolean, normalize_boolean_series
from post_docking_analysis.visualization_suite import (
    PlotContext,
    _as_relative,
    _plot_hit_class_matrix,
)
from workflow import execution
from workflow.state import write_meta_run_manifest


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_MANIFEST = REPO_ROOT / "test" / "fixtures" / "spec031" / "vina_1iep" / "fixture_manifest.json"
EVIDENCE_MANIFEST = REPO_ROOT / "specs" / "031-appraisal-remediation-p0" / "evidence" / "manifest.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def test_vina_1iep_fixture_hashes_and_direct_parser_rows() -> None:
    payload = json.loads(FIXTURE_MANIFEST.read_text(encoding="utf-8"))
    paths = [payload["receptor"], payload["prepared_receptor"], payload["ligand"], payload["pose_output"]]
    for record in paths:
        artifact = REPO_ROOT / record["path"]
        assert artifact.is_file(), record["path"]
        assert _sha256(artifact) == record["sha256"]

    pose_output = REPO_ROOT / payload["pose_output"]["path"]
    parsed = parse_vina_pdbqt(pose_output)
    assert len(parsed) == payload["pose_output"]["retained_parsed_rows"] == 2
    assert parsed["pose"].tolist() == [1, 2]
    assert parsed["vina_affinity"].tolist() == [-11.702, -8.898]
    assert payload["ligand"]["protonation_ph"] == "unknown"
    assert payload["ligand"]["tautomer_state"] == "unknown"
    assert payload["ligand"]["charge_method"] == "unknown"


def test_evidence_manifest_classifies_every_finding() -> None:
    payload = json.loads(EVIDENCE_MANIFEST.read_text(encoding="utf-8"))
    expected = {
        *(f"F-{index:02d}" for index in range(1, 8)),
        *(f"NMR-{index:02d}" for index in range(1, 10)),
        "M-01",
        "M-02",
    }
    findings = {row["id"]: row for row in payload["findings"]}
    assert set(findings) == expected
    assert all(row["classification"] in {"verified_code", "observed_artifact", "reported_unverified"} for row in findings.values())
    assert payload["superseded_claims"][0]["classification"] == "superseded_claim"
    assert payload["nmrbox_evidence"]["classification"] == "reported_unverified"


def test_hit_class_matrix_executes_and_writes_nonempty_figure(tmp_path: Path) -> None:
    classified = pd.DataFrame(
        [
            {"ligand": "L1", "protein": "P1", "docking_quality_class": "Strong", "best_affinity_kcal_mol": -9.0},
            {"ligand": "L1", "protein": "P2", "docking_quality_class": "Moderate", "best_affinity_kcal_mol": -8.0},
            {"ligand": "L2", "protein": "P1", "docking_quality_class": "Weak", "best_affinity_kcal_mol": -6.0},
            {"ligand": "L2", "protein": "P2", "docking_quality_class": "Inactive", "best_affinity_kcal_mol": -4.0},
        ]
    )
    ctx = PlotContext(
        output_root=tmp_path,
        classified_hits=classified,
        consensus_ranked=pd.DataFrame(),
        engine_agreement=pd.DataFrame(),
        normalized_scores=pd.DataFrame(),
        top_pose_global=pd.DataFrame(),
        engine_scope_config={},
        validation_gate={},
        validation_table=pd.DataFrame(),
        reference_baselines=pd.DataFrame(),
        is_multi_engine=True,
        engine_name=None,
        primary_score_col="affinity_kcal_mol",
        primary_score_label="affinity_kcal_mol",
        has_references=False,
        validation_passed=False,
        validation_status="unknown",
        engine_count=2,
    )
    output = tmp_path / "hit_class_matrix.png"
    _plot_hit_class_matrix(ctx, output)
    assert output.is_file()
    assert output.stat().st_size > 0


def test_manifest_relative_paths_are_posix(tmp_path: Path) -> None:
    nested = tmp_path / "0-Input" / "receptors" / "source.pdbqt"
    nested.parent.mkdir(parents=True)
    nested.write_text("ATOM\n", encoding="utf-8")
    manifest_path = write_meta_run_manifest(
        tmp_path,
        input_paths=[nested],
        tool_versions={"python": "test"},
    )
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["input_checksums"][0]["path"] == "0-Input/receptors/source.pdbqt"
    assert _as_relative(nested, tmp_path) == "0-Input/receptors/source.pdbqt"
    assert _safe_relative(nested, tmp_path) == "0-Input/receptors/source.pdbqt"


def test_boolean_normalization_rejects_false_like_values() -> None:
    values = [False, "False", "0", "", None, np.nan, "no", "off"]
    assert [normalize_boolean(value) for value in values] == [False] * len(values)
    assert [normalize_boolean(value) for value in [True, "True", "1", "yes", 1]] == [True] * 5
    normalized = normalize_boolean_series(pd.Series(["False", "true", np.nan, 0, 1]))
    assert normalized.tolist() == [False, True, False, False, True]


def test_pocket_outputs_disclose_heuristic_contract() -> None:
    pipeline = object.__new__(MolecularDockingPipeline)
    pipeline._load_structure = lambda *_args, **_kwargs: []
    result = pipeline.analyze_pocket_properties("unused.pdb", np.array([1.0, 2.0, 3.0]))
    assert result["descriptor_contract"] == "uncalibrated_pocket_heuristics_v1"
    assert result["pocket_volume_method"] == "fixed_sphere_radius_5A_heuristic"
    assert result["druggability_score_method"] == "uncalibrated_descriptor_mean_v1"
    assert result["druggability_score_calibrated"] is False
    assert result["druggability_interpretation"].startswith("Heuristic ")
    assert "not calibrated" in result["scientific_limitation"]


def test_pairlist_preserves_cocrystal_metadata(monkeypatch, tmp_path: Path) -> None:
    proteins = tmp_path / "prepared_proteins"
    ligands = tmp_path / "prepared_ligands"
    proteins.mkdir()
    ligands.mkdir()
    (proteins / "1IEP_receptor.pdbqt").write_text("ATOM\n", encoding="utf-8")
    (ligands / "1IEP_ligand_STI_A_201.pdbqt").write_text("ATOM\n", encoding="utf-8")

    site_catalog = pd.DataFrame(
        [
            {
                "pdb_id": "1IEP",
                "selected_ligand": "1IEP_ligand_STI_A_201.pdbqt",
                "center_x": 1.0,
                "center_y": 2.0,
                "center_z": 3.0,
            }
        ]
    )
    monkeypatch.setattr(
        pairlist_builder,
        "load_site_catalog_from_summary",
        lambda *_args, **_kwargs: (site_catalog, []),
    )
    monkeypatch.setattr(
        pairlist_builder,
        "ensure_project_alias_files",
        lambda **_kwargs: {"protein_alias_file": "", "ligand_alias_file": ""},
    )
    monkeypatch.setattr(pairlist_builder, "build_protein_alias_lookup", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(pairlist_builder, "build_ligand_alias_lookup", lambda *_args, **_kwargs: {})

    summary = pairlist_builder.build_pairlists(
        project_root=tmp_path / "project",
        prepared_proteins=proteins,
        prepared_ligands=ligands,
        excel_path=tmp_path / "unused.xlsx",
        mode="cocrystal_only",
    )
    pairlist = pd.read_csv(summary["pairlist_file"])
    assert pairlist_builder.PAIRLIST_COLUMNS == list(pairlist.columns)
    row = pairlist.iloc[0]
    assert normalize_boolean(row["is_cocrystal_benchmark"])
    assert row["pair_source"] == "cocrystal"
    assert row["pdb_id"] == "1IEP"
    assert row["cocrystal_ligand_name"] == "1IEP_ligand_STI_A_201.pdbqt"


def test_non_gnina_smoke_call_supplies_analysis_config() -> None:
    parameter = inspect.signature(execution._prepare_non_gnina_context).parameters["analysis_config"]
    assert parameter.default is inspect.Parameter.empty
    source = Path(REPO_ROOT / "test" / "test_dockforge_smoke.py").read_text(encoding="utf-8")
    call_start = source.index("result = workflow_execution._prepare_non_gnina_context(")
    call_end = source.index("\n        )", call_start)
    assert "analysis_config=None" in source[call_start:call_end]
