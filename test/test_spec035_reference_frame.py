"""Spec 035 S18 + D17: the reference frame guard, exercised through the preparation builder and the DAG.

The receptor frame comes from its prepared-asset sidecar. The reference frame is supplied in the pair
intent. The DAG then writes the merged pair metadata into reports/combined_engine_scores.csv, and that
file is what the guard classifies.
"""
from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import pandas as pd

from docking.preparation.project_builder import DockingPreparationConfig, DockingProjectBuilder
from docking.project_layout import ensure_engine_layout
from post_docking_analysis.reference_policy import classify_reference_row
from post_docking_analysis.replicate_names import pair_tag
from post_docking_analysis.unified_pipeline import UnifiedPostDockingPipeline
from test_dockforge_smoke import _stub_dag_optional_interactions
from test_preparation_correctness import ligand_file, pdbqt, sdf

RECEPTOR = "1ABC_cleaned.pdbqt"
LIGAND = "1ABC_ligand_LIG.pdbqt"
SITE = "native"


def _build_and_analyse(root: Path, *, receptor_frame: str, reference_frame: str) -> pd.DataFrame:
    receptors, ligands = root / "receptors", root / "ligands"
    receptors.mkdir(parents=True)
    ligands.mkdir(parents=True)
    receptor = receptors / RECEPTOR
    receptor.write_text(pdbqt(1, record="ATOM", name="CA", resname="ALA") + "\n")
    receptor.with_suffix(".pdbqt.preparation.json").write_text(json.dumps({"receptor_frame_id": receptor_frame}))
    ligand_file(ligands / LIGAND)
    native = sdf(root / "1ABC_native.sdf")
    # The reference is declared in the pairlist, so the builder keeps it; the frame guard then decides.
    pairlist = root / "input.csv"
    pd.DataFrame([dict(
        receptor=RECEPTOR, site_id=SITE, ligand=LIGAND,
        center_x=0.0, center_y=0.0, center_z=0.0, size_x=20.0, size_y=20.0, size_z=20.0,
        is_cocrystal_benchmark=True, reference_source="experimental:1ABC", reference_pdb_id="1ABC",
        reference_frame_id=reference_frame, reference_pose_file=str(native),
    )]).to_csv(pairlist, index=False)
    project = root / "project"
    built = DockingProjectBuilder(DockingPreparationConfig(
        prepared_proteins=receptors, prepared_ligands=ligands, pair_intent=None, excel_path=None,
        output_dir=project, pairlist_file=pairlist, asset_mode="copy", engines=["vina"],
    )).build()
    pair_row = pd.read_csv(built["pairlist_file"]).iloc[0]
    assert pair_row["receptor_frame_id"] == receptor_frame
    assert pair_row["reference_frame_id"] == reference_frame

    # One Vina score for the pair, seeded in the engine's normalized-score layout (as the smoke fixtures do).
    tag = pair_tag(RECEPTOR, SITE, LIGAND)
    layout = ensure_engine_layout(project, "vina")
    (layout["poses"] / f"{tag}.pdbqt").write_text("REMARK\n", encoding="utf-8")
    (layout["logs"] / f"{tag}.log").write_text("REMARK\n", encoding="utf-8")
    pd.DataFrame([{
        "engine": "vina", "tag": tag, "protein": RECEPTOR, "ligand": LIGAND, "site_id": SITE,
        "pose": 1, "affinity_kcal_mol": -8.0, "score_name_primary": "vina_affinity", "score_primary": -8.0,
        "score_name_secondary": "", "score_secondary": None, "rmsd_lb": None, "rmsd_ub": None,
        "pose_file": str(layout["poses"] / f"{tag}.pdbqt"), "log_file": str(layout["logs"] / f"{tag}.log"),
    }]).to_csv(layout["scores"] / "normalized_scores.csv", index=False)

    pipeline = UnifiedPostDockingPipeline(
        project_dir=str(project),
        output_dir=str(project / "analysis" / "sessions" / "s18"),
        analysis_mode="comparative_all_engines",
        analysis_scope="full",
        engines="vina",
        rmsd_workers=2,
    )
    _stub_dag_optional_interactions(pipeline)
    assert pipeline.run() is True
    report = pipeline.request_artifact_scope("full")
    statuses = {name: payload["status"] for name, payload in report["nodes"].items()}
    assert statuses["comparative"] in {"completed", "completed_with_warnings"}, statuses
    assert statuses["top_pose_atlas"] in {"completed", "completed_with_warnings", "cache_hit"}, statuses
    combined = pd.read_csv(pipeline.output_dir / "reports" / "combined_engine_scores.csv")
    rows = combined[combined["tag"] == tag]
    assert len(rows) == 1
    return rows


def test_cross_frame_reference_is_cleared_through_the_pipeline(tmp_path):
    root = Path(tempfile.mkdtemp(prefix="s18_cross_"))
    try:
        rows = _build_and_analyse(root, receptor_frame="frame2", reference_frame="frame1")
        row = rows.iloc[0].to_dict()
        assert row["reference_frame_id"] == "frame1" and row["receptor_frame_id"] == "frame2"
        classification = classify_reference_row(row)
        assert classification.is_reference is False
        assert classification.validation_anchor_eligible is False
        assert classification.reason == "reference_frame_differs_from_receptor"
        assert classification.identity_fields == {} and classification.lineage_fields == {}
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_same_frame_reference_is_kept_through_the_pipeline(tmp_path):
    # Control: the same frame must not be cleared, otherwise the guard would hide every reference.
    root = Path(tempfile.mkdtemp(prefix="s18_same_"))
    try:
        rows = _build_and_analyse(root, receptor_frame="frame1", reference_frame="frame1")
        row = rows.iloc[0].to_dict()
        classification = classify_reference_row(row)
        assert classification.is_reference is True
        assert classification.reason != "reference_frame_differs_from_receptor"
    finally:
        shutil.rmtree(root, ignore_errors=True)
