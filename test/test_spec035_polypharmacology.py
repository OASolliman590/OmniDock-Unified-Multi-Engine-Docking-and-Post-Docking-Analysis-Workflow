"""Spec 035 D10 (S26): cross-engine polypharmacology is not evaluable, because raw affinities are not comparable."""
import json
import logging

import pandas as pd

from post_docking_analysis.simplified_pipeline_impl import SimplifiedPostDockingPipeline


def _pipeline(tmp_path, best_poses):
    pipeline = object.__new__(SimplifiedPostDockingPipeline)
    pipeline.output_dir = tmp_path / "out"
    pipeline.logger = logging.getLogger("test_spec035_polypharmacology")
    pipeline.results = {}
    pipeline.polypharm_config = {"enabled": True}
    pipeline._get_best_poses_for_visualization = lambda: best_poses
    return pipeline


def _rows(engines):
    return pd.DataFrame(
        [
            {"engine": engine, "ligand": "L1", "protein": "P1", "protein_label": "P1",
             "vina_affinity": -8.0 - index, "tag": f"P1_site_1_L1_{engine}"}
            for index, engine in enumerate(engines)
        ]
    )


def test_cross_engine_polypharmacology_is_not_evaluable(tmp_path):
    pipeline = _pipeline(tmp_path, _rows(["vina", "smina"]))
    assert pipeline._analyze_polypharmacology() is True
    status = json.loads((tmp_path / "out" / "analysis" / "polypharmacology" / "status.json").read_text(encoding="utf-8"))
    assert status["status"] == "not_evaluable"
    assert status["reason"] == "cross_engine_raw_score_comparison_requires_calibration"
    assert pipeline.results["polypharmacology_analysis"] == status
    assert not (tmp_path / "out" / "analysis" / "polypharmacology" / "polypharmacology_pairs.csv").exists()
