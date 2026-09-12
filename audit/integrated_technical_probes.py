"""Exercise actual pipeline registration and selected existing smoke contracts."""
import contextlib
import io
import importlib.util
import json
import os
import sys
import tempfile
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pandas as pd
spec = importlib.util.spec_from_file_location("omnidock_smoke_audit", ROOT / "test" / "test_dockforge_smoke.py")
smoke = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smoke)


def actual_raw_score_cache():
    with tempfile.TemporaryDirectory(prefix="omnidock_integrated_audit_") as tmp:
        project = Path(tmp)
        # Use the existing-directory compatibility path on Windows, where
        # ordinary users cannot necessarily create symbolic links.
        (project / "5-Post_Docking_Analysis").mkdir()
        smoke.bootstrap_project_layout(project, ["vina"], layout_profile="docking_legacy")
        smoke._write_test_pairlist(project)
        smoke._seed_engine_outputs(project, ("vina",))
        output_dir = project / "audit_session"

        def make_pipeline():
            return smoke.MultiEngineAnalysisPipeline(project_dir=str(project), output_dir=str(output_dir), analysis_mode="single_engine", engine="vina", engines_in_scope=["vina"])

        pipeline = make_pipeline()
        paths = pipeline._dag_artifact_paths()
        graph = pipeline.build_artifact_graph()
        first = graph.request(str(paths["raw_scores"]))
        before = pd.read_csv(paths["raw_scores"])["affinity_kcal_mol"].tolist()
        pose = smoke.ensure_engine_layout(project, "vina")["poses"] / "R1_site_1_L1.pdbqt"
        pose.write_text(pose.read_text().replace("-7.20", "-11.20"))
        stat = pose.stat()
        os.utime(pose, (stat.st_atime + 2, stat.st_mtime + 2))
        second_pipeline = make_pipeline()
        second_graph = second_pipeline.build_artifact_graph()
        second = second_graph.request(str(paths["raw_scores"]))
        cached = pd.read_csv(paths["raw_scores"])["affinity_kcal_mol"].tolist()
        forced = second_graph.request(str(paths["raw_scores"]), force=True)
        fresh = pd.read_csv(paths["raw_scores"])["affinity_kcal_mol"].tolist()
        return {"first_status": first["nodes"]["raw_scores"]["status"], "second_status": second["nodes"]["raw_scores"]["status"], "before": before, "after_raw_pose_changed": cached, "forced_reparse": fresh, "registered_inputs": graph.nodes["raw_scores"].inputs}


if __name__ == "__main__":
    results = {"actual_pipeline_cache": actual_raw_score_cache(), "selected_existing_smokes": {}}
    names = [
        "_smoke_artifact_graph_cache_hit",
        "_smoke_artifact_graph_registration_contract",
        "_smoke_consensus_direction_and_single_engine_qc",
        "_smoke_top_pose_consensus_direction_contracts",
        "_smoke_geometric_consensus_soft_alignment_contract",
        "_smoke_geometric_consensus_hard_mismatch_contract",
        "_smoke_redocking_multi_pose_sdf_parser_contract",
        "_smoke_redocking_validation_multi_pose_end_to_end_contract",
    ]
    captured = io.StringIO()
    for name in names:
        try:
            with contextlib.redirect_stdout(captured), contextlib.redirect_stderr(captured):
                getattr(smoke, name)()
        except Exception as exc:
            results["selected_existing_smokes"][name] = {"status": "failed", "error": str(exc)}
            traceback.print_exc(file=captured)
        else:
            results["selected_existing_smokes"][name] = {"status": "passed"}
    text = json.dumps(results, indent=2)
    print(text)
    (ROOT / "audit" / "integrated-technical-results.json").write_text(text + "\n", encoding="utf-8")
    (ROOT / "audit" / "selected-smoke-run.log").write_text(captured.getvalue(), encoding="utf-8")
