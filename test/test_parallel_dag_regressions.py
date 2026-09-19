from __future__ import annotations

import threading
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

from post_docking_analysis.multi_engine_pipeline_impl import MultiEngineAnalysisPipeline
from post_docking_analysis.visualization_suite import _save_placeholder


def _pipeline(output_dir: Path) -> MultiEngineAnalysisPipeline:
    pipeline = MultiEngineAnalysisPipeline.__new__(MultiEngineAnalysisPipeline)
    pipeline.output_dir = output_dir
    pipeline.analysis_scope = "full"
    return pipeline


def _dag_paths(output_dir: Path) -> dict[str, Path]:
    reports = output_dir / "reports"
    return {
        "comparative": reports / "comparative",
        "polypharmacology": reports / "polypharmacology",
        "biology_correlation": reports / "biology_correlation",
    }


def _forbid_producer_reentry(*_args, **_kwargs):
    raise AssertionError("a DAG consumer re-entered the comparative producer")


def test_consumers_materialize_from_comparative_owned_bundles_after_output_retirement(tmp_path):
    output_dir = tmp_path / "session"
    paths = _dag_paths(output_dir)
    poly_source = paths["comparative"] / "polypharmacology_source"
    biology_source = paths["comparative"] / "biology_correlation_source"
    poly_source.mkdir(parents=True)
    biology_source.mkdir(parents=True)
    (poly_source / "polypharmacology_ligand_summary_with_biology.csv").write_text(
        "ligand,consensus_score\nL1,0.75\n",
        encoding="utf-8",
    )
    (biology_source / "biology_mapping_report.json").write_text(
        '{"mapping_mode":"none","matched_rows":0}',
        encoding="utf-8",
    )

    for published in (paths["polypharmacology"], paths["biology_correlation"]):
        published.mkdir(parents=True)
        (published / "stale.txt").write_text("stale", encoding="utf-8")
        retired = tmp_path / "retired" / published.name
        retired.parent.mkdir(parents=True, exist_ok=True)
        published.rename(retired)

    pipeline = _pipeline(output_dir)
    pipeline._load_or_build_scores = _forbid_producer_reentry
    pipeline._write_comparative_reports = _forbid_producer_reentry

    original_poly = (poly_source / "polypharmacology_ligand_summary_with_biology.csv").read_bytes()
    original_biology = (biology_source / "biology_mapping_report.json").read_bytes()
    pipeline._dag_compute_polypharmacology_node(paths)
    pipeline._dag_compute_biology_correlation_node(paths)

    assert (paths["polypharmacology"] / "polypharmacology_ligand_summary_with_biology.csv").read_bytes() == original_poly
    assert (paths["biology_correlation"] / "biology_mapping_report.json").read_bytes() == original_biology
    assert (poly_source / "polypharmacology_ligand_summary_with_biology.csv").read_bytes() == original_poly
    assert (biology_source / "biology_mapping_report.json").read_bytes() == original_biology
    assert not (paths["polypharmacology"] / "stale.txt").exists()
    assert not (paths["biology_correlation"] / "stale.txt").exists()


@pytest.mark.parametrize(
    ("node_method", "output_key", "missing_name"),
    (
        ("_dag_compute_polypharmacology_node", "polypharmacology", "polypharmacology_source"),
        ("_dag_compute_biology_correlation_node", "biology_correlation", "biology_correlation_source"),
    ),
)
def test_consumers_fail_closed_when_comparative_owned_bundle_is_missing(
    tmp_path,
    node_method: str,
    output_key: str,
    missing_name: str,
):
    output_dir = tmp_path / "session"
    paths = _dag_paths(output_dir)
    paths["comparative"].mkdir(parents=True)
    published = paths[output_key]
    published.mkdir(parents=True)
    (published / "stale.txt").write_text("stale", encoding="utf-8")
    retired = tmp_path / "retired" / published.name
    retired.parent.mkdir(parents=True, exist_ok=True)
    published.rename(retired)

    pipeline = _pipeline(output_dir)
    pipeline._load_or_build_scores = _forbid_producer_reentry
    pipeline._write_comparative_reports = _forbid_producer_reentry

    with pytest.raises(FileNotFoundError, match=missing_name):
        getattr(pipeline, node_method)(paths)
    assert not published.exists()


def test_forced_comparative_rerun_refreshes_owned_source_bundles(tmp_path):
    output_dir = tmp_path / "session"
    paths = _dag_paths(output_dir)
    pipeline = _pipeline(output_dir)
    pipeline._load_or_build_scores = lambda: pd.DataFrame()
    generation = 0

    def write_generation(_scores, analysis_scope="full", *, preserve_dag_atlas=False):
        nonlocal generation
        generation += 1
        reports = output_dir / "reports"
        poly = reports / "polypharmacology"
        poly.mkdir(parents=True, exist_ok=True)
        (poly / "generation.txt").write_text(str(generation), encoding="utf-8")
        (reports / "biology_mapping_report.json").write_text(
            '{"generation":%d}' % generation,
            encoding="utf-8",
        )
        (reports / "biology_correlation_global.csv").write_text(
            "generation\n%d\n" % generation,
            encoding="utf-8",
        )
        (reports / "biology_correlation_per_protein.csv").write_text(
            "generation\n%d\n" % generation,
            encoding="utf-8",
        )

    pipeline._write_comparative_reports = write_generation

    pipeline._dag_compute_comparative_node(paths)
    assert (paths["comparative"] / "polypharmacology_source" / "generation.txt").read_text(encoding="utf-8") == "1"
    assert '"generation":1' in (
        paths["comparative"] / "biology_correlation_source" / "biology_mapping_report.json"
    ).read_text(encoding="utf-8")

    retired = tmp_path / "retired-comparative"
    paths["comparative"].rename(retired)
    pipeline._dag_compute_comparative_node(paths)

    assert (paths["comparative"] / "polypharmacology_source" / "generation.txt").read_text(encoding="utf-8") == "2"
    assert '"generation":2' in (
        paths["comparative"] / "biology_correlation_source" / "biology_mapping_report.json"
    ).read_text(encoding="utf-8")
    assert (retired / "polypharmacology_source" / "generation.txt").read_text(encoding="utf-8") == "1"


def test_matplotlib_critical_section_is_reentrant_and_exception_safe():
    from post_docking_analysis.plotting_lock import matplotlib_critical_section

    with matplotlib_critical_section():
        with matplotlib_critical_section():
            pass

    with pytest.raises(RuntimeError, match="controlled failure"):
        with matplotlib_critical_section():
            raise RuntimeError("controlled failure")

    acquired = threading.Event()

    def acquire_after_exception():
        with matplotlib_critical_section():
            acquired.set()

    thread = threading.Thread(target=acquire_after_exception, daemon=True)
    thread.start()
    thread.join(timeout=2)
    assert not thread.is_alive()
    assert acquired.is_set()


def test_actual_plotting_helpers_serialize_global_pyplot_state(monkeypatch, tmp_path):
    pipeline = MultiEngineAnalysisPipeline.__new__(MultiEngineAnalysisPipeline)
    paired_rows = pd.DataFrame(
        [
            {"engine": "gnina", "affinity_kcal_mol": -8.1},
            {"engine": "vina", "affinity_kcal_mol": -7.3},
            {"engine": "smina", "affinity_kcal_mol": -7.8},
        ]
    )
    engine_order = ["gnina", "vina", "smina"]
    engine_palette = {"gnina": "#4c72b0", "vina": "#55a868", "smina": "#c17d11"}

    original_subplots = plt.subplots
    original_savefig = plt.savefig
    a_at_save = threading.Event()
    release_a = threading.Event()
    b_started = threading.Event()
    b_created = threading.Event()
    records: dict[str, int] = {}
    errors: list[BaseException] = []

    def controlled_subplots(*args, **kwargs):
        fig, axes = original_subplots(*args, **kwargs)
        name = threading.current_thread().name
        records[f"{name}_created"] = id(fig)
        if name == "B":
            b_created.set()
        return fig, axes

    def controlled_savefig(*_args, **_kwargs):
        name = threading.current_thread().name
        if name == "A":
            records["A_save_before"] = id(plt.gcf())
            a_at_save.set()
            if not release_a.wait(timeout=5):
                raise TimeoutError("thread A was not released")
            records["A_save_after"] = id(plt.gcf())

    monkeypatch.setattr(plt, "subplots", controlled_subplots)
    monkeypatch.setattr(plt, "savefig", controlled_savefig)

    def run_comparative_plot():
        try:
            pipeline._plot_cross_engine_distribution(
                paired_rows,
                engine_order,
                engine_palette,
                tmp_path / "distribution.png",
            )
        except BaseException as exc:  # captured for the parent assertion
            errors.append(exc)

    def run_suite_plot():
        b_started.set()
        try:
            _save_placeholder(tmp_path / "placeholder.png", "Placeholder", "B figure")
        except BaseException as exc:  # captured for the parent assertion
            errors.append(exc)

    thread_a = threading.Thread(target=run_comparative_plot, name="A", daemon=True)
    thread_b = threading.Thread(target=run_suite_plot, name="B", daemon=True)
    thread_a.start()
    assert a_at_save.wait(timeout=5)
    thread_b.start()
    assert b_started.wait(timeout=2)
    try:
        assert not b_created.wait(timeout=0.5), "visualization-suite plotting entered while comparative pyplot state was active"
    finally:
        release_a.set()
        thread_a.join(timeout=5)
        thread_b.join(timeout=5)
        plt.close("all")

    assert not thread_a.is_alive()
    assert not thread_b.is_alive()
    assert not errors
    assert records["A_created"] == records["A_save_before"] == records["A_save_after"]
    assert "B_created" in records
