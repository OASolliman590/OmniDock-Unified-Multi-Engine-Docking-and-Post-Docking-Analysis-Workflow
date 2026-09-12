"""Regression tests for scientific artifact and persistent-state contracts."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import multiprocessing
import os
import pytest

from post_docking_analysis.artifact_graph import ArtifactGraph, ArtifactNode
from workflow import state
from workflow.locking import file_lock


def test_graph_blocks_every_descendant_of_failed_node(tmp_path):
    graph = ArtifactGraph(cache_file=tmp_path / "cache.json")
    outputs = [str(tmp_path / f"{name}.txt") for name in "abc"]
    calls = []

    def fail():
        raise ValueError("invalid source")

    graph.register(ArtifactNode("a", [], [outputs[0]], fail))
    graph.register(ArtifactNode("b", [outputs[0]], [outputs[1]], lambda: calls.append("b")))
    graph.register(ArtifactNode("c", [outputs[1]], [outputs[2]], lambda: calls.append("c")))
    report = graph.request(outputs[2])
    assert [report["nodes"][name]["status"] for name in "abc"] == ["failed", "blocked_by_failure", "blocked_by_failure"]
    assert calls == []


def test_lock_timeout_also_bounds_competing_threads(tmp_path):
    path = tmp_path / "state.lock"
    def contend():
        with file_lock(path, timeout=0.05):
            return True
    with ThreadPoolExecutor(max_workers=1) as pool:
        with file_lock(path):
            future = pool.submit(contend)
            with pytest.raises(TimeoutError):
                future.result(timeout=2.0)


def test_graph_rejects_success_without_declared_output(tmp_path):
    graph = ArtifactGraph(cache_file=tmp_path / "cache.json")
    output = str(tmp_path / "missing.txt")
    graph.register(ArtifactNode("writer", [], [output], lambda: True))
    assert graph.request(output)["nodes"]["writer"]["status"] == "failed"


def test_graph_rejects_false_even_with_old_output(tmp_path):
    output = tmp_path / "old.txt"
    output.write_text("old data")
    graph = ArtifactGraph(cache_file=tmp_path / "cache.json")
    graph.register(ArtifactNode("writer", [], [str(output)], lambda: False))
    assert graph.request(str(output))["nodes"]["writer"]["status"] == "failed"


def test_graph_missing_external_input_blocks_callback(tmp_path):
    graph = ArtifactGraph(cache_file=tmp_path / "cache.json")
    output = str(tmp_path / "output.txt")
    calls = []
    graph.register(ArtifactNode("writer", [str(tmp_path / "missing-source.txt")], [output], lambda: calls.append(True)))
    report = graph.request(output)
    assert report["nodes"]["writer"]["status"] == "failed"
    assert not calls


def test_old_output_does_not_satisfy_successful_noop(tmp_path):
    output = tmp_path / "old.txt"
    output.write_text("previous scientific result")
    graph = ArtifactGraph(cache_file=tmp_path / "cache.json")
    graph.register(ArtifactNode("writer", [], [str(output)], lambda: True))
    assert graph.request(str(output))["nodes"]["writer"]["status"] == "failed"
    assert not output.exists()
    assert any(path.read_text() == "previous scientific result" for path in (tmp_path / ".artifact_history").rglob("old.txt"))


def test_failed_optional_recomputation_cannot_feed_stale_output(tmp_path):
    optional, result = tmp_path / "optional.txt", tmp_path / "report.txt"
    optional.write_text("old interactions")
    graph = ArtifactGraph(cache_file=tmp_path / "cache.json")
    def fail():
        optional.write_text("partial new interactions")
        raise ValueError("backend failed")
    graph.register(ArtifactNode("optional", [], [str(optional)], fail, optional=True))
    graph.register(ArtifactNode("report", [str(optional)], [str(result)],
                                lambda: result.write_text(optional.read_text() if optional.exists() else "unavailable"),
                                optional_inputs=[str(optional)]))
    report = graph.request(str(result))
    assert report["nodes"]["optional"]["status"] == "skipped_optional"
    assert report["nodes"]["report"]["status"] == "completed"
    assert result.read_text() == "unavailable"


def test_cache_detects_same_size_content_change_with_preserved_timestamp(tmp_path):
    source, output = tmp_path / "source.txt", tmp_path / "result.txt"
    source.write_text("AAAA")
    graph = ArtifactGraph(cache_file=tmp_path / "cache.json")
    graph.register(ArtifactNode("copy", [str(source)], [str(output)], lambda: output.write_text(source.read_text())))
    graph.request(str(output))
    previous = source.stat()
    source.write_text("BBBB")
    os.utime(source, ns=(previous.st_atime_ns, previous.st_mtime_ns))
    assert graph.request(str(output))["nodes"]["copy"]["status"] != "cache_hit"
    assert output.read_text() == "BBBB"


def test_cache_detects_nested_directory_content_change(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    nested = source / "nested.txt"
    nested.write_text("old")
    output = tmp_path / "result.txt"
    graph = ArtifactGraph(cache_file=tmp_path / "cache.json")
    graph.register(ArtifactNode("copy", [str(source)], [str(output)], lambda: output.write_text(nested.read_text())))
    graph.request(str(output))
    nested.write_text("new")
    assert graph.request(str(output))["nodes"]["copy"]["status"] != "cache_hit"
    assert output.read_text() == "new"


def test_cache_rejects_modified_output(tmp_path):
    source, output = tmp_path / "input.txt", tmp_path / "output.txt"
    source.write_text("expected")
    graph = ArtifactGraph(cache_file=tmp_path / "cache.json")
    graph.register(ArtifactNode("copy", [str(source)], [str(output)], lambda: output.write_text(source.read_text())))
    graph.request(str(output))
    output.write_text("corrupted")
    assert graph.request(str(output))["nodes"]["copy"]["status"] != "cache_hit"
    assert output.read_text() == "expected"


def test_input_mutation_during_computation_is_not_cached_as_success(tmp_path):
    source, output = tmp_path / "source.txt", tmp_path / "output.txt"
    source.write_text("initial")
    graph = ArtifactGraph(cache_file=tmp_path / "cache.json")

    def compute():
        output.write_text(source.read_text())
        source.write_text("changed while running")

    graph.register(ArtifactNode("unstable", [str(source)], [str(output)], compute))
    result = graph.request(str(output))
    assert result["nodes"]["unstable"]["status"] == "failed"
    assert "changed during" in result["nodes"]["unstable"]["details"]


def test_optional_missing_input_is_explicit_and_does_not_block_report(tmp_path):
    missing, output = str(tmp_path / "optional.txt"), str(tmp_path / "report.txt")
    graph = ArtifactGraph(cache_file=tmp_path / "cache.json")
    graph.register(ArtifactNode("report", [missing], [output], lambda: Path(output).write_text("optional feature unavailable"), optional_inputs=[missing]))
    assert graph.request(output)["nodes"]["report"]["status"] == "completed"


def test_state_concurrent_threads_preserve_all_updates(tmp_path):
    state.ensure_state(tmp_path)
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = [pool.submit(state.update_context, tmp_path, **{f"key_{i}": i}) for i in range(40)]
        for future in futures:
            future.result()
    context = state.load_state(tmp_path)["current_context"]
    assert all(context.get(f"key_{i}") == i for i in range(40))


def _state_process_writer(root, start):
    for i in range(start, start + 8):
        state.update_context(Path(root), **{f"process_{i}": i})


def test_state_concurrent_processes_preserve_all_updates(tmp_path):
    state.ensure_state(tmp_path)
    ctx = multiprocessing.get_context("spawn")
    processes = [ctx.Process(target=_state_process_writer, args=(str(tmp_path), offset)) for offset in (0, 8, 16)]
    for process in processes:
        process.start()
    for process in processes:
        process.join(30)
        assert process.exitcode == 0
    context = state.load_state(tmp_path)["current_context"]
    assert all(context.get(f"process_{i}") == i for i in range(24))


def test_checkpoint_order_is_newest_first_even_when_timestamps_tie(tmp_path, monkeypatch):
    monkeypatch.setattr(state, "_now_iso", lambda: "2026-09-12T10:00:00Z")
    for identifier in ("first", "second"):
        state.record_checkpoint_metadata(tmp_path, source_project_dir="a", target_project_dir="b", layout_profile="canonical", checkpoint_id=identifier)
    assert [row["checkpoint_id"] for row in state.list_checkpoint_metadata(tmp_path)] == ["second", "first"]
