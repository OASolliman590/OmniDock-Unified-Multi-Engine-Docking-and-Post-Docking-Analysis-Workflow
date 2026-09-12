"""Non-destructive audit counterexamples; run from the repository with Python.

These probes describe current defects rather than asserting desired behavior.
All generated workflow data is confined to temporary directories.
"""
from __future__ import annotations

import ast
import json
import sys
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from post_docking_analysis.artifact_graph import ArtifactGraph, ArtifactNode
from workflow import state


def transitive_failure(root):
    root.mkdir(parents=True)
    graph = ArtifactGraph(cache_file=root / "cache.json")
    a, b, c = [str(root / f"{name}.txt") for name in "abc"]

    def fail():
        raise RuntimeError("upstream parsing failed")

    graph.register(ArtifactNode("a", [], [a], fail))
    graph.register(ArtifactNode("b", [a], [b], lambda: Path(b).write_text("b")))
    graph.register(ArtifactNode("c", [b], [c], lambda: Path(c).write_text("report complete")))
    result = graph.request(c)
    return {k: v["status"] for k, v in result["nodes"].items()}


def missing_output(root):
    graph = ArtifactGraph(cache_file=root / "cache.json")
    out = str(root / "never-created.txt")
    graph.register(ArtifactNode("no_output", [], [out], lambda: False))
    report = graph.request(out)
    return {"status": report["nodes"]["no_output"]["status"], "output_exists": Path(out).exists()}


def undeclared_raw_input(root):
    """Use the exact raw_scores dependency shape registered by the pipeline.

    The compute callback is deliberately a tiny stand-in for the engine parser;
    the DAG executor and registration source are the actual repository code.
    """
    root.mkdir(parents=True)
    pairlist, scope, raw, output = [root / n for n in ("pairlist.csv", "scope.json", "pose.pdbqt", "scores.csv")]
    pairlist.write_text("same pair")
    scope.write_text("same scope")
    raw.write_text("-7")
    graph = ArtifactGraph(cache_file=root / "cache.json")
    graph.register(ArtifactNode("raw_scores", [str(pairlist), str(scope)], [str(output)], lambda: output.write_text(raw.read_text())))
    graph.request(str(output))
    raw.write_text("-12")
    report = graph.request(str(output))
    return {"status": report["nodes"]["raw_scores"]["status"], "actual_raw_score": raw.read_text(), "reported_score": output.read_text()}


def lost_state_update(root):
    state.ensure_state(root)
    original = state.ensure_state
    barrier = threading.Barrier(2)

    def synchronized_read(path):
        snapshot = original(path)
        barrier.wait(timeout=10)
        return snapshot

    state.ensure_state = synchronized_read
    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(state.update_context, root, **{key: True}) for key in ("audit_worker_a", "audit_worker_b")]
            for future in futures:
                future.result()
    finally:
        state.ensure_state = original
    context = state.load_state(root)["current_context"]
    return {"expected_keys": ["audit_worker_a", "audit_worker_b"], "persisted_keys": sorted(k for k in context if k.startswith("audit_worker_"))}


def inventory():
    files = [p for p in ROOT.rglob("*.py") if "audit" not in p.relative_to(ROOT).parts and ".git" not in p.parts]
    counts = sorted(((len(p.read_text(encoding="utf-8-sig").splitlines()), p.relative_to(ROOT).as_posix()) for p in files), reverse=True)
    tests = {}
    for name in ("test/test_dockforge_smoke.py", "test_pipeline.py"):
        module = ast.parse((ROOT / name).read_text(encoding="utf-8-sig"))
        functions = [n.name for n in module.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
        tests[name] = {"pytest_test_functions": sum(n.startswith("test_") for n in functions), "custom_smoke_functions": sum(n.startswith("_smoke_") for n in functions)}
    return {"python_files": len(files), "python_lines": sum(n for n, _ in counts), "largest_modules": counts[:8], "test_discovery": tests}


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="omnidock_audit_") as tmp:
        root = Path(tmp)
        findings = {"transitive_failure": transitive_failure(root / "failure"), "missing_output": missing_output(root / "missing"), "undeclared_raw_input": undeclared_raw_input(root / "cache"), "concurrent_state_update": lost_state_update(root / "state"), "inventory": inventory()}
    text = json.dumps(findings, indent=2)
    print(text)
    (ROOT / "audit" / "technical-probe-results.json").write_text(text + "\n", encoding="utf-8")
