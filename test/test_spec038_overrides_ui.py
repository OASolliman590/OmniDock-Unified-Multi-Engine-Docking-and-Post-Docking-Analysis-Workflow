"""Spec 038 R4: protonation overrides in the interactive shell and the web UI.

Both surfaces call the Spec 037 resolver in ``docking.project_layout``:

- the view (``workflow.interactive.build_protonation_view``) only reads, through ``resolve_entity_protonation``
  and ``check_preparation_policy_current``;
- the interactive set flow stores through ``set_protonation_override``, the setter ``workflow protonation-policy``
  uses (checked by a spy, not by a re-implementation);
- the web UI shows the same view, and its set action is a launched CLI job (``workflow protonation-policy``),
  because the web process itself must not write into a project (FR-035).

Clearing an override is not covered: the resolver has no clear operation yet (see the report that added this file).
"""

from __future__ import annotations

import importlib
import json
import time
from pathlib import Path

import pytest

from docking.project_layout import (
    POLICY_HASH_VERSION,
    ProtonationPolicyConflict,
    resolve_entity_protonation,
    set_protonation_override,
    set_protonation_policy,
)
from test_spec036_protonation_policy import _write_manifest
from workflow import interactive
from workflow.cli import main as cli_main

REPO_ROOT = Path(__file__).resolve().parents[1]
STALE = "stale_preparation_policy_changed"


# ----------------------------------------------------------------------------- helpers


class _Prompt:
    def __init__(self, value):
        self._value = value

    def ask(self):
        return self._value


class _FakeQuestionary:
    """Queued answers, in the order the flow asks for them. Unqueued prompts return their default."""

    class Choice:
        def __init__(self, title, value, checked=False):
            self.title = title
            self.value = value
            self.checked = checked

    def __init__(self, *, selects=(), texts=(), confirms=()):
        self.selects = list(selects)
        self.texts = list(texts)
        self.confirms = list(confirms)

    def select(self, _message, choices=None, default=None, **_kwargs):
        return _Prompt(self.selects.pop(0) if self.selects else default)

    def text(self, _message, default=""):
        return _Prompt(self.texts.pop(0) if self.texts else default)

    def confirm(self, _message, default=True):
        return _Prompt(self.confirms.pop(0) if self.confirms else default)


def _project(tmp_path: Path) -> Path:
    """Canonical project with a stored default, two prepared receptors, two prepared ligands and one
    pair-status block in the manifest.

    - 1IEP: prepared under the current policy (current).
    - 2ABC: receptor override pH 6.0 set after it was prepared under the default (stale).
    - NEUT: ligand under the project default, prepared under it (current).
    - STI: ligand override ``as_input`` set after it was prepared under the default (stale).
    """
    project = tmp_path / "project"
    _write_manifest(project, layout_profile="canonical")
    set_protonation_policy(
        project, receptor_ph=7.4, receptor_force_field="AMBER", ligand_policy="ph_model", ligand_ph=7.4
    )
    before = {
        ("receptor", "2ABC"): resolve_entity_protonation(project, "receptor", "2ABC")["policy_hash"],
        ("ligand", "STI"): resolve_entity_protonation(project, "ligand", "STI")["policy_hash"],
    }
    set_protonation_override(project, kind="receptor", name="2ABC", receptor_ph=6.0)
    set_protonation_override(project, kind="ligand", name="STI", ligand_policy="as_input")
    current = {
        ("receptor", "1IEP"): resolve_entity_protonation(project, "receptor", "1IEP")["policy_hash"],
        ("ligand", "NEUT"): resolve_entity_protonation(project, "ligand", "NEUT")["policy_hash"],
    }
    _stamp_receptor(project, "1IEP.pdbqt", current[("receptor", "1IEP")])
    _stamp_receptor(project, "2ABC.pdbqt", before[("receptor", "2ABC")])
    _stamp_ligand(project, "NEUT", current[("ligand", "NEUT")])
    _stamp_ligand(project, "STI", before[("ligand", "STI")])
    manifest_file = project / "project_manifest.json"
    payload = json.loads(manifest_file.read_text(encoding="utf-8"))
    payload["pair_protonation"] = {
        "pairlist_file": str(project / "pairlist.csv"),
        "pair_count": 4,
        "status_counts": {"coherent": 1, "pair_ph_mismatch": 1, "explicit_state_no_ph": 2},
    }
    manifest_file.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return project


def _stamp_receptor(project: Path, name: str, policy_hash: str) -> None:
    prepared = project / "prepared_proteins"
    prepared.mkdir(parents=True, exist_ok=True)
    (prepared / name).write_text("ATOM\n", encoding="utf-8")
    record = {"policy_hash": policy_hash, "policy_hash_version": POLICY_HASH_VERSION, "is_valid": True}
    (prepared / f"{name}.preparation.json").write_text(json.dumps(record), encoding="utf-8")


def _stamp_ligand(project: Path, stem: str, policy_hash) -> None:
    prepared = project / "prepared_ligands"
    steps = prepared / "preparation_steps"
    steps.mkdir(parents=True, exist_ok=True)
    output = prepared / f"{stem}.pdbqt"
    output.write_text("ATOM\n", encoding="utf-8")
    record = {"output_file": str(output), "policy_hash": policy_hash, "policy_hash_version": POLICY_HASH_VERSION}
    (steps / f"{stem}.json").write_text(json.dumps(record), encoding="utf-8")


def _snapshot(root: Path) -> dict:
    return {
        p.relative_to(root).as_posix(): (p.is_dir(), p.stat().st_mtime_ns)
        for p in sorted(root.rglob("*"))
    }


def _row(rows, name):
    return next(row for row in rows if row["name"] == name)


# ----------------------------------------------------------------------------- view (read-only, resolver-backed)


def test_view_levels_come_from_the_resolver(tmp_path):
    project = _project(tmp_path)
    view = interactive.build_protonation_view(project)

    receptors = {row["name"]: row for row in view["receptors"]}
    ligands = {row["name"]: row for row in view["ligands"]}
    assert set(receptors) == {"1IEP", "2ABC"}
    assert set(ligands) == {"NEUT", "STI"}

    # Each value equals what the resolver returns for that entity, and the level is the resolver's level.
    for kind, rows in (("receptor", receptors), ("ligand", ligands)):
        for name, row in rows.items():
            resolved = resolve_entity_protonation(project, kind, name)
            assert row["policy_level"] == resolved["policy_level"]
            assert row["field_levels"] == resolved["field_levels"]
            assert row["policy_hash"] == resolved["policy_hash"]

    assert receptors["1IEP"]["receptor_ph"] == 7.4
    assert receptors["1IEP"]["field_levels"]["receptor_ph"] == "project"
    assert receptors["1IEP"]["policy_level"] == "project"
    assert receptors["2ABC"]["receptor_ph"] == 6.0
    assert receptors["2ABC"]["field_levels"]["receptor_ph"] == "receptor"
    assert receptors["2ABC"]["field_levels"]["receptor_force_field"] == "project"
    assert receptors["2ABC"]["policy_level"] == "receptor"

    assert ligands["NEUT"]["ligand_policy"] == "ph_model"
    assert ligands["NEUT"]["field_levels"]["ligand_policy"] == "project"
    assert ligands["STI"]["ligand_policy"] == "as_input"
    assert ligands["STI"]["field_levels"]["ligand_policy"] == "ligand"
    assert ligands["STI"]["policy_level"] == "ligand"


def test_view_reports_pair_status_counts_from_the_manifest(tmp_path):
    view = interactive.build_protonation_view(_project(tmp_path))
    assert view["pair"]["status_counts"] == {"coherent": 1, "pair_ph_mismatch": 1, "explicit_state_no_ph": 2}
    assert view["pair"]["pair_count"] == 4


def test_view_marks_stale_entities_and_their_reason(tmp_path):
    project = _project(tmp_path)
    before = _snapshot(project)
    view = interactive.build_protonation_view(project)

    assert _row(view["receptors"], "1IEP")["prepared"]["state"] == "current"
    assert _row(view["ligands"], "NEUT")["prepared"]["state"] == "current"
    assert _row(view["receptors"], "2ABC")["prepared"]["reason"] == STALE
    assert _row(view["ligands"], "STI")["prepared"]["reason"] == STALE
    assert {(item["kind"], item["name"], item["reason"]) for item in view["stale"]} == {
        ("receptor", "2ABC", STALE),
        ("ligand", "STI", STALE),
    }
    assert _snapshot(project) == before, "the view must not write into the project"


def test_view_refuses_a_preparation_record_without_a_policy_hash(tmp_path):
    project = _project(tmp_path)
    record = project / "prepared_proteins" / "1IEP.pdbqt.preparation.json"
    record.write_text(json.dumps({"is_valid": True}), encoding="utf-8")
    view = interactive.build_protonation_view(project)
    state = _row(view["receptors"], "1IEP")["prepared"]
    assert state["state"] == "refused"
    assert state["reason"] == "preparation_policy_unrecorded"


def test_view_without_a_manifest_says_so(tmp_path):
    view = interactive.build_protonation_view(tmp_path / "empty")
    assert view["has_manifest"] is False
    assert view["error"] == "no_manifest"


# ----------------------------------------------------------------------------- interactive shell


def test_interactive_receptor_override_is_stored_through_the_resolver(tmp_path, monkeypatch, capsys):
    project = _project(tmp_path)
    calls = []
    real_setter = set_protonation_override

    def spy(project_root, **kwargs):
        calls.append(kwargs)
        return real_setter(project_root, **kwargs)

    monkeypatch.setattr(interactive, "set_protonation_override", spy)
    # Name 1IEP, pH 6.5, force field kept (empty choice), store.
    fake = _FakeQuestionary(texts=["1IEP", "6.5"], selects=[""], confirms=[True])
    interactive._set_protonation_override_flow(fake, project, "receptor")

    assert calls == [{"kind": "receptor", "name": "1IEP", "replace": False, "receptor_ph": 6.5}]
    resolved = resolve_entity_protonation(project, "receptor", "1IEP")
    assert resolved["receptor_ph"] == 6.5
    assert resolved["field_levels"]["receptor_ph"] == "receptor"

    out = capsys.readouterr().out
    assert STALE in out, "the stale warning must be shown before the change is stored"
    assert "Now stale" in out and "1IEP" in out, "the entity that became stale must be listed after the change"


def test_interactive_ligand_override_is_stored_through_the_resolver(tmp_path, monkeypatch):
    project = _project(tmp_path)
    calls = []
    real_setter = set_protonation_override

    def spy(project_root, **kwargs):
        calls.append(kwargs)
        return real_setter(project_root, **kwargs)

    monkeypatch.setattr(interactive, "set_protonation_override", spy)
    # Name NEUT, policy ph_model, pH 6.0, store.
    fake = _FakeQuestionary(texts=["NEUT", "6.0"], selects=["ph_model"], confirms=[True])
    interactive._set_protonation_override_flow(fake, project, "ligand")

    assert calls == [{"kind": "ligand", "name": "NEUT", "replace": False, "ligand_policy": "ph_model", "ligand_ph": 6.0}]
    resolved = resolve_entity_protonation(project, "ligand", "NEUT")
    assert resolved["ligand_ph"] == 6.0
    assert resolved["field_levels"]["ligand_ph"] == "ligand"


def test_interactive_conflict_is_reported_and_refused_unless_replaced(tmp_path, capsys):
    project = _project(tmp_path)
    # 2ABC already has pH 6.0 stored. Asking for 7.0 conflicts; declining the replace keeps 6.0.
    fake = _FakeQuestionary(texts=["2ABC", "7.0"], selects=[""], confirms=[True, False])
    interactive._set_protonation_override_flow(fake, project, "receptor")

    out = capsys.readouterr().out
    assert "protonation_policy_conflict" in out
    assert "The override was not changed." in out
    assert resolve_entity_protonation(project, "receptor", "2ABC")["receptor_ph"] == 6.0


def test_interactive_conflict_replace_stores_the_new_value(tmp_path):
    project = _project(tmp_path)
    fake = _FakeQuestionary(texts=["2ABC", "7.0"], selects=[""], confirms=[True, True])
    interactive._set_protonation_override_flow(fake, project, "receptor")
    assert resolve_entity_protonation(project, "receptor", "2ABC")["receptor_ph"] == 7.0


def test_interactive_refused_store_keeps_the_manifest(tmp_path, capsys):
    project = _project(tmp_path)
    before = (project / "project_manifest.json").read_text(encoding="utf-8")
    # The user declines the store after reading the stale warning.
    fake = _FakeQuestionary(texts=["1IEP", "6.5"], selects=[""], confirms=[False])
    interactive._set_protonation_override_flow(fake, project, "receptor")
    assert (project / "project_manifest.json").read_text(encoding="utf-8") == before
    assert "The override was not stored." in capsys.readouterr().out


def test_interactive_stage_shows_levels_and_stale_entities(tmp_path, capsys):
    project = _project(tmp_path)
    # Pick "back" at the first menu prompt: the stage prints the view, then the menu returns control.
    fake = _FakeQuestionary(selects=["back"])
    with pytest.raises(interactive._PromptNavigation):
        interactive._protonation_stage(fake, project)
    out = capsys.readouterr().out
    assert "Project default (level: project)" in out
    assert "from receptor" in out and "2ABC" in out
    assert "policy as_input (from ligand)" in out
    assert "Pair protonation status" in out and "pair_ph_mismatch=1" in out
    assert "2ABC: " in out and STALE in out
    assert "Clearing an override is not available yet" in out


# ----------------------------------------------------------------------------- web UI


@pytest.fixture()
def client(monkeypatch, tmp_path):
    monkeypatch.setenv("OMNIDOCK_WEBUI_HOME", str(tmp_path / "webui_home"))
    import webui.config
    import webui.registry
    import webui.jobs
    import webui.app

    for module in (webui.config, webui.registry, webui.jobs, webui.app):
        importlib.reload(module)
    webui.config.ensure_app_dirs()
    app = webui.app.create_app()
    app.config.update(TESTING=True)
    test_client = app.test_client()
    yield test_client
    import webui.jobs as J

    for job in J.list_jobs():
        if not job.is_terminal:
            try:
                J.cancel(job.job_id)
            except Exception:
                pass


def _register(client, project: Path) -> str:
    return client.post("/api/projects", json={"path": str(project)}).get_json()["id"]


def test_webui_api_and_page_show_the_levels(client, tmp_path):
    project = _project(tmp_path)
    pid = _register(client, project)

    body = client.get(f"/api/projects/{pid}/protonation").get_json()
    receptors = {row["name"]: row for row in body["receptors"]}
    assert receptors["2ABC"]["field_levels"]["receptor_ph"] == "receptor"
    assert receptors["1IEP"]["field_levels"]["receptor_ph"] == "project"
    assert {(item["kind"], item["name"]) for item in body["stale"]} == {("receptor", "2ABC"), ("ligand", "STI")}
    assert body["pair"]["status_counts"]["pair_ph_mismatch"] == 1

    page = client.get(f"/project/{pid}/protonation")
    assert page.status_code == 200
    html = page.get_data(as_text=True)
    assert "Project default" in html and "receptor" in html and "project" in html
    assert 'id="receptor-table"' in html and 'id="ligand-table"' in html
    assert "as_input" in html
    assert STALE in html, "the stale reason must appear on the page"
    assert "Stale prepared entities" in html
    assert "pair_ph_mismatch" in html
    assert "stale_preparation_policy_changed" in html and "Changing a policy makes" in html


def test_webui_read_routes_never_write_into_the_project(client, tmp_path):
    """The read-only contract (FR-035): the view and the page leave the project tree unchanged."""
    project = _project(tmp_path)
    pid = _register(client, project)
    before = _snapshot(project)
    assert client.get(f"/api/projects/{pid}/protonation").status_code == 200
    assert client.get(f"/project/{pid}/protonation").status_code == 200
    assert _snapshot(project) == before


def test_webui_unavailable_project_is_reported_not_written(client, tmp_path):
    project = tmp_path / "gone"
    project.mkdir()
    pid = _register(client, project)
    project.rmdir()
    assert client.get(f"/api/projects/{pid}/protonation").status_code == 409
    page = client.get(f"/project/{pid}/protonation")
    assert page.status_code == 200
    assert "no longer available" in page.get_data(as_text=True)


def test_webui_override_form_is_built_from_the_protonation_target(client, tmp_path):
    project = _project(tmp_path)
    pid = _register(client, project)
    schema = client.get("/api/targets/workflow.protonation-policy/form").get_json()
    names = {f["name"] for f in schema["fields"]}
    assert {"receptor", "receptor_ph", "receptor_force_field", "ligand", "ligand_policy", "ligand_ph",
            "ligand_smiles", "ligand_net_charge", "ligand_state_map", "replace"} <= names
    html = client.get(f"/project/{pid}/protonation").get_data(as_text=True)
    assert "workflow.protonation-policy" in html


def test_webui_set_action_is_the_cli_command_and_the_preview_matches(client, tmp_path):
    project = _project(tmp_path)
    pid = _register(client, project)
    values = {"receptor": "1IEP", "receptor_ph": "6.5"}
    preview = client.post(f"/api/projects/{pid}/preview", json={"target": "workflow.protonation-policy", "values": values})
    argv = preview.get_json()["argv"]
    assert argv[1:4] == ["main.py", "workflow", "protonation-policy"]
    assert argv[4] == "--project-dir" and Path(argv[5]).resolve() == project.resolve()
    assert argv[6:] == ["--receptor", "1IEP", "--receptor-ph", "6.5"]


def test_webui_scope_validation_refuses_receptor_and_ligand_together(client, tmp_path):
    import webui.targets as T

    errors = T.validate("workflow.protonation-policy", {"receptor": "1IEP", "receptor_ph": "6.5", "ligand": "STI", "ligand_policy": "as_input"})
    assert "ligand" in errors
    assert "receptor" not in T.validate("workflow.protonation-policy", {"receptor": "1IEP", "receptor_ph": "6.5"})
    assert "ligand_policy" in T.validate("workflow.protonation-policy", {"receptor": "1IEP", "receptor_ph": "6.5", "ligand_policy": "as_input"})
    assert "receptor" in T.validate("workflow.protonation-policy", {"receptor": "1IEP"})
    assert "ligand" in T.validate("workflow.protonation-policy", {"ligand": "STI"})


def _wait_for_job(client, job_id: str, timeout: float = 180.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = client.get(f"/api/jobs/{job_id}").get_json()
        if job["status"] in ("completed", "failed", "cancelled"):
            return job
        time.sleep(0.5)
    raise AssertionError(f"job {job_id} did not finish")


def test_webui_set_runs_the_cli_and_the_conflict_message_reaches_the_job_log(client, tmp_path):
    """End to end: the form's job runs `workflow protonation-policy`, which hits the same conflict check."""
    project = _project(tmp_path)
    pid = _register(client, project)
    res = client.post(f"/api/projects/{pid}/jobs", json={
        "target": "workflow.protonation-policy",
        "values": {"receptor": "2ABC", "receptor_ph": "7.0"},
    })
    assert res.status_code == 202, res.get_data(as_text=True)
    job = _wait_for_job(client, res.get_json()["job_id"])
    assert job["status"] == "failed" and job["exit_code"] == 2
    log = client.get(f"/api/jobs/{job['job_id']}/log").get_json()["text"]
    assert "protonation_policy_conflict" in log
    assert resolve_entity_protonation(project, "receptor", "2ABC")["receptor_ph"] == 6.0


def test_webui_set_stores_an_override_through_the_cli_job(client, tmp_path):
    project = _project(tmp_path)
    pid = _register(client, project)
    res = client.post(f"/api/projects/{pid}/jobs", json={
        "target": "workflow.protonation-policy",
        "values": {"receptor": "1IEP", "receptor_ph": "7.0"},
    })
    assert res.status_code == 202
    job = _wait_for_job(client, res.get_json()["job_id"])
    assert job["status"] == "completed" and job["exit_code"] == 0
    assert resolve_entity_protonation(project, "receptor", "1IEP")["receptor_ph"] == 7.0
    view = client.get(f"/api/projects/{pid}/protonation").get_json()
    assert {row["name"]: row for row in view["receptors"]}["1IEP"]["field_levels"]["receptor_ph"] == "receptor"


# ----------------------------------------------------------------------------- CLI conflict (what the job runs)


def test_cli_conflict_names_the_stored_value(tmp_path, capsys):
    project = _project(tmp_path)
    code = cli_main(["workflow", "protonation-policy", "--project-dir", str(project),
                     "--receptor", "2ABC", "--receptor-ph", "7.0"])
    assert code == 2
    err = capsys.readouterr().err
    assert "protonation_policy_conflict" in err and "6.0" in err


def test_conflict_is_a_value_error_for_callers(tmp_path):
    project = _project(tmp_path)
    with pytest.raises(ProtonationPolicyConflict):
        set_protonation_override(project, kind="receptor", name="2ABC", receptor_ph=7.0)
