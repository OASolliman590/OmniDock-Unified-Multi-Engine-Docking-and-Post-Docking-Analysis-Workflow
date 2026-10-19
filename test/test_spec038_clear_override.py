"""Spec 038 (coordinator item): clear a stored protonation override.

- ``clear_protonation_override`` removes ``receptors.<key>`` or ``ligands.<key>`` and the entity then
  inherits the project default.
- The effective policy hash changes when the override differed from the default, so the prepared-reuse
  guard (``reuse_block_reason``) reports a prepared entity as ``policy_changed``.
- ``workflow protonation-policy --clear-receptor NAME`` / ``--clear-ligand NAME`` drive the same function.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from docking.preparation.prepared_reuse import reuse_block_reason
from docking.project_layout import (
    ProtonationPolicyError,
    clear_protonation_override,
    load_protonation_policy,
    resolve_entity_protonation,
    set_protonation_override,
    set_protonation_policy,
)
from workflow.cli import build_parser, main as workflow_main


def _project(tmp_path: Path) -> Path:
    project = tmp_path / "project"
    project.mkdir()
    (project / "project_manifest.json").write_text(
        json.dumps(
            {
                "project_name": "spec038_clear",
                "project_root": str(project),
                "layout_profile": "canonical",
                "engines": ["vina"],
                "favorite_engine": "vina",
            }
        ),
        encoding="utf-8",
    )
    set_protonation_policy(
        project,
        receptor_ph=7.4,
        receptor_force_field="AMBER",
        ligand_policy="ph_model",
        ligand_ph=7.4,
    )
    return project


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_clear_receptor_override_restores_the_project_default(tmp_path):
    project = _project(tmp_path)
    set_protonation_override(project, kind="receptor", name="1IEP", receptor_ph=6.0)
    assert resolve_entity_protonation(project, "receptor", "1IEP")["policy_level"] == "receptor"

    result = clear_protonation_override(project, kind="receptor", name="1IEP")

    assert result["removed"]["receptor_ph"] == 6.0
    after = resolve_entity_protonation(project, "receptor", "1IEP")
    assert after["policy_level"] == "project"
    assert after["effective"]["receptor_ph"] == 7.4
    assert "1IEP" not in (load_protonation_policy(project).get("receptors") or {})
    assert result["policy_hash_after"] == after["policy_hash"]
    assert result["policy_hash_before"] != result["policy_hash_after"]


def test_clear_ligand_override_accepts_a_file_name(tmp_path):
    project = _project(tmp_path)
    set_protonation_override(project, kind="ligand", name="STI", ligand_ph=6.5)

    result = clear_protonation_override(project, kind="ligand", name="STI.sdf")

    assert result["entity"] == "STI"
    assert result["removed"]["ligand_ph"] == 6.5
    assert resolve_entity_protonation(project, "ligand", "STI")["policy_level"] == "project"
    assert "STI" not in (load_protonation_policy(project).get("ligands") or {})


def test_clearing_an_override_equal_to_the_default_changes_no_hash(tmp_path):
    project = _project(tmp_path)
    set_protonation_override(project, kind="receptor", name="1IEP", receptor_ph=7.4)

    result = clear_protonation_override(project, kind="receptor", name="1IEP")

    assert result["policy_hash_before"] == result["policy_hash_after"]


def test_clearing_a_missing_override_is_refused_and_changes_nothing(tmp_path):
    project = _project(tmp_path)
    manifest_before = (project / "project_manifest.json").read_text(encoding="utf-8")

    with pytest.raises(ProtonationPolicyError) as excinfo:
        clear_protonation_override(project, kind="receptor", name="NOPE")

    assert excinfo.value.reason == "protonation_override_missing"
    assert (project / "project_manifest.json").read_text(encoding="utf-8") == manifest_before


def test_unknown_kind_and_empty_name_are_refused(tmp_path):
    project = _project(tmp_path)
    with pytest.raises(ValueError, match="unknown protonation entity kind"):
        clear_protonation_override(project, kind="protein", name="1IEP")
    with pytest.raises(ValueError, match="protonation_override_name_required"):
        clear_protonation_override(project, kind="receptor", name="  ")


def test_cleared_entity_prepared_under_the_override_is_stale(tmp_path):
    """The hash guard: an output prepared under the override hash is refused once the override is cleared."""
    project = _project(tmp_path)
    set_protonation_override(project, kind="receptor", name="1IEP", receptor_ph=6.0)
    override_hash = resolve_entity_protonation(project, "receptor", "1IEP")["policy_hash"]

    source = tmp_path / "1IEP.pdb"
    output = tmp_path / "1IEP.pdbqt"
    source.write_text("ATOM\n", encoding="utf-8")
    output.write_text("PDBQT\n", encoding="utf-8")
    provenance = {
        "policy_hash": override_hash,
        "request_fingerprint": "fp",
        "input_sha256": _sha(source),
        "output_sha256": _sha(output),
    }
    # Still under the override: the record matches the current hash.
    assert reuse_block_reason(
        provenance, output=output, input_file=source, policy_hash=override_hash, fingerprint="fp", force=False
    ) == ""

    clear_protonation_override(project, kind="receptor", name="1IEP")
    current = resolve_entity_protonation(project, "receptor", "1IEP")["policy_hash"]
    assert current != override_hash
    assert reuse_block_reason(
        provenance, output=output, input_file=source, policy_hash=current, fingerprint="fp", force=False
    ) == "policy_changed"


def test_parser_accepts_the_clear_flags():
    parser = build_parser()
    args = parser.parse_args(["workflow", "protonation-policy", "--project-dir", "p", "--clear-receptor", "1IEP"])
    assert args.clear_receptor == "1IEP" and args.clear_ligand is None
    args = parser.parse_args(["workflow", "protonation-policy", "--project-dir", "p", "--clear-ligand", "STI"])
    assert args.clear_ligand == "STI" and args.clear_receptor is None


def test_cli_clear_receptor_and_clear_ligand_work_end_to_end(tmp_path, capsys):
    project = _project(tmp_path)
    set_protonation_override(project, kind="receptor", name="1IEP", receptor_ph=6.0)
    set_protonation_override(project, kind="ligand", name="STI", ligand_ph=6.5)

    assert workflow_main(["workflow", "protonation-policy", "--project-dir", str(project), "--clear-receptor", "1IEP"]) == 0
    assert "cleared" in capsys.readouterr().out
    assert resolve_entity_protonation(project, "receptor", "1IEP")["policy_level"] == "project"
    assert resolve_entity_protonation(project, "ligand", "STI")["policy_level"] == "ligand"

    assert workflow_main(["workflow", "protonation-policy", "--project-dir", str(project), "--clear-ligand", "STI"]) == 0
    assert resolve_entity_protonation(project, "ligand", "STI")["policy_level"] == "project"


def test_cli_refuses_combined_or_missing_clears(tmp_path, capsys):
    project = _project(tmp_path)
    set_protonation_override(project, kind="receptor", name="1IEP", receptor_ph=6.0)
    base = ["workflow", "protonation-policy", "--project-dir", str(project)]

    assert workflow_main(base + ["--clear-receptor", "1IEP", "--receptor-ph", "7.0"]) == 2
    assert "cannot be combined" in capsys.readouterr().err
    assert workflow_main(base + ["--clear-receptor", "1IEP", "--clear-ligand", "STI"]) == 2
    assert workflow_main(base + ["--clear-ligand", "NOPE"]) == 2
    assert "protonation_override_missing" in capsys.readouterr().err
    # The stored override is untouched by the refused calls.
    assert resolve_entity_protonation(project, "receptor", "1IEP")["policy_level"] == "receptor"
