"""Configuration seam tests that do not require a chemistry runtime."""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from docking.models import (
    DEFAULT_PREPARATION_PH,
    LIGAND_PREPARATION_PROFILES,
    normalize_ligand_preparation_profile,
    validate_ligand_preparation_profile,
    validate_preparation_ph,
)
from docking.preparation import ligand_preparation as preparation
from docking.preparation.ligand_identity import sha256_file


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf"), True, False])
def test_preparation_ph_rejects_nonfinite_and_boolean_values(value):
    valid, normalized, error = validate_preparation_ph(value)
    assert valid is False
    assert normalized == DEFAULT_PREPARATION_PH
    assert error


def test_unknown_profile_is_invalid_at_strict_validation_seam():
    result = validate_ligand_preparation_profile("made-up-profile", ["vina"])
    assert result.is_valid is False
    assert any("Unknown ligand preparation profile" in error for error in result.errors)
    assert all(profile in result.errors[0] for profile in LIGAND_PREPARATION_PROFILES)


@pytest.mark.parametrize(
    ("alias", "canonical"),
    [
        ("", "engine_aware_full"),
        ("auto", "engine_aware_full"),
        ("default", "engine_aware_full"),
        ("obabel", "openbabel_only"),
        ("openbabel", "openbabel_only"),
        ("meeko", "meeko_only"),
        ("adt", "autodocktools_only"),
    ],
)
def test_documented_profile_aliases_remain_valid(alias, canonical):
    assert normalize_ligand_preparation_profile(alias) == canonical
    result = validate_ligand_preparation_profile(alias, [])
    assert result.is_valid is True
    assert result.requested_profile == canonical


def test_direct_environment_seam_does_not_default_invalid_profile_or_ph(monkeypatch):
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "made-up-profile")
    with pytest.raises(ValueError, match="Unknown ligand preparation profile"):
        preparation._selected_profile()

    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PROFILE", "meeko")
    monkeypatch.setenv("PDBWIZARD_LIGAND_PREP_PH", "nan")
    with pytest.raises(ValueError, match="finite"):
        preparation._selected_ph()


def test_exact_output_validation_requires_is_valid_for_legacy_profiles(tmp_path, monkeypatch):
    source = tmp_path / "source.sdf"
    normalized = tmp_path / "normalized.sdf"
    output = tmp_path / "output.pdbqt"
    for path, payload in (
        (source, b"source"),
        (normalized, b"normalized"),
        (output, b"output"),
    ):
        path.write_bytes(payload)

    monkeypatch.setattr(
        preparation,
        "_validate_prepared_pdbqt_contract",
        lambda _path: {"is_valid": True, "errors": [], "warnings": []},
    )
    summary = {
        "input_file": str(source),
        "output_file": str(output),
        "input_sha256": sha256_file(source),
        "output_sha256": sha256_file(output),
        "requested_profile": "openbabel_autodocktools",
        "requested_backend": "openbabel_autodocktools",
        "effective_profile": "openbabel_autodocktools",
        "selected_engines": [],
        "preparation_method": "openbabel_then_autodocktools",
        "protonation_ph": 7.4,
        "charge_state_status": "backend_assigned_autodocktools",
        "normalization": {
            "normalization_backend": "openbabel",
            "protonation_applied": True,
            "normalized_sdf": str(normalized),
            "normalized_sdf_sha256": sha256_file(normalized),
            "identity_ledger": {"is_valid": True},
            "normalization_provenance": {"procedure": "test"},
        },
        "output_chemistry_validation": {
            "status": "exact",
            "is_valid": False,
            "scope": "test",
        },
    }
    result = preparation.validate_ligand_preparation_output_contract(summary)
    assert result["is_valid"] is False
    assert any("is_valid=true" in error for error in result["errors"])
