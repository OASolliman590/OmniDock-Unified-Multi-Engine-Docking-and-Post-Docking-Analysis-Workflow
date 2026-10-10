"""Spec 036 R5 (docking parameters from evidence) and R6 software fixes.

Expected values are derived here from closed-form geometry, not from the code under test.
No real docking is run: dry runs and a fake Vina executable only.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from openpyxl import Workbook

import docking.cli as docking_cli
import docking.preparation.excel_sites as excel_sites
import docking.preparation.pairlist_builder as pairlist_builder
from docking import box_policy
from docking.box_policy import (
    BOX_METHOD_RG,
    BOX_METHOD_USER_FIXED,
    CONTAINMENT_CONTAINED,
    CONTAINMENT_NOT_CONTAINED,
    EDGE_WARNING_ANGSTROM,
    RG_EDGE_FACTOR,
    box_warnings,
    compute_ligand_box,
    containment_status,
)
from docking.models import PairlistRow, classify_preparation_status, validate_ligand_preparation_profile
from docking.parameter_schema import (
    DEFAULT_PARAMETER_PRESET,
    apply_schema_to_runtime,
    resolve_parameter_schema,
    transform_pairlist_rows,
    validate_parameter_schema,
)
from docking.runners.base import build_pair_index, replicate_job_tag, split_replicate_stem
from docking.runners.smina import SminaRunner
from docking.runners.vina import VinaRunner
from workflow import execution


REPO_ROOT = Path(__file__).resolve().parents[1]
STI_SDF = REPO_ROOT / "test" / "fixtures" / "spec032" / "phase2_1iep_inputs" / "STI_1IEP_A201.sdf"
EXCEL_CENTRE = (15.61389196241224, 53.38013499491924, 15.45483767019736)  # 1IEP STI centre in the scratch project
LIGAND_NAME = "1IEP_ligand_STI_A_201.pdbqt"
RECEPTOR_NAME = "1IEP_receptor.pdbqt"


# ----------------------------------------------------------------------------- fixtures


def _pdbqt_text(atoms):
    """atoms: (x, y, z, type). Columns follow the PDB layout that the box reader slices."""
    lines = []
    for serial, (x, y, z, atom_type) in enumerate(atoms, start=1):
        head = f"ATOM  {serial:5d} {'C':<4s} LIG A   1"
        body = f"{x:8.3f}{y:8.3f}{z:8.3f}" + "  1.00  0.00    0.000"
        lines.append((head.ljust(30) + body).ljust(77) + f"{atom_type:<2s}")
    return "\n".join(lines) + "\n"


def _hexagon_atoms(centre=(0.0, 0.0, 0.0)):
    """Six heavy carbons at 1.39 A from the centre, six hydrogens at 2.48 A (excluded from Rg)."""
    cx, cy, cz = centre
    atoms = []
    for k in range(6):
        angle = math.radians(60.0 * k)
        atoms.append((cx + 1.39 * math.cos(angle), cy + 1.39 * math.sin(angle), cz, "C"))
    for k in range(6):
        angle = math.radians(60.0 * k)
        atoms.append((cx + 2.48 * math.cos(angle), cy + 2.48 * math.sin(angle), cz, "HD"))
    return atoms


def _chain_atoms(n=10, spacing=4.0):
    return [(spacing * k, 0.0, 0.0, "C") for k in range(n)]


def _write_excel(path: Path, *, method=None, selected=LIGAND_NAME, centre=(0.0, 0.0, 0.0)) -> Path:
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Summary"
    sheet.append(["PDB_ID", "Property", "Value"])
    rows = [
        ("1IEP", "selected_ligand", selected),
        ("1IEP", "active_site_center_x", centre[0]),
        ("1IEP", "active_site_center_y", centre[1]),
        ("1IEP", "active_site_center_z", centre[2]),
    ]
    if method is not None:
        rows.append(("1IEP", "binding_site_center_method", method))
    for row in rows:
        sheet.append(list(row))
    workbook.save(path)
    return path


@pytest.fixture
def pair_patches(monkeypatch):
    monkeypatch.setattr(
        pairlist_builder,
        "ensure_project_alias_files",
        lambda **_kwargs: {"protein_alias_file": "", "ligand_alias_file": ""},
    )
    monkeypatch.setattr(pairlist_builder, "build_protein_alias_lookup", lambda *_a, **_k: {})
    monkeypatch.setattr(pairlist_builder, "build_ligand_alias_lookup", lambda *_a, **_k: {})
    return monkeypatch


def _build_pairlist(tmp_path: Path, *, ligand_atoms, box_size=None, method="binding_site_center_v1",
                    centre=(0.0, 0.0, 0.0), reference_atoms=None, reference_name=LIGAND_NAME):
    proteins = tmp_path / "prepared_proteins"
    ligands = tmp_path / "prepared_ligands"
    proteins.mkdir(exist_ok=True)
    ligands.mkdir(exist_ok=True)
    (proteins / RECEPTOR_NAME).write_text("ATOM\n", encoding="utf-8")
    (ligands / LIGAND_NAME).write_text(_pdbqt_text(ligand_atoms), encoding="utf-8")
    raw_dir = None
    if reference_atoms is not None:
        raw_dir = tmp_path / "raw_ligands"
        raw_dir.mkdir(exist_ok=True)
        (raw_dir / reference_name).write_text(_pdbqt_text(reference_atoms), encoding="utf-8")
    excel = _write_excel(tmp_path / "multi_pdb_analysis.xlsx", method=method, centre=centre)
    summary = pairlist_builder.build_pairlists(
        project_root=tmp_path / "project",
        prepared_proteins=proteins,
        prepared_ligands=ligands,
        excel_path=excel,
        mode="cocrystal_only",
        default_box_size=box_size,
        raw_ligands=raw_dir,
    )
    return summary, pd.read_csv(summary["pairlist_file"])


# ----------------------------------------------------------------------------- R5a box


def test_rg_scaled_edge_on_known_geometry(tmp_path: Path) -> None:
    ligand = tmp_path / "hexagon.pdbqt"
    ligand.write_text(_pdbqt_text(_hexagon_atoms()), encoding="utf-8")
    box = compute_ligand_box(ligand)
    # Heavy atoms only: six carbons at 1.39 A from the centroid -> Rg = 1.39 A exactly.
    expected_rg = 1.39
    expected_edge = 2.9 * expected_rg  # 4.031 A
    assert box.box_method == BOX_METHOD_RG
    assert box.ligand_heavy_atom_count == 6
    assert box.ligand_rg_angstrom == pytest.approx(expected_rg, abs=1e-3)
    assert box.edge_angstrom == pytest.approx(expected_edge, abs=1e-3)
    assert RG_EDGE_FACTOR == 2.9


def test_rg_of_chain_matches_closed_form(tmp_path: Path) -> None:
    ligand = tmp_path / "chain.pdbqt"
    ligand.write_text(_pdbqt_text(_chain_atoms(n=10, spacing=4.0)), encoding="utf-8")
    box = compute_ligand_box(ligand)
    # Collinear equally spaced atoms: Rg = spacing * sqrt((n^2 - 1) / 12).
    expected_rg = 4.0 * math.sqrt((10**2 - 1) / 12.0)
    assert box.ligand_rg_angstrom == pytest.approx(expected_rg, abs=1e-9)
    assert box.edge_angstrom == pytest.approx(2.9 * expected_rg, abs=1e-9)
    assert box.edge_angstrom > EDGE_WARNING_ANGSTROM


def test_sti_fixture_rg_matches_independent_parse() -> None:
    # Independent reader: fixed-width V2000 atom block, parsed here without the box module.
    lines = STI_SDF.read_text(encoding="utf-8").splitlines()
    natoms = int(lines[3][:3])
    atoms = [line.split() for line in lines[4 : 4 + natoms]]
    heavy = np.array([[float(a[0]), float(a[1]), float(a[2])] for a in atoms if a[3] != "H"])
    expected_rg = float(np.sqrt(np.mean(np.sum((heavy - heavy.mean(axis=0)) ** 2, axis=1))))
    assert heavy.shape[0] == 37
    box = compute_ligand_box(STI_SDF)
    assert box.ligand_rg_angstrom == pytest.approx(expected_rg, abs=1e-9)
    assert box.edge_angstrom == pytest.approx(2.9 * expected_rg, abs=1e-9)
    # Frame check: the Excel centre equals the crystal ligand heavy-atom centroid.
    assert np.allclose(heavy.mean(axis=0), EXCEL_CENTRE, atol=5e-3)


def test_user_fixed_is_recorded(tmp_path: Path) -> None:
    ligand = tmp_path / "hexagon.pdbqt"
    ligand.write_text(_pdbqt_text(_hexagon_atoms()), encoding="utf-8")
    box = compute_ligand_box(ligand, fixed_edge=24.0)
    assert box.box_method == BOX_METHOD_USER_FIXED
    assert box.edge_angstrom == 24.0
    assert box.ligand_rg_angstrom == pytest.approx(1.39, abs=1e-3)


def test_pairlist_records_method_rg_and_edge_per_pair(pair_patches, tmp_path: Path) -> None:
    _summary, pairlist = _build_pairlist(tmp_path, ligand_atoms=_hexagon_atoms())
    row = pairlist.iloc[0]
    assert row["box_method"] == BOX_METHOD_RG
    assert float(row["ligand_rg_angstrom"]) == pytest.approx(1.39, abs=1e-3)
    assert float(row["edge_angstrom"]) == pytest.approx(2.9 * 1.39, abs=1e-3)
    assert float(row["size_x"]) == pytest.approx(2.9 * 1.39, abs=1e-3)
    # Old columns are kept.
    for column in ("receptor", "site_id", "ligand", "center_x", "size_x", "pdb_id", "pair_source"):
        assert column in pairlist.columns


def test_pairlist_user_fixed_box_is_recorded(pair_patches, tmp_path: Path) -> None:
    _summary, pairlist = _build_pairlist(tmp_path, ligand_atoms=_hexagon_atoms(), box_size=20.0)
    row = pairlist.iloc[0]
    assert row["box_method"] == BOX_METHOD_USER_FIXED
    assert float(row["edge_angstrom"]) == 20.0
    assert float(row["size_x"]) == 20.0


def test_containment_warning_fires_when_reference_is_outside_box(pair_patches, tmp_path: Path) -> None:
    # Reference ligand sits 10 A away from the site centre; the 4.03 A box cannot hold it.
    summary, pairlist = _build_pairlist(
        tmp_path,
        ligand_atoms=_hexagon_atoms(),
        centre=(0.0, 0.0, 0.0),
        reference_atoms=_hexagon_atoms(centre=(10.0, 0.0, 0.0)),
    )
    row = pairlist.iloc[0]
    assert row["box_containment_status"] == CONTAINMENT_NOT_CONTAINED
    assert "does not contain all reference heavy atoms" in str(row["box_warnings"])
    assert any("does not contain all reference heavy atoms" in message for message in summary["warnings"])


def test_no_containment_warning_when_reference_is_inside_box(pair_patches, tmp_path: Path) -> None:
    summary, pairlist = _build_pairlist(
        tmp_path,
        ligand_atoms=_hexagon_atoms(),
        centre=(0.0, 0.0, 0.0),
        reference_atoms=_hexagon_atoms(centre=(0.0, 0.0, 0.0)),
    )
    assert pairlist.iloc[0]["box_containment_status"] == CONTAINMENT_CONTAINED
    assert not any("reference heavy atoms" in message for message in summary["warnings"])


def test_containment_status_helper_counts_atoms_outside() -> None:
    reference = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    result = containment_status((0.0, 0.0, 0.0), 4.0, reference)  # half-edge 2 A
    assert result["status"] == CONTAINMENT_NOT_CONTAINED
    assert result["outside_count"] == 1
    assert result["reference_atom_count"] == 2
    assert containment_status((0.0, 0.0, 0.0), 4.0, None)["status"] == "not_evaluated_no_reference"


def test_edge_over_30_angstrom_warning_fires_for_rg_box(pair_patches, tmp_path: Path) -> None:
    # Chain of 10 carbons at 4 A spacing: Rg = 11.489 A, edge = 33.32 A > 30 A.
    summary, pairlist = _build_pairlist(tmp_path, ligand_atoms=_chain_atoms(n=10, spacing=4.0))
    row = pairlist.iloc[0]
    assert float(row["edge_angstrom"]) > 30.0
    assert "exceeds 30 A" in str(row["box_warnings"])
    assert any("exceeds 30 A" in message for message in summary["warnings"])


def test_edge_over_30_angstrom_warning_fires_for_user_fixed(pair_patches, tmp_path: Path) -> None:
    summary, pairlist = _build_pairlist(tmp_path, ligand_atoms=_hexagon_atoms(), box_size=35.0)
    assert pairlist.iloc[0]["box_method"] == BOX_METHOD_USER_FIXED
    assert any("exceeds 30 A" in message for message in summary["warnings"])


def test_box_warnings_helper_levels() -> None:
    not_contained = {"status": CONTAINMENT_NOT_CONTAINED, "outside_count": 2, "reference_atom_count": 37}
    messages = box_warnings(pair_label="p", edge_angstrom=31.0, containment=not_contained, reference_label="STI")
    assert len(messages) == 2
    assert box_warnings(pair_label="p", edge_angstrom=20.0, containment={"status": CONTAINMENT_CONTAINED}) == []


def test_each_ligand_gets_its_own_edge(pair_patches, tmp_path: Path) -> None:
    proteins = tmp_path / "prepared_proteins"
    ligands = tmp_path / "prepared_ligands"
    proteins.mkdir()
    ligands.mkdir()
    (proteins / RECEPTOR_NAME).write_text("ATOM\n", encoding="utf-8")
    (ligands / LIGAND_NAME).write_text(_pdbqt_text(_hexagon_atoms()), encoding="utf-8")
    (ligands / "1IEP_ligand_CHAIN.pdbqt").write_text(_pdbqt_text(_chain_atoms(n=10, spacing=4.0)), encoding="utf-8")
    excel = _write_excel(tmp_path / "multi_pdb_analysis.xlsx", method="binding_site_center_v1")
    summary = pairlist_builder.build_pairlists(
        project_root=tmp_path / "project",
        prepared_proteins=proteins,
        prepared_ligands=ligands,
        excel_path=excel,
        mode="curated_per_protein",
        curated_mapping={RECEPTOR_NAME: [LIGAND_NAME, "1IEP_ligand_CHAIN.pdbqt"]},
    )
    pairlist = pd.read_csv(summary["pairlist_file"]).set_index("ligand")
    assert float(pairlist.loc[LIGAND_NAME, "edge_angstrom"]) == pytest.approx(2.9 * 1.39, abs=1e-3)
    assert float(pairlist.loc["1IEP_ligand_CHAIN.pdbqt", "edge_angstrom"]) == pytest.approx(
        2.9 * 4.0 * math.sqrt(99 / 12.0), abs=1e-9
    )


def test_unreadable_conformer_is_refused_not_defaulted(tmp_path: Path) -> None:
    empty = tmp_path / "empty.pdbqt"
    empty.write_text("ATOM\n", encoding="utf-8")
    with pytest.raises(box_policy.BoxError):
        compute_ligand_box(empty)


# ----------------------------------------------------------------------------- dock-run level


def _dock_project(tmp_path: Path, monkeypatch) -> Path:
    root = tmp_path / "project"
    execution.run_workflow_init(str(root), layout_profile="canonical")
    from docking.project_layout import ensure_project_layout

    layout = ensure_project_layout(root, "canonical")
    (layout["prepared_proteins"] / RECEPTOR_NAME).write_text("ATOM\n", encoding="utf-8")
    (layout["prepared_ligands"] / LIGAND_NAME).write_text(_pdbqt_text(_hexagon_atoms()), encoding="utf-8")
    # `prep project` materialises shared receptors/ligands; preflight checks those directories.
    (layout["receptors"] / RECEPTOR_NAME).write_text("ATOM\n", encoding="utf-8")
    (layout["ligands"] / LIGAND_NAME).write_text(_pdbqt_text(_hexagon_atoms()), encoding="utf-8")
    excel = _write_excel(tmp_path / "multi_pdb_analysis.xlsx", method="binding_site_center_v1")
    monkeypatch.setattr(
        pairlist_builder,
        "ensure_project_alias_files",
        lambda **_kwargs: {"protein_alias_file": "", "ligand_alias_file": ""},
    )
    monkeypatch.setattr(pairlist_builder, "build_protein_alias_lookup", lambda *_a, **_k: {})
    monkeypatch.setattr(pairlist_builder, "build_ligand_alias_lookup", lambda *_a, **_k: {})
    pairlist_builder.build_pairlists(
        project_root=root,
        prepared_proteins=layout["prepared_proteins"],
        prepared_ligands=layout["prepared_ligands"],
        excel_path=excel,
        mode="cocrystal_only",
        layout_profile="canonical",
    )
    return root


def _dock_args(root: Path, *extra: str):
    return [
        "--project-dir",
        str(root),
        "--engines",
        "vina",
        "--vina-binary",
        "vina",
        "--dry-run",
        "--no-ligand-qc-gate",
        "--no-receptor-qc-gate",
        *extra,
    ]


def _run_manifest(root: Path, engine: str = "vina") -> dict:
    candidates = [
        root / "engines" / engine / "run_manifest.json",
        root / "3-Docking" / engine / "run_manifest.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return json.loads(candidate.read_text(encoding="utf-8"))
    raise AssertionError(f"no run manifest found for {engine}: {candidates}")


def test_missing_seed_fails_dock_run(monkeypatch, tmp_path: Path, pair_patches) -> None:
    root = _dock_project(tmp_path, monkeypatch)
    assert docking_cli.dock_main(_dock_args(root)) != 0


def test_replicates_produce_seeds_and_paths_in_dry_run_manifest(monkeypatch, tmp_path: Path, pair_patches) -> None:
    root = _dock_project(tmp_path, monkeypatch)
    assert docking_cli.dock_main(_dock_args(root, "--seed", "100", "--replicates", "3")) == 0
    manifest = _run_manifest(root)
    assert manifest["effective"]["replicates"] == 3
    assert manifest["effective"]["seeds"] == [100, 101, 102]
    jobs = manifest["jobs"]
    assert [job["seed"] for job in jobs] == [100, 101, 102]
    assert [job["replicate_id"] for job in jobs] == [1, 2, 3]
    assert all(job["pair_tag"] == jobs[0]["pair_tag"] for job in jobs)
    for job in jobs:
        assert [item["seed"] for item in job["replicates"]] == [100, 101, 102]
        assert [item["output_path"] for item in job["replicates"]] == [
            item["output_path"] for item in jobs[0]["replicates"]
        ]
        assert job["box_method"] == BOX_METHOD_RG
        assert job["edge_angstrom"] == pytest.approx(2.9 * 1.39, abs=1e-3)
    names = [Path(job["pose_file"]).name for job in jobs]
    assert all(name.endswith(f"__rep0{index}.pdbqt") for index, name in enumerate(names, start=1))
    # Each seed is on the command line of its own replicate.
    for job in jobs:
        index = job["command"].index("--seed")
        assert job["command"][index + 1] == str(job["seed"])


def test_single_replicate_is_recorded_as_one(monkeypatch, tmp_path: Path, pair_patches) -> None:
    root = _dock_project(tmp_path, monkeypatch)
    assert docking_cli.dock_main(_dock_args(root, "--seed", "7", "--replicates", "1")) == 0
    manifest = _run_manifest(root)
    assert manifest["effective"]["replicates"] == 1
    assert [job["seed"] for job in manifest["jobs"]] == [7]


def test_exhaustiveness_defaults_to_32_and_is_recorded(monkeypatch, tmp_path: Path, pair_patches) -> None:
    root = _dock_project(tmp_path, monkeypatch)
    assert docking_cli.dock_main(_dock_args(root, "--seed", "1")) == 0
    manifest = _run_manifest(root)
    assert manifest["effective"]["exhaustiveness"] == 32
    assert manifest["jobs"][0]["exhaustiveness"] == 32
    project_manifest = json.loads((root / "project_manifest.json").read_text(encoding="utf-8"))
    assert project_manifest["engine_effective_parameters"]["vina"]["exhaustiveness"] == 32
    assert project_manifest["parameter_schema"]["common"]["exhaustiveness"] == 32


def test_default_parameter_preset_is_exhaustive() -> None:
    assert DEFAULT_PARAMETER_PRESET == "exhaustive"
    schema = resolve_parameter_schema(
        mode="basic",
        preset=DEFAULT_PARAMETER_PRESET,
        exhaustiveness=16,
        num_modes=20,
        seed=1,
        box_scale=1.0,
        box_padding=0.0,
        runtime_by_engine={"vina": {}},
    )
    assert schema.common.exhaustiveness == 32
    assert schema.common.exhaustiveness_overridden is False


def test_advanced_exhaustiveness_override_is_recorded() -> None:
    schema = resolve_parameter_schema(
        mode="advanced",
        preset="exhaustive",
        exhaustiveness=8,
        num_modes=20,
        seed=1,
        box_scale=1.0,
        box_padding=0.0,
        runtime_by_engine={"vina": {}},
    )
    assert schema.common.exhaustiveness == 8
    assert schema.common.exhaustiveness_overridden is True


def test_energy_range_is_passed_on_the_vina_command_line(tmp_path: Path) -> None:
    row = PairlistRow(
        receptor=RECEPTOR_NAME,
        site_id="site_1",
        ligand=LIGAND_NAME,
        center_x=0.0,
        center_y=0.0,
        center_z=0.0,
        size_x=4.0,
        size_y=4.0,
        size_z=4.0,
    )
    runtime = {"exhaustiveness": 32, "num_modes": 40, "seed": 5, "energy_range": 3.0}
    for runner_cls in (VinaRunner, SminaRunner):
        runner = runner_cls(tmp_path / runner_cls.name, runtime=dict(runtime))
        command = runner.plan_jobs([row])[0].command
        index = command.index("--energy_range")
        assert command[index + 1] == "3"
    default_runner = VinaRunner(tmp_path / "default", runtime={"seed": 5})
    command = default_runner.plan_jobs([row])[0].command
    assert command[command.index("--energy_range") + 1] == "3"


def test_energy_range_recorded_in_schema_and_runtime() -> None:
    schema = resolve_parameter_schema(
        mode="basic",
        preset="exhaustive",
        exhaustiveness=32,
        num_modes=20,
        seed=3,
        box_scale=1.0,
        box_padding=0.0,
        runtime_by_engine={"vina": {}},
    )
    runtime = apply_schema_to_runtime(schema, {"vina": {}})
    assert runtime["vina"]["energy_range"] == 3.0
    assert runtime["vina"]["replicates"] == 3


def test_poses_returned_and_warning_when_fewer_than_num_modes(tmp_path: Path) -> None:
    fake = tmp_path / "bin" / "fake_vina"
    fake.parent.mkdir()
    fake.write_text(
        "#!/bin/bash\n"
        "out=\"\"\n"
        "while [ $# -gt 0 ]; do if [ \"$1\" = \"--out\" ]; then out=\"$2\"; fi; shift; done\n"
        "cat > \"$out\" <<'EOF'\n"
        "MODEL 1\nREMARK VINA RESULT:     -7.5      0.000      0.000\nENDMDL\n"
        "MODEL 2\nREMARK VINA RESULT:     -7.1      1.000      2.000\nENDMDL\n"
        "EOF\n",
        encoding="utf-8",
    )
    fake.chmod(0o755)
    row = PairlistRow(
        receptor=RECEPTOR_NAME,
        site_id="site_1",
        ligand=LIGAND_NAME,
        center_x=0.0,
        center_y=0.0,
        center_z=0.0,
        size_x=4.0,
        size_y=4.0,
        size_z=4.0,
    )
    runner = VinaRunner(
        tmp_path / "proj",
        runtime={"binary": str(fake), "exhaustiveness": 32, "num_modes": 5, "seed": 9, "replicates": 1,
                 "energy_range": 3.0},
    )
    payload = runner.run([row], dry_run=False)
    job = payload["jobs"][0]
    assert job["status"] == "completed"
    assert job["poses_returned"] == 2
    assert job["seed"] == 9
    assert any("poses_returned=2 is fewer than num_modes=5" in warning for warning in job["warnings"])
    assert payload["effective"]["energy_range"] == 3.0


def test_seed_is_required_by_schema_validation() -> None:
    schema = resolve_parameter_schema(
        mode="basic",
        preset="exhaustive",
        exhaustiveness=32,
        num_modes=20,
        seed=None,
        box_scale=1.0,
        box_padding=0.0,
        runtime_by_engine={"vina": {}},
    )
    errors, _warnings = validate_parameter_schema(schema, ["vina"])
    assert any("`seed` is required" in error for error in errors)


def test_replicates_without_base_seed_is_refused(tmp_path: Path) -> None:
    runner = VinaRunner(tmp_path / "p", runtime={"replicates": 3})
    with pytest.raises(ValueError, match="base seed"):
        runner.replicate_seeds()


def test_transform_keeps_box_provenance_and_records_effective_edge() -> None:
    row = PairlistRow(
        receptor=RECEPTOR_NAME,
        site_id="site_1",
        ligand=LIGAND_NAME,
        center_x=1.0,
        center_y=2.0,
        center_z=3.0,
        size_x=10.0,
        size_y=10.0,
        size_z=10.0,
        box_method=BOX_METHOD_RG,
        ligand_rg_angstrom=3.0,
        edge_angstrom=10.0,
    )
    schema = resolve_parameter_schema(
        mode="basic",
        preset="exhaustive",
        exhaustiveness=32,
        num_modes=20,
        seed=1,
        box_scale=1.5,
        box_padding=0.0,
        runtime_by_engine={"vina": {}},
    )
    [out] = transform_pairlist_rows([row], schema)
    assert out.box_method == BOX_METHOD_RG
    assert out.ligand_rg_angstrom == 3.0
    assert out.edge_angstrom == pytest.approx(15.0)
    assert out.size_x == pytest.approx(15.0)


# ----------------------------------------------------------------------------- R6 tags


def test_tag_uses_receptor_stem_without_pdbqt() -> None:
    row = PairlistRow(
        receptor="1IEP_A_protein.pdbqt",
        site_id="site_1",
        ligand="STI_A_201.pdbqt",
        center_x=0.0,
        center_y=0.0,
        center_z=0.0,
        size_x=20.0,
        size_y=20.0,
        size_z=20.0,
    )
    assert row.tag == "1IEP_A_protein_site_1_STI_A_201.pdbqt"
    assert ".pdbqt_site" not in row.tag
    assert row.legacy_tag == "1IEP_A_protein.pdbqt_site_1_STI_A_201.pdbqt"
    assert row.receptor == "1IEP_A_protein.pdbqt"  # file lookup key is unchanged


def test_legacy_tag_pose_files_still_resolve_to_their_pair() -> None:
    row = PairlistRow(
        receptor="1IEP_A_protein.pdbqt",
        site_id="site_1",
        ligand="STI_A_201.pdbqt",
        center_x=0.0,
        center_y=0.0,
        center_z=0.0,
        size_x=20.0,
        size_y=20.0,
        size_z=20.0,
    )
    index = build_pair_index([row])
    assert index[row.tag] is row
    assert index[row.legacy_tag] is row
    assert split_replicate_stem(replicate_job_tag(row.tag, 2)) == (row.tag, 2)
    assert split_replicate_stem(row.legacy_tag) == (row.legacy_tag, None)


def test_pairlist_tag_in_pairlist_csv_has_no_pdbqt_receptor_suffix(pair_patches, tmp_path: Path) -> None:
    _summary, pairlist = _build_pairlist(tmp_path, ligand_atoms=_hexagon_atoms())
    receptor = pairlist.iloc[0]["receptor"]
    assert receptor == RECEPTOR_NAME  # the file name is kept for lookup
    row = PairlistRow(**{key: pairlist.iloc[0][key] for key in ["receptor", "site_id", "ligand", "center_x", "center_y",
                                                             "center_z", "size_x", "size_y", "size_z"]})
    assert row.tag == "1IEP_receptor_site_1_1IEP_ligand_STI_A_201.pdbqt"


# ----------------------------------------------------------------------------- R6 status and Excel


def test_informational_only_notes_give_completed_exit_zero() -> None:
    compatibility = validate_ligand_preparation_profile("engine_aware_full", selected_engines=[])
    assert compatibility.errors == []
    assert compatibility.warnings == []
    assert compatibility.info, "the no-engine note must be classified as informational"
    payload = compatibility.to_dict()
    assert payload["status"] == "completed"
    assert payload["exit_code"] == 0


def test_real_warning_gives_completed_with_warnings_exit_zero() -> None:
    compatibility = validate_ligand_preparation_profile("autodocktools_only", selected_engines=["vina"])
    assert compatibility.warnings
    assert compatibility.status == "completed_with_warnings"
    assert compatibility.exit_code == 0


def test_error_gives_failed_exit_one() -> None:
    compatibility = validate_ligand_preparation_profile("openbabel_only", selected_engines=["autodock4"])
    assert compatibility.errors
    assert classify_preparation_status(compatibility.errors, []) == "failed"
    assert compatibility.exit_code == 1


def test_excel_centre_without_method_is_refused(tmp_path: Path) -> None:
    path = _write_excel(tmp_path / "no_method.xlsx", method=None)
    with pytest.raises(ValueError, match="without a recorded method"):
        excel_sites.load_site_catalog_from_summary(path, allow_missing_coordinates=True)


def test_excel_centre_with_unknown_method_is_refused(tmp_path: Path) -> None:
    path = _write_excel(tmp_path / "bad_method.xlsx", method="guessed_centre")
    with pytest.raises(ValueError, match="without a recorded method"):
        excel_sites.load_site_catalog_from_summary(path, allow_missing_coordinates=True)


@pytest.mark.parametrize("method", ["binding_site_center_v1", "user_explicit"])
def test_excel_centre_with_recorded_method_is_accepted(tmp_path: Path, method: str) -> None:
    path = _write_excel(tmp_path / f"{method}.xlsx", method=method)
    catalog = excel_sites.load_site_catalog_from_summary(path, allow_missing_coordinates=True)
    assert catalog.iloc[0]["center_method"] == method


def test_excel_cell_value_serialises_dict_provenance_as_json() -> None:
    from core_pipeline import excel_cell_value

    provenance = {"status": "completed", "center_method": "binding_site_center_v1"}
    text = excel_cell_value(provenance)
    assert isinstance(text, str)
    assert json.loads(text) == provenance
    assert excel_cell_value(1.5) == 1.5


def test_pairlist_refuses_excel_centre_without_method(pair_patches, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="without a recorded method"):
        _build_pairlist(tmp_path, ligand_atoms=_hexagon_atoms(), method=None)
