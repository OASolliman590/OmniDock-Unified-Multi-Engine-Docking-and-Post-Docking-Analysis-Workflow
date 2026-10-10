from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import pytest

from post_docking_analysis.complex_validation import validate_md_receptor_pdb, validate_md_system_pdb
from post_docking_analysis.md_chemistry import (
    MDAtom,
    MDBond,
    PreparedLigand,
    get_openbabel_backend,
    openbabel_capability,
)
from post_docking_analysis.md_inputs import (
    CONSUMER_PROFILE,
    MDInputsRequest,
    build_md_inputs_config_payload,
    export_md_inputs,
    load_exact_tags,
    sha256_file,
)


def _pdb_atom(
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
    *,
    record: str = "ATOM",
    altloc: str = " ",
) -> str:
    chars = [" "] * 80
    chars[0:6] = list(f"{record:<6s}"[:6])
    chars[6:11] = list(f"{serial:5d}")
    chars[12:16] = list(f"{name:>4s}"[-4:])
    chars[16] = altloc
    chars[17:20] = list(f"{resname:>3s}"[-3:])
    chars[21] = chain
    chars[22:26] = list(f"{resseq:4d}")
    chars[30:38] = list(f"{x:8.3f}")
    chars[38:46] = list(f"{y:8.3f}")
    chars[46:54] = list(f"{z:8.3f}")
    chars[54:60] = list(f"{1.0:6.2f}")
    chars[60:66] = list(f"{10.0:6.2f}")
    chars[76:78] = list(f"{element:>2s}"[-2:])
    return "".join(chars)


def _assert_schema_subset(instance: object, schema: dict[str, object], location: str = "$") -> None:
    """Validate the JSON-Schema keywords used by the versioned Spec 032 schemas."""
    expected_type = schema.get("type")
    type_map = {
        "object": dict,
        "array": list,
        "string": str,
        "integer": int,
        "number": (int, float),
    }
    if expected_type:
        assert isinstance(instance, type_map[str(expected_type)]), f"{location}: expected {expected_type}"
        if expected_type in {"integer", "number"}:
            assert not isinstance(instance, bool), f"{location}: booleans are not numbers"
    if "const" in schema:
        assert instance == schema["const"], f"{location}: const mismatch"
    if "enum" in schema:
        assert instance in schema["enum"], f"{location}: enum mismatch"
    if isinstance(instance, str):
        assert len(instance) >= int(schema.get("minLength", 0)), f"{location}: too short"
        if schema.get("pattern"):
            assert re.fullmatch(str(schema["pattern"]), instance), f"{location}: pattern mismatch"
    if isinstance(instance, (int, float)) and not isinstance(instance, bool):
        if "minimum" in schema:
            assert instance >= schema["minimum"], f"{location}: below minimum"
        if "maximum" in schema:
            assert instance <= schema["maximum"], f"{location}: above maximum"
    if isinstance(instance, dict):
        for key in schema.get("required", []):
            assert key in instance, f"{location}: missing required field {key}"
        properties = schema.get("properties", {})
        for key, child_schema in properties.items():
            if key in instance:
                _assert_schema_subset(instance[key], child_schema, f"{location}.{key}")
    if isinstance(instance, list) and schema.get("items"):
        for index, value in enumerate(instance):
            _assert_schema_subset(value, schema["items"], f"{location}[{index}]")


def _write_protein(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                _pdb_atom(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
                _pdb_atom(2, "CA", "ALA", "A", 1, 1.4, 0.0, 0.0, "C"),
                _pdb_atom(3, "C", "ALA", "A", 1, 2.1, 1.2, 0.0, "C"),
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_dna(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        _pdb_atom(1, "P", "DA", "A", 1, 0.0, 0.0, 0.0, "P") + "\nEND\n",
        encoding="utf-8",
    )


def _sdf_record(name: str, affinity: float, x_offset: float) -> str:
    return f"""{name}
  DockForge

  2  1  0  0  0  0  0  0  0  0999 V2000
{x_offset:10.4f}    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
{x_offset + 1.2:10.4f}    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
>  <minimizedAffinity>
{affinity}

$$$$
"""


class FakeChemistryBackend:
    name = "fake_contract_backend"
    version = "1.0"

    def __init__(self) -> None:
        self.source_poses = []

    def prepare(
        self,
        source_file,
        *,
        source_pose,
        pH,
        protonate,
        pose_heavy_coordinates=None,
        pose_heavy_elements=None,
        pose_to_topology=None,
    ):
        self.source_poses.append(source_pose)
        if pose_heavy_coordinates:
            c_xyz, o_xyz = pose_heavy_coordinates
        else:
            # The test's selected SDF record is pose 2 at x=10 and x=11.2.
            c_xyz, o_xyz = (10.0, 0.0, 0.0), (11.2, 0.0, 0.0)
        atoms = [
            MDAtom(1, "C", *c_xyz, name="C1"),
            MDAtom(2, "O", *o_xyz, name="O1"),
            MDAtom(3, "H", c_xyz[0], 1.0, 0.0, name="H1"),
            MDAtom(4, "H", c_xyz[0], -1.0, 0.0, name="H2"),
        ]
        return PreparedLigand(
            atoms=atoms,
            bonds=[MDBond(1, 2, "1"), MDBond(1, 3, "1"), MDBond(1, 4, "1")],
            net_charge=0,
            backend=self.name,
            backend_version=self.version,
            protonated=protonate,
            pH=pH,
            implicit_hydrogen_count=0,
            heavy_coordinate_max_delta=0.0,
        )

    def write(self, ligand, output_dir, formats):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        written = {}
        requested = set(formats) | {"mol2"}
        for fmt in requested:
            path = output_dir / f"ligand.{fmt}"
            if fmt == "mol2":
                atom_lines = [
                    f"{atom.index:7d} {atom.name:<4s} {atom.x:10.4f} {atom.y:10.4f} {atom.z:10.4f} {atom.element} 1 LIG {atom.formal_charge:.4f}"
                    for atom in ligand.atoms
                ]
                bond_lines = [
                    f"{index:6d} {bond.begin:4d} {bond.end:4d} {bond.order}"
                    for index, bond in enumerate(ligand.bonds, start=1)
                ]
                path.write_text(
                    "@<TRIPOS>MOLECULE\nLIG\n"
                    f"{len(ligand.atoms)} {len(ligand.bonds)} 0 0 0\nSMALL\nUSER_CHARGES\n\n"
                    "@<TRIPOS>ATOM\n" + "\n".join(atom_lines) + "\n"
                    "@<TRIPOS>BOND\n" + "\n".join(bond_lines) + "\n",
                    encoding="utf-8",
                )
            else:
                path.write_text(
                    f"fake {fmt}\natoms={len(ligand.atoms)} bonds={len(ligand.bonds)} charge={ligand.net_charge}\n",
                    encoding="utf-8",
                )
            written[fmt] = path
        return written


def _project(tmp_path: Path, *, receptor_writer=_write_protein):
    project = tmp_path / "project"
    selection_dir = project / "4-Working" / "scores" / "unified"
    selection_dir.mkdir(parents=True)
    pose = project / "3-Docking" / "gnina" / "poses" / "prot_site_lig.sdf"
    pose.parent.mkdir(parents=True)
    pose.write_text(_sdf_record("pose1", -7.0, 0.0) + _sdf_record("pose2", -8.5, 10.0), encoding="utf-8")
    pd.DataFrame(
        [
            {
                "engine": "gnina",
                "tag": "prot_site_lig",
                "protein": "prot",
                "ligand": "lig",
                "pose": 2,
                "pose_file": str(pose),
                "affinity_kcal_mol": -8.5,
            }
        ]
    ).to_csv(selection_dir / "best_pose_per_tag_by_engine.csv", index=False)
    pd.DataFrame(
        [
            {
                "receptor": "prot.pdbqt",
                "site_id": "site",
                "ligand": "lig.pdbqt",
                "center_x": 0.0,
                "center_y": 0.0,
                "center_z": 0.0,
                "size_x": 20.0,
                "size_y": 20.0,
                "size_z": 20.0,
            }
        ]
    ).to_csv(project / "pairlist.csv", index=False)
    (project / "project_manifest.json").write_text(
        json.dumps(
            {
                "project_name": "spec032",
                "project_root": str(project),
                "layout_profile": "canonical",
                "engines": ["gnina", "vina", "smina"],
                "pairlist_file": str(project / "pairlist.csv"),
                "favorite_engine": "gnina",
            }
        ),
        encoding="utf-8",
    )
    tags = project / "tags.csv"
    pd.DataFrame([{"tag": "prot_site_lig"}]).to_csv(tags, index=False)
    receptor = project / "prepared" / "prot.pdb"
    receptor_writer(receptor)
    receptor_map = project / "receptor_map.csv"
    pd.DataFrame([{"tag": "prot_site_lig", "receptor_file": str(receptor)}]).to_csv(receptor_map, index=False)
    return project, pose, tags, receptor, receptor_map


def _request(project, tags, receptor_map, **kwargs):
    payload = {
        "project_dir": str(project),
        "engine": "gnina",
        "tags_file": str(tags),
        "receptor_map": str(receptor_map),
        "pH": 7.4,
        "protonation_policy": "openbabel_predicted",
        "consumer_profile": CONSUMER_PROFILE,
        "ligand_formats": ["mol2", "sdf"],
        "protonate": True,
        **kwargs,
    }
    return MDInputsRequest.from_dict(payload)


def test_exact_tags_reject_duplicates(tmp_path):
    tags = tmp_path / "tags.txt"
    tags.write_text("a\na\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate_requested_tag"):
        load_exact_tags(tags)


def test_non_first_gnina_pose_exports_strict_system_and_provenance(tmp_path):
    project, pose, tags, receptor, receptor_map = _project(tmp_path)
    backend = FakeChemistryBackend()
    result = export_md_inputs(_request(project, tags, receptor_map), chemistry_backend=backend)

    assert result["status"] == "completed"
    assert backend.source_poses == [2]
    row = result["rows"][0]
    assert row["status"] == "completed"
    destination = Path(row["output_dir"])
    provenance = json.loads((destination / "provenance.json").read_text(encoding="utf-8"))
    assert provenance["schema_version"] == "md_inputs_provenance_v1"
    assert provenance["pose_no"] == 2
    assert provenance["affinity_kcal_mol"] == pytest.approx(-8.5)
    assert provenance["net_charge"] == 0
    assert provenance["charge_authority"] == "predicted"
    assert provenance["scientific_review"] == "human_review_required"
    assert provenance["inputs"]["pose"]["sha256"] == sha256_file(pose)
    assert provenance["gates"]["G5"]["amino_acid_residue_count"] == 1
    assert provenance["docking_provenance"]["engine"] == "gnina"
    assert provenance["docking_provenance"]["engine_version"] == "unknown"
    assert set(provenance["gates"]) == {"G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8"}
    system = destination / "system.pdb"
    strict = validate_md_system_pdb(
        system,
        expected_receptor_atoms=3,
        expected_ligand_atoms=4,
        ligand_chain="Z",
        expected_ligand_bonds=[(1, 2), (1, 3), (1, 4)],
    )
    assert strict["is_valid"], strict
    ligand_lines = [line for line in system.read_text(encoding="utf-8").splitlines() if line.startswith("HETATM")]
    assert all(len(line) >= 78 and line[16] == " " and line[21] == "Z" for line in ligand_lines)
    assert Path(result["run_tracking_manifest_file"]).is_file()
    assert Path(result["config_file"]) == project / ".meta" / "md_inputs_config.json"


def test_versioned_schemas_and_all_manifest_hashes_validate(tmp_path):
    project, _, tags, _, receptor_map = _project(tmp_path)
    result = export_md_inputs(_request(project, tags, receptor_map), chemistry_backend=FakeChemistryBackend())
    manifest = json.loads(Path(result["manifest_file"]).read_text(encoding="utf-8"))
    provenance = json.loads(Path(result["rows"][0]["provenance_file"]).read_text(encoding="utf-8"))
    schema_root = Path(__file__).parents[1] / "post_docking_analysis" / "schemas"
    manifest_schema = json.loads((schema_root / "md_inputs_manifest_v1.schema.json").read_text(encoding="utf-8"))
    provenance_schema = json.loads((schema_root / "md_inputs_provenance_v1.schema.json").read_text(encoding="utf-8"))
    _assert_schema_subset(manifest, manifest_schema)
    _assert_schema_subset(provenance, provenance_schema)
    for section in ("inputs", "outputs"):
        for artifact in provenance[section].values():
            path = Path(artifact["path"])
            assert path.is_file()
            assert artifact["sha256"] == sha256_file(path)
            assert artifact["size_bytes"] == path.stat().st_size


def test_content_hash_cache_reuse_and_pose_change_invalidation(tmp_path):
    project, pose, tags, _, receptor_map = _project(tmp_path)
    backend = FakeChemistryBackend()
    request = _request(project, tags, receptor_map)
    first = export_md_inputs(request, chemistry_backend=backend)
    second = export_md_inputs(request, chemistry_backend=backend)
    assert first["rows"][0]["cache_reused"] is False
    assert second["rows"][0]["cache_reused"] is True
    assert backend.source_poses == [2]

    pose.write_text(pose.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    third = export_md_inputs(request, chemistry_backend=backend)
    assert third["rows"][0]["status"] == "completed"
    assert third["rows"][0]["cache_reused"] is False
    assert backend.source_poses == [2, 2]


def test_dna_only_receptor_fails_protein_profile_without_outputs(tmp_path):
    project, _, tags, _, receptor_map = _project(tmp_path, receptor_writer=_write_dna)
    backend = FakeChemistryBackend()
    result = export_md_inputs(_request(project, tags, receptor_map), chemistry_backend=backend)
    row = result["rows"][0]
    assert row["status"] == "failed"
    assert "no_recognized_amino_acid" in row["reason"]
    assert backend.source_poses == []
    assert not (project / "5-Analysis" / "md_inputs" / "prot" / "lig" / "gnina").exists()


def test_prepared_mmcif_receptor_is_deterministically_materialized(tmp_path):
    from Bio.PDB import MMCIFIO, PDBParser

    project, _, tags, receptor, receptor_map = _project(tmp_path)
    structure = PDBParser(QUIET=True).get_structure("prot", str(receptor))
    mmcif = receptor.with_suffix(".cif")
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(mmcif))
    pd.DataFrame([{"tag": "prot_site_lig", "receptor_file": str(mmcif)}]).to_csv(receptor_map, index=False)
    result = export_md_inputs(_request(project, tags, receptor_map), chemistry_backend=FakeChemistryBackend())
    assert result["status"] == "completed"
    provenance = json.loads(Path(result["rows"][0]["provenance_file"]).read_text(encoding="utf-8"))
    assert provenance["gates"]["G5"]["conversion"] == "biopython_mmcif_to_pdb"
    materialized = Path(provenance["gates"]["G5"]["materialized_pdb"])
    assert materialized.is_file()
    assert sha256_file(materialized) == provenance["gates"]["G5"]["materialized_pdb_sha256"]


def test_pdbqt_only_vina_is_not_comparable(tmp_path):
    project, _, tags, receptor, receptor_map = _project(tmp_path)
    pose = project / "3-Docking" / "vina" / "poses" / "prot_site_lig.pdbqt"
    pose.parent.mkdir(parents=True)
    atom = _pdb_atom(1, "C1", "LIG", "B", 1, 1.0, 2.0, 3.0, "C", record="HETATM")
    pose.write_text("MODEL 1\nREMARK VINA RESULT: -8.500 0.000 0.000\n" + atom + "  0.000 C\nENDMDL\n", encoding="utf-8")
    selection = project / "4-Working" / "scores" / "unified" / "best_pose_per_tag_by_engine.csv"
    frame = pd.read_csv(selection)
    frame.loc[0, ["engine", "pose", "pose_file"]] = ["vina", 1, str(pose)]
    frame.to_csv(selection, index=False)
    request = MDInputsRequest.from_dict(
        {
            "project_dir": str(project),
            "engine": "vina",
            "tags_file": str(tags),
            "receptor_map": str(receptor_map),
            "pH": 7.4,
            "protonation_policy": "openbabel_predicted",
        }
    )
    result = export_md_inputs(request, chemistry_backend=FakeChemistryBackend())
    assert result["rows"][0]["status"] == "not_comparable"
    assert result["rows"][0]["reason"] == "pdbqt_requires_explicit_topology_map"


def test_vina_with_explicit_topology_and_atom_map_exports_pose_coordinates(tmp_path):
    project, pose, tags, _, receptor_map = _project(tmp_path)
    pose = project / "3-Docking" / "vina" / "poses" / "prot_site_lig.pdbqt"
    pose.parent.mkdir(parents=True)
    carbon = _pdb_atom(1, "C1", "LIG", "B", 1, 20.0, 0.0, 0.0, "C", record="HETATM")
    oxygen = _pdb_atom(2, "O1", "LIG", "B", 1, 21.2, 0.0, 0.0, "O", record="HETATM")
    pose.write_text(
        "MODEL 1\nREMARK VINA RESULT: -8.500 0.000 0.000\n"
        + carbon + "  0.000 C\n" + oxygen + "  0.000 OA\nENDMDL\n",
        encoding="utf-8",
    )
    selection = project / "4-Working" / "scores" / "unified" / "best_pose_per_tag_by_engine.csv"
    frame = pd.read_csv(selection)
    frame.loc[0, ["engine", "pose", "pose_file"]] = ["vina", 1, str(pose)]
    frame.to_csv(selection, index=False)
    topology = project / "topologies" / "lig.sdf"
    topology.parent.mkdir(parents=True)
    topology.write_text(_sdf_record("topology", 0.0, 0.0), encoding="utf-8")
    topology_map = project / "topology_map.csv"
    pd.DataFrame(
        [{"tag": "prot_site_lig", "topology_file": str(topology), "topology_pose": 1, "atom_map": "[0, 1]"}]
    ).to_csv(topology_map, index=False)
    request = MDInputsRequest.from_dict(
        {
            "project_dir": str(project),
            "engine": "vina",
            "tags_file": str(tags),
            "receptor_map": str(receptor_map),
            "topology_map": str(topology_map),
            "pH": 7.4,
            "protonation_policy": "openbabel_predicted",
        }
    )
    result = export_md_inputs(request, chemistry_backend=FakeChemistryBackend())
    assert result["status"] == "completed"
    provenance = json.loads(Path(result["rows"][0]["provenance_file"]).read_text(encoding="utf-8"))
    assert provenance["gates"]["G3"]["atom_map"] == [0, 1]
    ligand_lines = [
        line for line in (Path(result["rows"][0]["output_dir"]) / "system.pdb").read_text(encoding="utf-8").splitlines()
        if line.startswith("HETATM")
    ]
    assert float(ligand_lines[0][30:38]) == pytest.approx(20.0)
    assert float(ligand_lines[1][30:38]) == pytest.approx(21.2)


def test_approved_charge_map_changes_authority_and_rejects_mismatch(tmp_path):
    project, _, tags, _, receptor_map = _project(tmp_path)
    charge_map = project / "charge_map.csv"
    pd.DataFrame([{"tag": "prot_site_lig", "net_charge": 0}]).to_csv(charge_map, index=False)
    success = export_md_inputs(
        _request(project, tags, receptor_map, charge_map=str(charge_map)),
        chemistry_backend=FakeChemistryBackend(),
    )
    provenance = json.loads(Path(success["rows"][0]["provenance_file"]).read_text(encoding="utf-8"))
    assert provenance["charge_authority"] == "human_approved"
    assert provenance["scientific_review"] == "completed"

    pd.DataFrame([{"tag": "prot_site_lig", "net_charge": -1}]).to_csv(charge_map, index=False)
    mismatch = export_md_inputs(
        _request(project, tags, receptor_map, charge_map=str(charge_map), force=True),
        chemistry_backend=FakeChemistryBackend(),
    )
    assert mismatch["rows"][0]["status"] == "failed"
    assert mismatch["rows"][0]["reason"] == "approved_charge_mismatch:0!=-1"


def test_missing_openbabel_is_typed_and_emits_no_row_outputs(tmp_path, monkeypatch):
    project, _, tags, _, receptor_map = _project(tmp_path)
    monkeypatch.setattr("post_docking_analysis.md_inputs.get_openbabel_backend", lambda: None)
    monkeypatch.setattr(
        "post_docking_analysis.md_inputs.openbabel_capability",
        lambda: {"status": "skipped_missing_dependency", "backend": "openbabel", "backend_version": "", "reason": "test_missing"},
    )
    result = export_md_inputs(_request(project, tags, receptor_map))
    assert result["rows"][0]["status"] == "skipped_missing_dependency"
    assert result["rows"][0]["gates"]["G6"]["reason"] == "test_missing"
    assert not (project / "5-Analysis" / "md_inputs" / "prot" / "lig" / "gnina").exists()


def test_shifted_or_blank_chain_receptor_fails_strict_columns(tmp_path):
    receptor = tmp_path / "shifted.pdb"
    valid = _pdb_atom(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C")
    shifted = valid[:21] + " " + valid[22:]
    receptor.write_text(shifted + "\nEND\n", encoding="utf-8")
    result = validate_md_receptor_pdb(receptor)
    assert result["is_valid"] is False
    assert any("missing_chain_id_column_22" in error for error in result["errors"])


def test_config_payload_hashes_all_resolved_inputs(tmp_path):
    project, pose, tags, receptor, receptor_map = _project(tmp_path)
    payload = build_md_inputs_config_payload(_request(project, tags, receptor_map))
    by_path = {row["path"]: row for row in payload["dependencies"]}
    assert by_path[str(pose)]["sha256"] == sha256_file(pose)
    assert by_path[str(receptor)]["sha256"] == sha256_file(receptor)
    assert len(payload["request_fingerprint"]) == 64


def test_md_inputs_dag_node_is_explicit_optional_scope(tmp_path):
    from post_docking_analysis.multi_engine_pipeline import MultiEngineAnalysisPipeline

    project, _, tags, _, receptor_map = _project(tmp_path)
    request = _request(project, tags, receptor_map)
    pipeline = MultiEngineAnalysisPipeline(
        project_dir=str(project),
        output_dir=str(project / "5-Analysis" / "sessions" / "test_md"),
        analysis_mode="single_engine",
        engine="gnina",
        favorite_engine="gnina",
        engines_in_scope=["gnina"],
        analysis_scope="md_inputs",
        md_inputs_request=request.to_dict(),
    )
    graph = pipeline.build_artifact_graph()
    assert pipeline.resolve_dag_scope_artifact("md_inputs") == "md_inputs"
    assert graph.nodes["md_inputs"].optional is True
    assert str(project / "4-Working" / "scores" / "unified" / "best_pose_per_tag_by_engine.csv") in graph.nodes["md_inputs"].inputs


def test_workflow_dispatch_executes_md_inputs_through_dag(tmp_path, monkeypatch):
    import post_docking_analysis.multi_engine_pipeline_impl as pipeline_module
    from workflow.execution import run_analysis_md_inputs

    project, pose, tags, _, receptor_map = _project(tmp_path)
    backend = FakeChemistryBackend()
    real_export = export_md_inputs
    monkeypatch.setattr(
        pipeline_module,
        "export_md_inputs",
        lambda request: real_export(request, chemistry_backend=backend),
    )
    result = run_analysis_md_inputs(
        str(project),
        engine="gnina",
        tags_file=str(tags),
        receptor_map=str(receptor_map),
        ph=7.4,
        protonation_policy="openbabel_predicted",
        ligand_formats=["mol2"],
    )
    assert result.status == "completed"
    assert result.outputs["counts"]["completed"] == 1
    node = result.outputs["dag_report"]["nodes"]["md_inputs"]
    assert node["status"] == "completed"
    repeated = run_analysis_md_inputs(
        str(project),
        engine="gnina",
        tags_file=str(tags),
        receptor_map=str(receptor_map),
        ph=7.4,
        protonation_policy="openbabel_predicted",
        ligand_formats=["mol2"],
    )
    assert repeated.outputs["dag_report"]["nodes"]["md_inputs"]["status"] == "cache_hit"
    assert backend.source_poses == [2]

    pose.write_text(pose.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    invalidated = run_analysis_md_inputs(
        str(project),
        engine="gnina",
        tags_file=str(tags),
        receptor_map=str(receptor_map),
        ph=7.4,
        protonation_policy="openbabel_predicted",
        ligand_formats=["mol2"],
    )
    assert invalidated.outputs["dag_report"]["nodes"]["md_inputs"]["status"] == "completed"
    assert backend.source_poses == [2, 2]


def test_cli_default_is_docked_state_and_a_bare_ph_is_refused(tmp_path):
    # Spec 036 R1b: the default export uses the docked microspecies, so no pH is required. A bare --ph is
    # ambiguous and is refused; re-protonation needs --reprotonate-at-ph (or the legacy flag with --ph).
    from workflow.cli import build_parser, main

    base = ["analyze", "md-inputs", "--project-dir", str(tmp_path), "--engine", "gnina", "--tags-file", "tags.csv", "--receptor-map", "receptors.csv"]
    args = build_parser().parse_args(base)
    assert args.ph is None and args.protonation_policy is None
    with pytest.raises(SystemExit):
        main(base + ["--ph", "7.4"])


@pytest.mark.skipif(
    openbabel_capability()["status"] != "completed",
    reason="skipped_missing_dependency:openbabel",
)
def test_openbabel_backend_protonates_and_roundtrips_mol2(tmp_path):
    source = tmp_path / "ligand.sdf"
    source.write_text(_sdf_record("ligand", 0.0, 0.0), encoding="utf-8")
    backend = get_openbabel_backend()
    assert backend is not None
    ligand = backend.prepare(source, source_pose=1, pH=7.4, protonate=True)
    assert ligand.implicit_hydrogen_count == 0
    assert ligand.heavy_coordinate_max_delta <= 1.0e-4
    assert ligand.net_charge == 0
    assert sum(atom.element.upper() == "H" for atom in ligand.atoms) == 4
    written = backend.write(ligand, tmp_path / "out", ["mol2"])
    assert written["mol2"].is_file()

    project, _, tags, _, receptor_map = _project(tmp_path / "end_to_end")
    result = export_md_inputs(_request(project, tags, receptor_map))
    assert result["status"] == "completed", result
    provenance = json.loads(Path(result["rows"][0]["provenance_file"]).read_text(encoding="utf-8"))
    assert provenance["chemistry_backend"] == {"name": "openbabel", "version": "3.2.1"}
    assert provenance["net_charge"] == 0
    assert provenance["charge_authority"] == "predicted"
    assert provenance["scientific_review"] == "human_review_required"
    assert provenance["gates"]["G6"]["implicit_hydrogen_count"] == 0
    assert provenance["gates"]["G6"]["heavy_coordinate_max_delta"] <= 1.0e-4
