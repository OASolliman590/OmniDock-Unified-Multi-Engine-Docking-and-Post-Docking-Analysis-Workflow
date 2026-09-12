"""Scientific invariants that must hold independently of implementation details."""
from pathlib import Path
import json
import os
import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from post_docking_analysis.pose_geometry import load_pose_molecule, pose_rmsd, selected_record, PoseGeometryError
from post_docking_analysis.geometric_consensus import compute_geometric_consensus
from post_docking_analysis.redocking_validation import run_redocking_validation, _is_reference_candidate_row, _find_reference_pose_file, _parse_pdb_like_heavy_atoms
from post_docking_analysis.consensus import normalize_engine_scores, build_consensus_rankings, classify_hits_target_aware
from post_docking_analysis.top_pose_selector import build_top_pose_atlas
from post_docking_analysis.pose_extractor import _is_better_pose
from post_docking_analysis.biology_integration import attach_biology_annotations
from post_docking_analysis.docking_parser import parse_vina_pdbqt
from post_docking_analysis.generate_scores_csv import parse_gnina_log, _find_log_files
from post_docking_analysis.multi_engine_pipeline_impl import MultiEngineAnalysisPipeline

def molecule(smiles="CC(=O)[O-]"):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=23) == 0
    return Chem.RemoveHs(mol)

def write_sdf(path, molecules):
    with Chem.SDWriter(str(path)) as writer:
        for mol in molecules:
            writer.write(mol)
    return path

def moved(mol, offset):
    result = Chem.Mol(mol)
    coordinates = mol.GetConformer().GetPositions() + np.asarray(offset)
    for i, xyz in enumerate(coordinates):
        result.GetConformer().SetAtomPosition(i, xyz)
    return result

def row(path, *, engine="vina", tag="P_site_L", pose=1, **kwargs):
    return {"engine":engine,"tag":tag,"protein":"P","ligand":"L","site_id":"site","pose_file":str(path),"pose":pose,"affinity_kcal_mol":-8.,"receptor_frame_id":"experimental-frame", **kwargs}

def test_docking_rmsd_preserves_translation_and_rotation(tmp_path):
    mol = molecule()
    a = write_sdf(tmp_path/"a.sdf",[mol])
    b = write_sdf(tmp_path/"b.sdf",[moved(mol,[100,0,0])])
    assert pose_rmsd(a,b) == pytest.approx(100,abs=1e-3)
    rotated = Chem.Mol(mol)
    rotation = np.array([[0,-1,0],[1,0,0],[0,0,1]])
    for i,xyz in enumerate(mol.GetConformer().GetPositions() @ rotation):
        rotated.GetConformer().SetAtomPosition(i,xyz)
    c=write_sdf(tmp_path/"c.sdf",[rotated])
    assert pose_rmsd(a,c) > .5

def test_graph_permutation_preserves_pose_rmsd(tmp_path):
    mol=molecule()
    permuted=Chem.RenumberAtoms(mol,list(reversed(range(mol.GetNumAtoms()))))
    a=write_sdf(tmp_path/"a.sdf",[mol])
    b=write_sdf(tmp_path/"b.sdf",[permuted])
    assert pose_rmsd(a,b) < 1e-6

def test_same_size_different_chemistry_is_rejected(tmp_path):
    a=write_sdf(tmp_path/"a.sdf",[molecule("CCCO")])
    b=write_sdf(tmp_path/"b.sdf",[molecule("CCOC")])
    with pytest.raises(PoseGeometryError,match="chemical_identity_mismatch"):
        pose_rmsd(a,b)

def test_sdf_selects_exact_record_with_empty_title(tmp_path):
    mol=molecule()
    path=write_sdf(tmp_path/"poses.sdf",[mol,moved(mol,[7,0,0])])
    assert pose_rmsd(path,path,pose_a=1,pose_b=2) == pytest.approx(7,abs=.001)
    with pytest.raises(PoseGeometryError,match="out_of_range"):
        load_pose_molecule(path,3)

def test_pdbqt_selects_single_model_and_never_guesses_graph(tmp_path):
    atom="HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\n"
    path=tmp_path/"poses.pdbqt"
    path.write_text("MODEL 1\n"+atom+"ENDMDL\nMODEL 2\n"+atom+"ENDMDL\n")
    coords,_,error=_parse_pdb_like_heavy_atoms(path,pose_index=2)
    assert coords.shape == (1,3) and not error
    with pytest.raises(PoseGeometryError,match="missing_authoritative"):
        load_pose_molecule(path,1)

def test_meeko_smiles_mapping_preserves_atom_identity(tmp_path):
    path=tmp_path/"pose.pdbqt"
    path.write_text("REMARK SMILES CO\nREMARK SMILES IDX 1 2 2 1\nMODEL 1\n"
        "HETATM    1  O1  LIG A   1       1.400   0.000   0.000  1.00  0.00     0.000 OA\n"
        "HETATM    2  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\nENDMDL\n")
    mol=load_pose_molecule(path)
    assert Chem.MolToSmiles(mol) == "CO"
    assert mol.GetConformer().GetAtomPosition(0).x == 0.
    assert mol.GetConformer().GetAtomPosition(1).x == 1.4


def test_meeko_macrocycle_pseudoatoms_do_not_change_chemical_atom_count(tmp_path):
    path = tmp_path / "closure.pdbqt"
    path.write_text("REMARK SMILES CC\nREMARK SMILES IDX 1 1 2 2\nMODEL 1\n"
                    "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 CG0\n"
                    "HETATM    2  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00     0.000 CG0\n"
                    "HETATM    3  G0  LIG A   1       1.500   0.000   0.000  1.00  0.00     0.000 G0\nENDMDL\n")
    mol = load_pose_molecule(path)
    assert mol.GetNumAtoms() == 2 and Chem.MolToSmiles(mol) == "CC"

@pytest.mark.parametrize("missing",["file","frame","engine"])
def test_geometric_missingness_cannot_pass(tmp_path,missing):
    path=write_sdf(tmp_path/"pose.sdf",[molecule()])
    rows=[row(path,engine="vina"),row(path,engine="smina"),row(path,engine="gnina")]
    if missing == "file": rows[-1]["pose_file"]=str(tmp_path/"absent.sdf")
    if missing == "frame": rows[-1]["receptor_frame_id"]=None
    if missing == "engine": rows.pop()
    result=compute_geometric_consensus(pd.DataFrame(rows),expected_engines=["vina","smina","gnina"]).iloc[0]
    assert not result.geometric_agreement
    assert result.geometric_expected_pairs == 3
    assert result.geometric_valid_pairs < 3

def test_geometric_identical_graphs_common_frame_pass(tmp_path):
    path=write_sdf(tmp_path/"pose.sdf",[molecule()])
    result=compute_geometric_consensus(pd.DataFrame([row(path),row(path,engine="smina")])).iloc[0]
    assert result.geometric_agreement

@pytest.mark.parametrize("flag",[np.nan,None,"False","false",0])
def test_missing_or_false_reference_flags_are_not_evidence(flag):
    assert not _is_reference_candidate_row(pd.Series({"is_cocrystal_benchmark":flag,"ligand":"reference_named_compound","site_id":"reference"}))

def test_reference_discovery_never_uses_raw_conformer(tmp_path):
    raw=tmp_path/"0-Input"/"ligands"
    raw.mkdir(parents=True)
    (raw/"L.pdb").write_text("HETATM\n")
    assert _find_reference_pose_file(tmp_path,pd.Series({"ligand":"L"})) is None

def test_redocking_uses_selected_pose_and_requires_provenance(tmp_path):
    mol=molecule()
    reference=write_sdf(tmp_path/"reference.sdf",[mol])
    poses=write_sdf(tmp_path/"poses.sdf",[mol,moved(mol,[100,0,0])])
    record=row(poses,pose=2,is_cocrystal_benchmark=True,reference_pose_file=str(reference),reference_source="cocrystal",reference_pdb_id="1ABC",reference_frame_id="experimental-frame")
    result=run_redocking_validation(project_dir=tmp_path,best_by_engine=pd.DataFrame([record]),output_dir=tmp_path/"valid")
    assert result["validation_df"].iloc[0].redocking_classification == "fail"
    assert not result["validation_gate"]["allow_reference_anchor"]
    record["pose"]=1
    record["reference_source"]=None
    result=run_redocking_validation(project_dir=tmp_path,best_by_engine=pd.DataFrame([record]),output_dir=tmp_path/"missing")
    assert result["validation_df"].iloc[0].redocking_classification == "not_evaluable"

def test_sparse_validation_does_not_claim_high_confidence(tmp_path):
    path=write_sdf(tmp_path/"reference.sdf",[molecule()])
    records=[row(path,tag=f"P_ref_L{i}",ligand=f"L{i}",is_cocrystal_benchmark=True,reference_pose_file=str(path),reference_source="cocrystal",reference_pdb_id="1ABC",reference_frame_id="experimental-frame") for i in range(100)]
    for record in records[1:]: record["pose_file"]=str(tmp_path/"missing.sdf")
    result=run_redocking_validation(project_dir=tmp_path,best_by_engine=pd.DataFrame(records),output_dir=tmp_path/"validation")
    gate=result["validation_gate"]
    assert gate["evaluable_rows"] == 1 and gate["coverage"] == .01
    assert not gate["allow_reference_anchor"] and gate["confidence_signal"] != "HIGH"


def test_failed_reference_jobs_remain_in_validation_coverage(tmp_path):
    path = write_sdf(tmp_path / "reference.sdf", [molecule()])
    present = row(path, tag="P_site_L1", ligand="L1", is_cocrystal_benchmark=True,
                  reference_pose_file=str(path), reference_source="cocrystal",
                  reference_pdb_id="1ABC", reference_frame_id="experimental-frame")
    expected = pd.DataFrame([{**present, "tag":f"P_site_L{i}", "ligand":f"L{i}"} for i in range(1,4)])
    result = run_redocking_validation(project_dir=tmp_path, best_by_engine=pd.DataFrame([present]),
                                      expected_reference_rows=expected, expected_engines=["vina", "gnina"], output_dir=tmp_path / "validation")
    assert result["summary"]["total_reference_rows"] == 6
    assert result["summary"]["evaluable_rows"] == 1
    assert not result["validation_gate"]["allow_reference_anchor"]


def test_explicit_import_is_hash_checked_and_marked_unverified(tmp_path):
    from post_docking_analysis.score_import import read_explicit_score_import
    from post_docking_analysis.pose_geometry import content_hash
    path = tmp_path / "normalized_scores.csv"
    pd.DataFrame([row("absent")]).to_csv(path, index=False)
    assert read_explicit_score_import(tmp_path) is None
    (tmp_path / "normalized_scores.import.json").write_text(json.dumps({"kind":"imported_scores", "sha256":content_hash(path)}))
    assert read_explicit_score_import(tmp_path).iloc[0].score_provenance == "explicit_import_unverified_execution"
    path.write_text(path.read_text().replace("-8.0", "-9.0"))
    with pytest.raises(ValueError, match="matching CSV sha256"):
        read_explicit_score_import(tmp_path)


@pytest.mark.parametrize("field,value", [("pose",1.5), ("pose",np.nan), ("affinity_kcal_mol",np.inf), ("ligand","")])
def test_explicit_import_rejects_invalid_identity_score_and_pose(tmp_path, field, value):
    from post_docking_analysis.score_import import read_explicit_score_import
    from post_docking_analysis.pose_geometry import content_hash
    record = row("absent")
    record[field] = value
    path = tmp_path / "normalized_scores.csv"
    pd.DataFrame([record]).to_csv(path,index=False)
    (tmp_path / "normalized_scores.import.json").write_text(json.dumps({"kind":"imported_scores", "sha256":content_hash(path)}))
    with pytest.raises(ValueError):
        read_explicit_score_import(tmp_path)

@pytest.mark.parametrize("method",["per_engine_rank","per_engine_minmax","per_engine_zscore"])
@pytest.mark.parametrize("metric,values",[("affinity_kcal_mol",[-10.,-7.,-4.]),("cnn_affinity",[9.,6.,3.])])
def test_all_normalizers_preserve_scientific_direction(method,metric,values):
    frame=pd.DataFrame({"engine":["gnina"]*3,"protein":["P"]*3,metric:values})
    output=normalize_engine_scores(frame,method=method,score_column=metric)
    assert output.normalized_affinity_score.tolist() == pytest.approx([1.,.5,0.])

def test_gnina_selectors_prefer_higher_predicted_pk():
    data=pd.DataFrame([row("unused",pose=1,cnn_affinity=8.),row("unused",pose=2,cnn_affinity=4.)])
    pipeline=MultiEngineAnalysisPipeline.__new__(MultiEngineAnalysisPipeline)
    selected,_=pipeline._select_best_pose_rows(data,requested_metric="cnn_affinity")
    assert selected.iloc[0].pose == 1
    assert _is_better_pose({"cnn_affinity":8,"vina_affinity":-8,"mode":1},{"cnn_affinity":4,"vina_affinity":-8,"mode":2},criterion="cnn_affinity")

@pytest.mark.parametrize("mode",["weighted_hybrid","strict_consensus","dockbox_geometric","favorite_guardrails"])
def test_atlas_agrees_with_consensus_score_direction(mode):
    data=pd.DataFrame([row("unused",tag="good",site_id="site1"),row("unused",tag="poor",site_id="site2")])
    scores=pd.DataFrame([{"tag":"good","consensus_score":.9,"agreement_count":1},{"tag":"poor","consensus_score":.1,"agreement_count":1}])
    result=build_top_pose_atlas(data,scores,selection_policy="best_consensus",consensus_mode=mode)
    assert result["top_pose_per_ligand_per_protein"].iloc[0].tag == "good"

def test_engine_score_offset_does_not_change_consensus_winner_or_rank():
    data=pd.DataFrame([row("unused",engine=engine,tag=tag,ligand=tag,affinity_kcal_mol=score) for engine in ["vina","autodock4"] for tag,score in [("good",-9.),("poor",-4.)]])
    before=build_consensus_rankings(data,consensus_mode="weighted_hybrid")
    data.loc[data.engine.eq("autodock4"),"affinity_kcal_mol"] -= 100
    after=build_consensus_rankings(data,consensus_mode="weighted_hybrid")
    assert before.tag.tolist() == after.tag.tolist()
    assert before.winner_engine.tolist() == after.winner_engine.tolist()
    assert before.consensus_score.tolist() == after.consensus_score.tolist()

def test_reference_anchor_requires_same_engine_and_scoring_function():
    frame=pd.DataFrame([{"protein":"P","tag":"T","consensus_score":.9,"best_affinity_kcal_mol":-9.,"agreement_count":1,"single_engine":True,"winner_engine":"vina","winner_scoring_function":"vina"}])
    baseline=pd.DataFrame([{"protein":"P","engine":"autodock4","scoring_function":"ad4","reference_affinity":-8.,"reference_tag":"R","redocking_classification":"pass"}])
    result=classify_hits_target_aware(frame,policy="reference_anchor",reference_baselines=baseline)
    assert not result.iloc[0].reference_anchor_available


@pytest.mark.parametrize("bad_index", [np.nan, "bad", 0, 1.5])
def test_invalid_reference_pose_index_only_disables_that_row(tmp_path, bad_index):
    path = write_sdf(tmp_path / "reference.sdf", [molecule()])
    records = [row(path, tag=f"T{i}", pose=index, is_cocrystal_benchmark=True,
                   reference_pose_file=str(path), reference_source="experimental:1ABC",
                   reference_pdb_id="1ABC", reference_frame_id="experimental-frame")
               for i, index in enumerate([bad_index, 1])]
    result = run_redocking_validation(project_dir=tmp_path, best_by_engine=pd.DataFrame(records), output_dir=tmp_path / "validation")
    assert result["validation_df"].redocking_classification.tolist() == ["not_evaluable", "pass"]


def test_default_dag_carries_real_molecular_reference_to_gate_ranking_report(tmp_path):
    from docking.project_layout import bootstrap_project_layout, ensure_engine_layout, pairlist_path, shared_receptors_dir
    from post_docking_analysis.unified_pipeline import UnifiedPostDockingPipeline
    bootstrap_project_layout(tmp_path, ["gnina"], layout_profile="docking_legacy")
    layout = ensure_engine_layout(tmp_path, "gnina")
    records = []
    for number in range(4):
        ligand = f"L{number}"
        mol = molecule()
        mol.SetDoubleProp("minimizedAffinity", -8.0 - number)
        mol.SetDoubleProp("CNNaffinity", 6.0 + number)
        mol.SetDoubleProp("CNNscore", .9 - number * .1)
        write_sdf(layout["poses"] / f"P_site_{ligand}.sdf", [mol])
        reference = write_sdf(tmp_path / f"reference{number}.sdf", [mol])
        records.append({"receptor":"P", "site_id":"site", "ligand":ligand,
                        "center_x":0., "center_y":0., "center_z":0., "size_x":20., "size_y":20., "size_z":20.,
                        "is_cocrystal_benchmark":number < 3, "reference_pose_file":str(reference),
                        "reference_source":f"experimental:{number+1}ABC", "reference_pdb_id":f"{number+1}ABC",
                        "reference_frame_id":"source-frame", "receptor_frame_id":"source-frame"})
    pd.DataFrame(records).to_csv(pairlist_path(tmp_path), index=False)
    (shared_receptors_dir(tmp_path) / "P.pdbqt").write_text("ATOM      1  CA  GLY X   1      20.000  20.000  20.000  1.00  0.00     0.000 C\nEND\n")
    pipeline = UnifiedPostDockingPipeline(project_dir=str(tmp_path), output_dir=str(tmp_path / "analysis"),
                                          engines="gnina", analysis_scope="report_only", dag_scope="report_only", hit_class_policy="reference_anchor")
    # Optional rendering availability must not govern scientific table correctness.
    def unavailable_visualizations(paths):
        raise RuntimeError("Optional renderer intentionally unavailable in integration test")
    pipeline._dag_compute_visualizations_node = unavailable_visualizations
    assert pipeline.run(), pipeline.dag_execution_report
    paths = pipeline._dag_artifact_paths()
    raw = pd.read_csv(paths["normalized_scores"])
    assert raw.reference_frame_id.eq("source-frame").all()
    gate = json.loads(paths["validation_gate"].read_text())
    assert gate["evaluable_rows"] == 3 and gate["allow_reference_anchor"]
    ranked = pd.read_csv(paths["consensus_ranked"])
    assert ranked.iloc[0].ligand == "L3"
    classified = pd.read_csv(paths["classified_hits"])
    assert classified.qc_status.eq("not_evaluated_physical_validity").all()
    assert paths["reports"].is_file() and paths["reports"].stat().st_size > 0


def test_atlas_cross_engine_identity_and_rank_ignore_raw_score_offset():
    records = [row("unused", engine=engine, tag=tag, ligand=tag, affinity_kcal_mol=score)
               for engine in ["vina", "autodock4"] for tag, score in [("good", -9.), ("poor", -4.)]]
    before = pd.DataFrame(records)
    after = before.copy()
    after.loc[after.engine.eq("autodock4"), "affinity_kcal_mol"] -= 100.
    def atlas(data):
        ranks = build_consensus_rankings(data, consensus_mode="weighted_hybrid")
        return build_top_pose_atlas(data, ranks)["top_pose_per_ligand_per_protein"]
    first, second = atlas(before), atlas(after)
    assert first.engine.tolist() == second.engine.tolist()
    assert first.tag.tolist() == second.tag.tolist()


def test_requested_missing_engine_remains_in_consensus_denominator():
    result = build_consensus_rankings(pd.DataFrame([row("missing")]), expected_engines=["vina", "gnina", "smina"])
    assert result.iloc[0].agreement_fraction == pytest.approx(1 / 3)
    assert result.iloc[0].geometric_expected_pairs == 3
    assert not result.iloc[0].geometric_agreement


def test_sole_autodock4_route_keeps_and_ranks_its_own_scores():
    pipeline = MultiEngineAnalysisPipeline.__new__(MultiEngineAnalysisPipeline)
    pipeline._load_engine_scope_config = lambda paths: {"engines_in_scope":["autodock4"]}
    strategy = pipeline._resolve_dag_consensus_strategy({})
    assert strategy == "single/autodock4"
    data = pd.DataFrame([row("unused", engine="autodock4", pose=1, affinity_kcal_mol=-9.),
                         row("unused", engine="autodock4", pose=2, affinity_kcal_mol=-4.)])
    best, metadata = pipeline._prepare_dag_consensus_inputs(data, strategy)
    assert len(best) == 1 and best.iloc[0].pose == 1
    assert metadata["primary_score_name"] == "autodock4_affinity"


def test_pdbqt_heavy_atom_denominator_uses_one_pose_and_excludes_pseudoatoms(tmp_path):
    lines = [f"HETATM{i:5d}  X   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 {kind}"
             for i, kind in enumerate(["C", "OA", "CG0", "H", "HD", "HS", "G0", "G3"],1)]
    model = "\n".join(lines)
    path = tmp_path / "poses.pdbqt"
    path.write_text(f"MODEL 1\n{model}\nENDMDL\nMODEL 2\n{model}\nENDMDL\n")
    assert MultiEngineAnalysisPipeline._count_heavy_atoms_from_pdbqt(path) == 3


def test_agreement_table_reports_named_consensus_representative_score(tmp_path):
    data = pd.DataFrame([row("unused", engine="vina", affinity_kcal_mol=-8.),
                         row("unused", engine="autodock4", affinity_kcal_mol=-108.)])
    input_file = tmp_path / "scores.csv"
    data.to_csv(input_file,index=False)
    pipeline = MultiEngineAnalysisPipeline.__new__(MultiEngineAnalysisPipeline)
    pipeline._resolve_dag_consensus_strategy = lambda paths: "multi"
    pipeline._prepare_dag_consensus_inputs = lambda scores,strategy: (scores,{})
    pipeline._run_dag_consensus_strategy = lambda scores,strategy: pd.DataFrame([{"tag":"P_site_L", "winner_engine":"vina", "best_affinity_kcal_mol":-8.}])
    pipeline._finalize_dag_consensus_output = lambda scores,strategy: scores
    paths = {"normalized_scores":input_file, "consensus_ranked":tmp_path / "ranked.csv", "engine_agreement":tmp_path / "agreement.csv"}
    pipeline._dag_compute_consensus_ranked_node(paths)
    result = pd.read_csv(paths["engine_agreement"]).iloc[0]
    assert result.engine_support_count == 2
    assert result.winner_engine == "vina" and result.best_affinity_kcal_mol == -8.

def test_unchecked_geometry_is_not_qc_pass():
    frame=pd.DataFrame([{"protein":"P","tag":"T","consensus_score":.9,"best_affinity_kcal_mol":-9.,"agreement_count":1,"single_engine":True}])
    result=classify_hits_target_aware(frame)
    assert result.iloc[0].qc_status == "not_evaluated_physical_validity"

def test_biology_replicates_cannot_duplicate_candidates():
    docking=pd.DataFrame([{"protein":"P","ligand":"L","tag":"T"}])
    biology=pd.DataFrame([{"protein":"P","ligand":"L","IC50":10},{"protein":"P","ligand":"L","IC50":20}])
    with pytest.raises(ValueError,match="not unique"):
        attach_biology_annotations(docking,biology)

def test_missing_smina_rmsd_is_not_zero(tmp_path):
    path=tmp_path/"smina.pdbqt"
    path.write_text("MODEL 1\nREMARK minimizedAffinity -8.0\nENDMDL\n")
    result=parse_vina_pdbqt(path).iloc[0]
    assert pd.isna(result.rmsd_lb) and pd.isna(result.rmsd_ub)

def test_runner_and_unlisted_logs_do_not_enter_scores(tmp_path):
    content="mode | affinity | intramol | CNNscore | cnn_affinity\n 1 -8.0 -0.1 0.9 7.0\n"
    runner=tmp_path/"T.runner.log"
    runner.write_text(content)
    unknown=tmp_path/"unknown.log"
    unknown.write_text(content)
    assert runner not in _find_log_files(tmp_path)
    assert not parse_gnina_log(runner)[0]
    assert not parse_gnina_log(unknown,{"T":"T"})[0]

def test_source_cache_detects_same_timestamp_edits_and_deletions(tmp_path):
    layout={name:tmp_path/name for name in ("poses","logs","scores")}
    for path in layout.values(): path.mkdir()
    source=layout["poses"]/"T.pdbqt"
    source.write_text("old")
    old_stat=source.stat()
    before=MultiEngineAnalysisPipeline._source_fingerprint(layout,{"T":{"ligand":"L"}})
    source.write_text("new")
    os.utime(source,ns=(old_stat.st_atime_ns,old_stat.st_mtime_ns))
    edited=MultiEngineAnalysisPipeline._source_fingerprint(layout,{"T":{"ligand":"L"}})
    assert edited != before
    source.unlink()
    assert MultiEngineAnalysisPipeline._source_fingerprint(layout,{"T":{"ligand":"L"}}) != edited

def test_receptor_export_preserves_chains_and_cofactor_records(tmp_path):
    path=tmp_path/"receptor.pdb"
    path.write_text("ATOM      1  CA  ALA X   1       0.000   0.000   0.000  1.00  0.00           C\nHETATM    2 ZN    ZN Z 100       3.000   0.000   0.000  1.00  0.00          ZN\nEND\n")
    pipeline=MultiEngineAnalysisPipeline.__new__(MultiEngineAnalysisPipeline)
    lines=pipeline._read_receptor_lines(path)
    assert [line[21] for line in lines] == ["X","Z"]
    assert lines[1].startswith("HETATM")
