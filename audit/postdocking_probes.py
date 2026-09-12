"""Small adversarial scientific audit probes; run from repository root.

Requires only installed numpy/pandas. Does not run docking or modify pipeline code.
The tiny point sets isolate algorithmic invariants; they are not real ligands.
"""
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from post_docking_analysis import geometric_consensus as geo
from post_docking_analysis import redocking_validation as redock
from post_docking_analysis.consensus import normalize_engine_scores, build_consensus_rankings
from post_docking_analysis.top_pose_selector import build_top_pose_atlas
from post_docking_analysis.biology_integration import attach_biology_annotations

OUT = ROOT / "audit" / "postdocking-probe-artifacts"
OUT.mkdir(parents=True, exist_ok=True)
results = {}

def pdb_atoms(coords, record="HETATM"):
    return "".join(f"{record:<6}{i:5d}  C{i:<2} LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n" for i, (x,y,z) in enumerate(coords,1))

def sdf_block(coords, title="probe"):
    return (f"{title}\n  audit\n\n{len(coords):3d}  0  0  0  0  0            999 V2000\n" +
            "".join(f"{x:10.4f}{y:10.4f}{z:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0\n" for x,y,z in coords) + "M  END\n$$$$\n")

coords = np.array([[0.,0.,0.],[2.,0.,0.],[0.,3.,0.],[0.,0.,4.]])
translated = coords + [100.,0.,0.]
reference = OUT / "reference.pdb"
docked = OUT / "translated.pdb"
reference.write_text(pdb_atoms(coords) + "END\n")
docked.write_text(pdb_atoms(translated) + "END\n")
results["translation"] = {
    "fixed_frame_rmsd": float(np.sqrt(np.mean(np.sum((coords-translated)**2,axis=1)))),
    "redocking_rmsd": redock._compute_pose_rmsd(docked,reference),
    "geometric_rmsd": geo._kabsch_rmsd(coords,translated),
}
rotation = np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
results["rotation_alignment_implementation"] = {"expected_shape_rmsd":0., "actual":geo._kabsch_rmsd(coords,coords@rotation)}

multi_pdb = OUT / "two_models.pdbqt"
multi_pdb.write_text("MODEL        1\n"+pdb_atoms(coords)+"ENDMDL\nMODEL        2\n"+pdb_atoms(translated)+"ENDMDL\n")
results["redocking_pdbqt_model_selection"] = {"requested_model":2,"atom_count":int(redock._parse_pdb_like_heavy_atoms(multi_pdb,pose_index=2)[0].shape[0]),"expected_atom_count":4}
multi_sdf = OUT / "two_records.sdf"
multi_sdf.write_text(sdf_block(coords)+sdf_block(translated))
sdf_coords, _, sdf_error = geo._extract_sdf_pose(multi_sdf,2)
results["sdf_second_record"] = {"atom_count":len(sdf_coords),"error":sdf_error,"expected_atom_count":4}

pose_rows = pd.DataFrame([{"engine":"vina","tag":"T","pose_file":str(reference),"pose":1},{"engine":"smina","tag":"T","pose_file":str(docked),"pose":1},{"engine":"gnina","tag":"T","pose_file":str(OUT/"absent.sdf"),"pose":1}])
results["missing_engine_geometry"] = geo.compute_geometric_consensus(pose_rows).to_dict(orient="records")

score_rows = pd.DataFrame([{"engine":"vina","protein":"P","affinity_kcal_mol":x} for x in [-10.,-7.,-4.]])
results["normalization_direction"] = {m:normalize_engine_scores(score_rows,method=m)["normalized_affinity_score"].tolist() for m in ["per_engine_rank","per_engine_minmax","per_engine_zscore"]}
results["reference_boolean_missingness"] = {repr(v):redock._is_reference_candidate_row(pd.Series({"is_cocrystal_benchmark":v,"pair_source":"user","site_id":"site1","tag":"P_L","ligand":"L"})) for v in [np.nan,"False",False]}

refs = pd.DataFrame([{"protein":"P","ligand":"L0","tag":"reference0","engine":"vina","pose":1,"pose_file":str(docked),"reference_pose_file":str(reference),"is_cocrystal_benchmark":True,"affinity_kcal_mol":-8.}] + [{"protein":"P","ligand":f"L{i}","tag":f"reference{i}","engine":"vina","pose":1,"pose_file":str(OUT/f"missing{i}.pdb"),"is_cocrystal_benchmark":True,"affinity_kcal_mol":-8.} for i in range(1,100)])
refs["reference_pose_file"] = str(reference)
results["validation_coverage"] = redock.run_redocking_validation(project_dir=OUT,best_by_engine=refs,output_dir=OUT/"validation")["summary"]

best = pd.DataFrame([{"tag":"high","protein":"P","ligand":"L","site_id":"a","engine":"vina","pose":1,"pose_file":str(reference),"affinity_kcal_mol":-8.},{"tag":"low","protein":"P","ligand":"L","site_id":"b","engine":"vina","pose":1,"pose_file":str(reference),"affinity_kcal_mol":-7.}])
cons = pd.DataFrame([{"tag":"high","consensus_score":0.9,"agreement_count":1},{"tag":"low","consensus_score":0.1,"agreement_count":1}])
results["top_pose_atlas_direction"] = {mode:build_top_pose_atlas(best,cons,selection_policy="best_consensus",consensus_mode=mode)["top_pose_per_ligand_per_protein"]["tag"].tolist() for mode in ["weighted_hybrid","strict_consensus","dockbox_geometric"]}

# Execute unmodified AST method bodies without importing optional visualization/tool stacks.
pipeline_ast = ast.parse((ROOT/"post_docking_analysis/multi_engine_pipeline_impl.py").read_text(encoding="utf-8"))
pipeline_class = next(n for n in pipeline_ast.body if isinstance(n,ast.ClassDef) and n.name=="MultiEngineAnalysisPipeline")
method_names = {"_select_best_pose_rows","_resolve_best_pose_metric_for_frame"}
probe_class = ast.ClassDef(name="ProbePipeline",bases=[],keywords=[],body=[n for n in pipeline_class.body if isinstance(n,ast.FunctionDef) and n.name in method_names],decorator_list=[])
module = ast.fix_missing_locations(ast.Module(body=[ast.ImportFrom(module="__future__",names=[ast.alias(name="annotations")],level=0),probe_class],type_ignores=[]))
ns={"pd":pd,"np":np,"normalize_best_pose_selection_metric":lambda x:x}
exec(compile(module,"<unmodified pipeline method AST>","exec"),ns)
instance=ns["ProbePipeline"]()
gnina = pd.DataFrame([{"tag":"T","pose":1,"engine":"gnina","affinity_kcal_mol":-8.,"cnn_affinity":8.},{"tag":"T","pose":2,"engine":"gnina","affinity_kcal_mol":-8.,"cnn_affinity":4.}])
selected,_=instance._select_best_pose_rows(gnina,requested_metric="cnn_affinity")
results["gnina_direction"]={"selected_pose":int(selected.iloc[0]["pose"]),"selected_cnn_affinity":float(selected.iloc[0]["cnn_affinity"]),"expected_pose_for_highest_cnn_affinity":1}

annotated, report=attach_biology_annotations(pd.DataFrame([{"protein":"P","ligand":"L","tag":"T"}]),pd.DataFrame([{"protein":"P","ligand":"L","IC50":10.},{"protein":"P","ligand":"L","IC50":20.}]))
results["biology_duplicate_keys"]={"input_docking_rows":1,"output_rows":len(annotated),"reported_match_rate":report["match_rate"]}

(OUT/"results.json").write_text(json.dumps(results,indent=2,default=str),encoding="utf-8")
print(json.dumps(results,indent=2,default=str))
