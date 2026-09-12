"""Lightweight audit reproductions; no docking or chemistry tools are invoked.

Run from repository root: python audit/preparation-probes.py
AST extraction isolates core functions whose module imports require unavailable Bio.
These probes demonstrate behavior; the synthetic fixtures are not docking benchmarks.
"""
from __future__ import annotations

import ast
import contextlib
import io
import json
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from docking.preparation.ligand_quality import (
    sanitize_prepared_ligand_pdbqt, validate_prepared_ligand_pdbqt,
    sanitize_ligand_pdb_file,
)
from docking.preparation.receptor_quality import validate_receptor_file
from docking.preparation.project_builder import DockingProjectBuilder


def extracted(filename, name, namespace):
    tree = ast.parse((ROOT / filename).read_text(encoding="utf-8-sig"))
    node = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == name)
    module = ast.Module(body=[ast.ImportFrom(module="__future__", names=[ast.alias(name="annotations")], level=0), node], type_ignores=[])
    exec(compile(ast.fix_missing_locations(module), filename, "exec"), namespace)
    return namespace[name]


def atom(serial, x=0.0, y=0.0, z=0.0, kind="C", name="C1", alt="", occ=1.0):
    return f"HETATM{serial:5d} {name:>4s}{alt:1s}LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{0.0:6.2f}          {kind:>2s}  "


def main():
    result = {}
    with tempfile.TemporaryDirectory(prefix="omnidock-preparation-audit-") as temporary:
        folder = Path(temporary)
        macro = folder / "macrocycle.pdbqt"
        macro.write_text("ROOT\n" + atom(1)[:66] + "     0.000 CG0\n" + atom(2)[:66] + "     0.000 G0\nENDROOT\nTORSDOF 1\n")
        before = [i.reason for i in validate_prepared_ligand_pdbqt(macro)]
        repaired = sanitize_prepared_ligand_pdbqt(macro, macro)
        after = [i.reason for i in validate_prepared_ligand_pdbqt(macro)]
        result["macrocycle_tokens"] = {"issues_before": before, "dropped": repaired["dropped_atoms"], "retyped": repaired["replaced_atoms"], "issues_after": after}
        assert before == ["invalid_atom_types"] and repaired["dropped_atoms"] == 1 and after == []

        receptor = folder / "nonfinite.pdb"
        receptor.write_text("\n".join(atom(n, x=float("nan")) for n in range(1, 101)) + "\n")
        result["nonfinite_receptor"] = [i.reason for i in validate_receptor_file(receptor)]
        assert result["nonfinite_receptor"] == []

        ligand = folder / "connectivity.pdb"
        ligand.write_text(atom(1, kind="N", name="N1")[:78] + "1+\n" + atom(2, x=1.3) + "\nCONECT    1    2\nCONECT    2    1\nEND\n")
        sanitize_ligand_pdb_file(ligand, ligand)
        clean = ligand.read_text()
        result["ligand_metadata"] = {"conect_retained": "CONECT" in clean, "formal_charge_retained": "1+" in clean}
        assert result["ligand_metadata"] == {"conect_retained": False, "formal_charge_retained": False}

        # If an atom occurs only in B, the A-preferring sanitizer still keeps it.
        mixed = folder / "mixed_altloc.pdb"
        mixed.write_text(atom(1, x=0, name="C1", alt="A") + "\n" + atom(2, x=10, name="C1", alt="B") + "\n" + atom(3, x=11, name="C2", alt="B") + "\n")
        sanitize_ligand_pdb_file(mixed, mixed)
        result["mixed_altloc_x"] = [float(line[30:38]) for line in mixed.read_text().splitlines() if line.startswith("HETATM")]
        assert result["mixed_altloc_x"] == [0.0, 11.0]

        assets = folder / "assets"
        assets.mkdir()
        (assets / "1ABC.pdbqt").write_text("original")
        (assets / "1ABC_cleaned.pdbqt").write_text("cleaned")
        builder = DockingProjectBuilder(SimpleNamespace())
        index = builder._index_assets(assets, {".pdbqt"})
        chosen = builder._resolve_asset("1ABC_cleaned", index, "receptor")
        result["resolved_explicit_cleaned_stem"] = chosen.name
        # An exact cleaned stem currently resolves correctly; PDB-only identity chooses original.
        result["resolved_pdb_id"] = builder._resolve_asset("1ABC", index, "receptor").name
        assert result["resolved_pdb_id"] == "1ABC.pdbqt"

        class Pipeline:
            def fetch_pdb(self, *args): return "fake.pdb"
            def enumerate_hetatms(self, *args): return [("LIG", "A", 1, None)], ["LIG", "HOH"]
            def save_hetatm_as_pdb(self, *args, **kwargs): return "ligand.pdb"
            def clean_pdb(self, _file, remove, **kwargs): self.removed = list(remove); return "clean.pdb"
            def analyze_pocket_properties(self, *args): return {}
            def generate_summary_report(self, *args): return "report"
        namespace = {"NON_ANCHOR_DEFAULT": set(), "_select_ligand_instance": lambda *a, **kw: ("LIG", "A", 1), "extract_residue_level_coordinates": lambda *a: {"overall_center": [0,0,0], "num_interacting_residues": 1, "num_interacting_atoms": 1}}
        run_cli = extracted("cli_pipeline.py", "run_single_pdb_cli", namespace)
        fake = Pipeline()
        with contextlib.redirect_stdout(io.StringIO()):
            okay = run_cli(fake, "1ABC", {"cleaning": {"default": "common"}})
        result["default_cleaning"] = {"success": okay, "removed": fake.removed, "selected_ligand_retained": "LIG" not in fake.removed}
        assert result["default_cleaning"]["selected_ligand_retained"] and okay

        # Call the actual pocket method against an empty structural iterator.
        pocket = extracted("core_pipeline.py", "analyze_pocket_properties", {"np": np})
        fake_structure = SimpleNamespace(_load_structure=lambda *a: [])
        with contextlib.redirect_stdout(io.StringIO()):
            pocket_result = pocket(fake_structure, "empty.pdb", np.zeros(3))
        result["empty_pocket"] = {k: pocket_result[k] for k in ("pocket_volume_A3", "druggability_score", "druggability_interpretation")}
        assert result["empty_pocket"]["druggability_interpretation"] == "Good"

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
