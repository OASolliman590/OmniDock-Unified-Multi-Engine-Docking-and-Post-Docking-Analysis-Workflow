"""Focused standard-library receptor and PDB-selection contract tests."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from docking.preparation import receptor_preparation
from docking.preparation.structure_contract import select_pdb_lines, write_selection


def atom(
    serial: int,
    name: str = "CA",
    *,
    x: float = 0.0,
    y: float = 0.0,
    z: float = 0.0,
    element: str = "C",
    record: str = "ATOM",
    resname: str = "ALA",
    chain: str = "A",
    resid: int = 1,
    alt: str = "",
) -> str:
    return (
        f"{record:<6}{serial:5d} {name:>4}{alt:1}{resname:>3} {chain:1}"
        f"{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}"
        f"          {element:>2}"
    )


def pdbqt_atom(serial: int, name: str = "CA", *, x: float = 0.0, element: str = "C", resname: str = "ALA", chain: str = "A", resid: int = 1) -> str:
    # Keep coordinates in the fixed PDB/PDBQT columns; the final token is the
    # AutoDock atom type consumed by the current receptor parser.
    return atom(serial, name, x=x, element=element, resname=resname, chain=chain, resid=resid)[:66] + f"    {0.0:6.3f} {element}"


def _fixed_record(record: str, fields: dict[tuple[int, int], str]) -> str:
    chars = [" "] * 80
    chars[: len(record)] = record
    for (start, end), value in fields.items():
        text = str(value)
        chars[start:end] = list(text[: end - start].ljust(end - start))
    return "".join(chars)


def link_line(atom1: str, res1: str, chain1: str, seq1: int, atom2: str, res2: str, chain2: str, seq2: int) -> str:
    return _fixed_record(
        "LINK  ",
        {
            (12, 16): f"{atom1:>4}",
            (17, 20): f"{res1:>3}",
            (21, 22): chain1,
            (22, 26): f"{seq1:4d}",
            (42, 46): f"{atom2:>4}",
            (47, 50): f"{res2:>3}",
            (51, 52): chain2,
            (52, 56): f"{seq2:4d}",
        },
    )


def ssbond_line(res1: str, chain1: str, seq1: int, res2: str, chain2: str, seq2: int) -> str:
    return _fixed_record(
        "SSBOND ",
        {
            (7, 10): "1",
            (11, 14): f"{res1:>3}",
            (15, 16): chain1,
            (17, 21): f"{seq1:4d}",
            (25, 28): f"{res2:>3}",
            (29, 30): chain2,
            (31, 35): f"{seq2:4d}",
            (73, 78): "2.03",
        },
    )


def _write_lines(path: Path, lines: list[str]) -> Path:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_selection_preserves_retained_chain_and_connectivity_records(tmp_path: Path):
    source = _write_lines(
        tmp_path / "complex.pdb",
        [
            atom(1, "SG", element="S", resname="CYS", resid=1),
            atom(2, "CA", x=2.0, resname="CYS", resid=2),
            atom(3, "SG", x=3.0, element="S", resname="CYS", resid=2),
            "TER       4      CYS A   2",
            link_line("SG", "CYS", "A", 1, "CA", "CYS", "A", 2),
            ssbond_line("CYS", "A", 1, "CYS", "A", 2),
            "CONECT    1    2",
            "CONECT    2    1    3",
        ],
    )

    output, _ = select_pdb_lines(source.read_text().splitlines())
    records = [line[:6].strip() for line in output]
    assert "TER" in records
    assert "LINK" in records
    assert "SSBOND" in records
    assert output.index("TER       4      CYS A   2") == 3
    assert output[-3:] == ["CONECT    1    2", "CONECT    2    1    3", "END"]


def test_selection_rejects_link_crossing_selected_boundary(tmp_path: Path):
    source = _write_lines(
        tmp_path / "covalent.pdb",
        [
            atom(1, "C1", record="HETATM", resname="LIG", chain="L", resid=7),
            atom(2, "SG", element="S", resname="CYS", chain="A", resid=1),
            link_line("C1", "LIG", "L", 7, "SG", "CYS", "A", 1),
        ],
    )
    with pytest.raises(ValueError, match="LINK.*selection"):
        select_pdb_lines(source.read_text().splitlines(), residue=("L", "7", "", "LIG"))


def test_selection_rejects_conect_crossing_selected_boundary(tmp_path: Path):
    source = _write_lines(
        tmp_path / "linked.pdb",
        [
            atom(1, "C1", record="HETATM", resname="LIG", chain="L", resid=7),
            atom(2, "SG", element="S", resname="CYS", chain="A", resid=1),
            "CONECT    1    2",
        ],
    )
    with pytest.raises(ValueError, match="CONECT.*selection"):
        select_pdb_lines(source.read_text().splitlines(), residue=("L", "7", "", "LIG"))


def test_selection_rejects_disulfide_crossing_selected_boundary(tmp_path: Path):
    source = _write_lines(
        tmp_path / "disulfide.pdb",
        [
            atom(1, "SG", element="S", resname="CYS", chain="A", resid=1),
            atom(2, "SG", element="S", resname="CYS", chain="A", resid=2),
            ssbond_line("CYS", "A", 1, "CYS", "A", 2),
        ],
    )
    with pytest.raises(ValueError, match="SSBOND.*selected"):
        select_pdb_lines(source.read_text().splitlines(), residue=("A", "1", "", "CYS"))


def test_ligand_selection_does_not_inherit_same_chain_polymer_ter(tmp_path: Path):
    source = _write_lines(
        tmp_path / "ligand_after_ter.pdb",
        [
            atom(1, "CA", record="ATOM", resname="ALA", chain="A", resid=1),
            "TER       2      ALA A   1",
            atom(3, "C1", record="HETATM", resname="LIG", chain="A", resid=2),
        ],
    )
    output, _ = select_pdb_lines(source.read_text().splitlines(), residue=("A", "2", "", "LIG"))
    assert [line[:6].strip() for line in output] == ["HETATM", "END"]


def test_retained_link_altloc_is_normalized_with_selected_atom(tmp_path: Path):
    source = _write_lines(
        tmp_path / "alt_link.pdb",
        [
            atom(1, "C1", alt="A", record="HETATM", resname="LIG", chain="A", resid=1),
            atom(2, "SG", alt="A", element="S", resname="CYS", chain="A", resid=2),
            link_line("C1", "LIG", "A", 1, "SG", "CYS", "A", 2),
        ],
    )
    # Put the selected locations into the fixed-column LINK record.
    link = link_line("C1", "LIG", "A", 1, "SG", "CYS", "A", 2)
    link = link[:16] + "A" + link[17:46] + "A" + link[47:]
    source.write_text("\n".join([atom(1, "C1", alt="A", record="HETATM", resname="LIG", chain="A", resid=1), atom(2, "SG", alt="A", element="S", resname="CYS", chain="A", resid=2), link]) + "\n", encoding="utf-8")
    output, _ = select_pdb_lines(source.read_text().splitlines())
    retained = next(line for line in output if line.startswith("LINK"))
    assert retained[16] == " " and retained[46] == " "


def test_write_selection_frame_ids_are_source_scoped(tmp_path: Path):
    first = _write_lines(tmp_path / "first.pdb", [atom(1, "CA", resname="ALA")])
    second = _write_lines(tmp_path / "second.pdb", [atom(1, "CA", resname="GLY")])
    first_meta = write_selection(first, tmp_path / "first_out.pdb")
    second_meta = write_selection(second, tmp_path / "second_out.pdb")
    assert first_meta["reference_frame_id"].startswith("source-pdb-v1:")
    assert first_meta["reference_frame_id"] != second_meta["reference_frame_id"]
    assert first_meta["frame_source_sha256"] != second_meta["frame_source_sha256"]


def test_receptor_conservation_allows_serial_renumbering_and_records_map(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    source = _write_lines(tmp_path / "receptor.pdb", [atom(1, "CA"), atom(2, "CB", x=1.0)])
    output = tmp_path / "receptor.pdbqt"

    def fake_run(command, **kwargs):
        prepared = Path(command[command.index("-p") + 1])
        prepared.write_text("\n".join([pdbqt_atom(91, "CA"), pdbqt_atom(92, "CB", x=1.0), "ROOT", "ENDROOT", "TORSDOF 0"]) + "\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(receptor_preparation.subprocess, "run", fake_run)
    report = receptor_preparation.prepare_receptor(source, output, {"preparation": {}})
    assert report["heavy_atom_conservation"] == "passed"
    assert [(row["source_serial"], row["output_serial"]) for row in report["atom_map"]] == [(1, 91), (2, 92)]


def test_receptor_conservation_rejects_swapped_atom_labels(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    source = _write_lines(tmp_path / "receptor.pdb", [atom(1, "CA"), atom(2, "CB", x=1.0)])
    output = tmp_path / "receptor.pdbqt"

    def fake_run(command, **kwargs):
        prepared = Path(command[command.index("-p") + 1])
        prepared.write_text("\n".join([pdbqt_atom(91, "CB"), pdbqt_atom(92, "CA", x=1.0), "ROOT", "ENDROOT", "TORSDOF 0"]) + "\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(receptor_preparation.subprocess, "run", fake_run)
    with pytest.raises(ValueError, match="identity"):
        receptor_preparation.prepare_receptor(source, output, {"preparation": {}})


def test_receptor_conservation_rejects_moved_atom_at_record_precision(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    source = _write_lines(tmp_path / "receptor.pdb", [atom(1, "CA")])
    output = tmp_path / "receptor.pdbqt"

    def fake_run(command, **kwargs):
        prepared = Path(command[command.index("-p") + 1])
        prepared.write_text("\n".join([pdbqt_atom(91, "CA", x=0.010), "ROOT", "ENDROOT", "TORSDOF 0"]) + "\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(receptor_preparation.subprocess, "run", fake_run)
    with pytest.raises(ValueError, match="coordinate"):
        receptor_preparation.prepare_receptor(source, output, {"preparation": {}})


def test_receptor_conservation_ignores_polar_hydrogen_atom_type(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    source = _write_lines(tmp_path / "receptor.pdb", [atom(1, "N", element="N"), atom(2, "H", element="H")])
    output = tmp_path / "receptor.pdbqt"

    def fake_run(command, **kwargs):
        prepared = Path(command[command.index("-p") + 1])
        prepared.write_text("\n".join([pdbqt_atom(91, "N", element="N"), pdbqt_atom(92, "H", element="HD"), "ROOT", "ENDROOT", "TORSDOF 0"]) + "\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(receptor_preparation.subprocess, "run", fake_run)
    report = receptor_preparation.prepare_receptor(source, output, {"preparation": {}})
    assert [row["atom_name"] for row in report["atom_map"]] == ["N"]


def test_receptor_rejects_ambiguous_missing_ca_element(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    source = _write_lines(tmp_path / "receptor.pdb", [atom(1, "CA", element="")])
    output = tmp_path / "receptor.pdbqt"
    monkeypatch.setattr(receptor_preparation.subprocess, "run", lambda *args, **kwargs: pytest.fail("backend should not run"))
    with pytest.raises(ValueError, match="unambiguous element"):
        receptor_preparation.prepare_receptor(source, output, {"preparation": {}})


def test_receptor_rejects_duplicate_source_serials(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    source = _write_lines(tmp_path / "receptor.pdb", [atom(1, "CA"), atom(1, "CB", x=1.0)])
    output = tmp_path / "receptor.pdbqt"
    monkeypatch.setattr(receptor_preparation.subprocess, "run", lambda *args, **kwargs: pytest.fail("backend should not run"))
    with pytest.raises(ValueError, match="serial identifiers must be unique"):
        receptor_preparation.prepare_receptor(source, output, {"preparation": {}})


def test_pdb2pqr_titration_is_explicit_and_not_final_certification(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    source = _write_lines(tmp_path / "receptor.pdb", [atom(1, "CA")])
    output = tmp_path / "receptor.pdbqt"
    commands: list[list[str]] = []

    def fake_run(command, **kwargs):
        commands.append(list(command))
        if command[0] == "pdb2pqr30":
            Path(command[-1]).write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
        elif command[0] == "obabel":
            Path(command[command.index("-O") + 1]).write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
        else:
            prepared = Path(command[command.index("-p") + 1])
            prepared.write_text("\n".join([pdbqt_atom(91, "CA"), "ROOT", "ENDROOT", "TORSDOF 0"]) + "\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(receptor_preparation.subprocess, "run", fake_run)
    report = receptor_preparation.prepare_receptor(
        source,
        output,
        {"preparation": {"receptor_use_pdb2pqr": True, "ph": 6.5}},
    )
    pqr = commands[0]
    assert "--titration-state-method" in pqr
    assert pqr[pqr.index("--titration-state-method") + 1] == "propka"
    assert "--keep-chain" in pqr
    assert report["upstream_titration_status"] == "pdb2pqr_propka_applied"
    assert report["final_pH_validated"] is False
    assert report["protonation_status"] == "backend_state_not_ph_validated"


def test_receptor_rejects_nonfinite_pH_before_backend(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    source = _write_lines(tmp_path / "receptor.pdb", [atom(1, "CA")])
    output = tmp_path / "receptor.pdbqt"
    called = False

    def fake_run(command, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("backend should not run")

    monkeypatch.setattr(receptor_preparation.subprocess, "run", fake_run)
    with pytest.raises(ValueError, match="finite"):
        receptor_preparation.prepare_receptor(source, output, {"preparation": {"ph": float("nan")}})
    assert not called
