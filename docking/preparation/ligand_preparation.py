from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from docking.models import (
    DEFAULT_PREPARATION_PH,
    LIGAND_PREPARATION_PROFILES,
    MAX_PREPARATION_PH,
    MIN_PREPARATION_PH,
    normalize_engine_names,
    normalize_ligand_preparation_profile,
    resolve_effective_ligand_preparation_profile,
    validate_preparation_ph,
    validate_ligand_preparation_profile,
)


def _resolve_profile_alias(raw: str) -> str:
    normalized = normalize_ligand_preparation_profile(raw)
    return normalized if normalized in LIGAND_PREPARATION_PROFILES else "engine_aware_full"


def _run_command(
    command: list[str],
    env: Optional[Dict[str, str]] = None,
    cwd: Optional[Path] = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        env=env,
        cwd=str(cwd) if cwd else None,
    )


def _selected_profile() -> str:
    raw = str(
        os.environ.get("PDBWIZARD_LIGAND_PREP_PROFILE")
        or os.environ.get("PDBWIZARD_LIGAND_PREP_BACKEND")
        or "engine_aware_full"
    )
    return _resolve_profile_alias(raw)


def _selected_engines() -> Tuple[str, ...]:
    raw = str(os.environ.get("PDBWIZARD_SELECTED_ENGINES", "") or "").strip()
    if not raw:
        return ()
    engines = normalize_engine_names(raw.split(","))
    return tuple(engines)


def _resolve_prepare_ligand4_script() -> Optional[str]:
    repo_root = Path(__file__).resolve().parents[2]
    candidates = [
        os.environ.get("AUTODOCKTOOLS_PREPARE_LIGAND4"),
        os.environ.get("ADT_PREPARE_LIGAND4"),
        shutil.which("prepare_ligand4.py"),
        str(repo_root / ".workflow" / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"),
        str(repo_root / "tools" / "autodocktools-prepare-py3k" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"),
    ]
    for candidate in candidates:
        if not candidate:
            continue
        script = Path(candidate).expanduser()
        if script.exists():
            return str(script.resolve())
        if candidate.endswith(".py"):
            continue
        possible_py = Path(f"{candidate}.py").expanduser()
        if possible_py.exists():
            return str(possible_py.resolve())
    return None


def _is_zero_molecule_obabel_result(result: Optional[subprocess.CompletedProcess[str]]) -> bool:
    if result is None:
        return False
    stderr = str(result.stderr or "").lower()
    stdout = str(result.stdout or "").lower()
    return "0 molecules converted" in stderr or "0 molecules converted" in stdout


def _file_has_meaningful_content(path: Path) -> bool:
    file_path = Path(path)
    if not file_path.exists() or file_path.stat().st_size == 0:
        return False
    try:
        content = file_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return file_path.stat().st_size > 0
    return bool(content.strip())


def _mol2_has_atoms(path: Path) -> bool:
    file_path = Path(path)
    if not _file_has_meaningful_content(file_path):
        return False
    try:
        content = file_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return False
    return "@<TRIPOS>ATOM" in content


def _prepare_with_autodocktools(input_sdf: Path, output_pdbqt: Path) -> Dict[str, object]:
    script = _resolve_prepare_ligand4_script()
    if not script:
        raise FileNotFoundError(
            "AutoDockTools prepare_ligand4.py was not found. "
            "Set AUTODOCKTOOLS_PREPARE_LIGAND4 (or ADT_PREPARE_LIGAND4) to enable this fallback."
        )

    python_exe = os.environ.get("AUTODOCKTOOLS_PYTHON") or os.environ.get("ADT_PYTHON") or sys.executable
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    script_path = Path(script).expanduser().resolve()
    adt_repo_root = script_path.parent.parent.parent
    env = dict(os.environ)
    existing_pythonpath = str(env.get("PYTHONPATH", "") or "").strip()
    env["PYTHONPATH"] = str(adt_repo_root) if not existing_pythonpath else f"{adt_repo_root}:{existing_pythonpath}"
    input_file = Path(input_sdf).expanduser().resolve()
    output_file = Path(output_pdbqt).expanduser().resolve()
    work_dir = input_file.parent
    mol2_input = work_dir / f"{input_file.stem}_adt_input.mol2"
    local_output = work_dir / f"{input_file.stem}_adt.pdbqt"

    ligand_arg = input_file.name
    try:
        mol2_result = _run_command(["obabel", str(input_file), "-O", str(mol2_input)])
        if not _is_zero_molecule_obabel_result(mol2_result) and _mol2_has_atoms(mol2_input):
            ligand_arg = mol2_input.name
    except Exception:
        # Fall back to the original SDF when conversion is unavailable.
        ligand_arg = input_file.name

    command = [
        python_exe,
        str(script_path),
        "-l",
        ligand_arg,
        "-o",
        local_output.name,
        "-p",
        "Pd",
    ]
    result = _run_command(command, env=env, cwd=work_dir)
    if not local_output.exists():
        raise FileNotFoundError(f"AutoDockTools did not produce expected output file: {local_output}")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if local_output != output_file:
        shutil.copy2(local_output, output_file)
    return {
        "preparation_method": "autodocktools",
        "autodocktools_script": script,
        "autodocktools_python": python_exe,
        "autodocktools_command": command,
        "autodocktools_workdir": str(work_dir),
        "autodocktools_stdout": (result.stdout or "").strip(),
        "autodocktools_stderr": (result.stderr or "").strip(),
    }


def _load_rdkit() -> Tuple[Optional[Any], Optional[Any]]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return None, None
    return Chem, AllChem


def _load_rdkit_molecule(input_path: Path) -> Any:
    Chem, _ = _load_rdkit()
    if Chem is None:
        raise ImportError("RDKit is not available")

    suffix = input_path.suffix.lower()
    if suffix == ".sdf":
        supplier = Chem.SDMolSupplier(str(input_path), removeHs=False)
        for mol in supplier:
            if mol is not None:
                return mol
        return None
    if suffix == ".mol":
        return Chem.MolFromMolFile(str(input_path), removeHs=False)
    if suffix == ".mol2":
        return Chem.MolFromMol2File(str(input_path), removeHs=False)
    if suffix == ".pdb":
        return Chem.MolFromPDBFile(str(input_path), removeHs=False)
    raise ValueError(f"Unsupported ligand format for RDKit normalization: {input_path.suffix}")


def _has_3d_coordinates(mol: Any) -> bool:
    if mol is None or mol.GetNumConformers() == 0:
        return False
    conformer = mol.GetConformer()
    if conformer.Is3D():
        return True
    positions = conformer.GetPositions()
    if len(positions) == 0:
        return False
    return max(abs(float(position[2])) for position in positions) > 1e-3


def _embed_and_optimize_3d(mol: Any, all_chem: Any) -> str:
    params = all_chem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    status = all_chem.EmbedMolecule(mol, params)
    if status != 0:
        params.useRandomCoords = True
        status = all_chem.EmbedMolecule(mol, params)
    if status != 0:
        raise ValueError("RDKit could not embed a 3D conformer")

    try:
        if all_chem.MMFFHasAllMoleculeParams(mol):
            all_chem.MMFFOptimizeMolecule(mol)
            return "MMFF94"
        all_chem.UFFOptimizeMolecule(mol)
        return "UFF"
    except Exception:
        return ""


def _normalize_with_rdkit(input_path: Path, output_sdf: Path) -> Dict[str, object]:
    Chem, all_chem = _load_rdkit()
    if Chem is None or all_chem is None:
        raise ImportError("RDKit is not available")

    molecule = _load_rdkit_molecule(input_path)
    if molecule is None:
        raise ValueError(f"RDKit could not parse ligand file: {input_path}")

    Chem.SanitizeMol(molecule)
    had_3d = _has_3d_coordinates(molecule)
    molecule = Chem.AddHs(molecule, addCoords=had_3d)

    optimized_forcefield = ""
    if not had_3d:
        optimized_forcefield = _embed_and_optimize_3d(molecule, all_chem)

    molecule.SetProp("_Name", molecule.GetProp("_Name") if molecule.HasProp("_Name") else input_path.stem)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_sdf))
    try:
        writer.write(molecule)
    finally:
        writer.close()

    return {
        "normalization_backend": "rdkit",
        "source_file": str(input_path),
        "normalized_sdf": str(output_sdf),
        "had_3d_input": had_3d,
        "generated_3d": not had_3d,
        "optimized_forcefield": optimized_forcefield,
    }


class ProtonationStateError(RuntimeError):
    """A protonation policy could not be satisfied. ``reason`` is machine-readable (Spec 034 R2)."""

    def __init__(self, reason: str, details: str) -> None:
        super().__init__(f"{reason}: {details}")
        self.reason = reason
        self.details = details


def _rdkit_or_error() -> Tuple[Any, Any]:
    Chem, all_chem = _load_rdkit()
    if Chem is None or all_chem is None:
        raise ProtonationStateError("rdkit_unavailable", "RDKit is required to measure or build protonation states")
    return Chem, all_chem


# ---------------------------------------------------------------- Spec 034 R2b: measurement

def canonical_microspecies(mol: Any) -> str:
    """Canonical SMILES with hydrogens removed and charges kept.

    Stereo is not part of the protonation-state contract: 3D perception can add spurious
    marks on charged amines, and the state map carries no stereo. The stereo-free comparison
    is recorded in provenance as ``stereo_compared: false``.
    """
    Chem, _ = _rdkit_or_error()
    plain = Chem.RemoveHs(mol)
    Chem.RemoveStereochemistry(plain)
    return Chem.MolToSmiles(plain)


def canonical_smiles_from_text(smiles: str) -> str:
    Chem, _ = _rdkit_or_error()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ProtonationStateError("explicit_state_invalid_smiles", f"not a valid SMILES: {smiles}")
    return canonical_microspecies(mol)


def measure_molecule_state(mol: Any) -> Dict[str, object]:
    """Measured net formal charge, charged atoms and microspecies SMILES of one RDKit molecule."""
    Chem, _ = _rdkit_or_error()
    charged_atoms = [
        {"element": atom.GetSymbol(), "index": int(atom.GetIdx()) + 1, "charge": int(atom.GetFormalCharge())}
        for atom in mol.GetAtoms()
        if atom.GetFormalCharge()
    ]
    return {
        "net_formal_charge": int(Chem.GetFormalCharge(mol)),
        "charged_atoms": charged_atoms,
        "charged_atom_index_basis": "1-based atom order of the measured file",
        "microspecies_smiles": canonical_microspecies(mol),
        "stereo_compared": False,
        "hydrogen_count": sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1),
        "heavy_atom_count": sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1),
    }


def measure_ligand_file(path: Path) -> Dict[str, object]:
    """Measure the protonation state of a ligand file (SDF, MOL, MOL2 or PDB) with RDKit."""
    mol = _load_rdkit_molecule(Path(path))
    if mol is None:
        raise ValueError(f"RDKit could not read a ligand state from {Path(path).name}")
    state = measure_molecule_state(mol)
    state["measured_file"] = Path(path).name
    return state


def _try_measure_ligand(path: Path) -> Tuple[Optional[Dict[str, object]], str]:
    try:
        return measure_ligand_file(path), ""
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"


def measure_pdbqt_remark_smiles(path: Path) -> Optional[Dict[str, object]]:
    """Measure the docking input: the ``REMARK SMILES`` written by Meeko. None when absent."""
    smiles = ""
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX"):
            smiles = line[len("REMARK SMILES "):].strip()
            break
    if not smiles:
        return None
    Chem, _ = _rdkit_or_error()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ProtonationStateError("pdbqt_smiles_unreadable", f"REMARK SMILES could not be parsed: {smiles}")
    return {
        "smiles": smiles,
        "net_formal_charge": int(Chem.GetFormalCharge(mol)),
        "microspecies_smiles": canonical_microspecies(mol),
        "stereo_compared": False,
    }


def verify_microspecies_matches(state: Dict[str, object], smiles: str, *, stage: str) -> None:
    """Fail with ``explicit_state_mismatch`` unless the measured state is the requested microspecies."""
    expected = canonical_smiles_from_text(smiles)
    measured = str(state.get("microspecies_smiles") or "")
    if measured != expected:
        raise ProtonationStateError(
            "explicit_state_mismatch",
            f"{stage} microspecies {measured!r} is not the requested {expected!r}",
        )


def _sdf_hydrogen_count(path: Path) -> int:
    """Count hydrogen atom records in the first V2000 molecule of an SDF (no RDKit needed)."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 4:
        return 0
    try:
        atom_count = int(lines[3][:3])
    except ValueError:
        return 0
    symbols = [line.split()[3] for line in lines[4:4 + atom_count] if len(line.split()) > 3]
    return sum(1 for symbol in symbols if symbol == "H")


# ---------------------------------------------------------------- Spec 034 R2c: policies

PROTONATION_PH_MODEL = "ph_model"
PROTONATION_EXPLICIT_STATE = "explicit_state"
PROTONATION_AS_INPUT = "as_input"
PROTONATION_POLICIES = (PROTONATION_PH_MODEL, PROTONATION_EXPLICIT_STATE, PROTONATION_AS_INPUT)
DEFAULT_PROTONATION_POLICY = PROTONATION_PH_MODEL
# Profiles that convert the input file directly and never run the Open Babel pH model.
DIRECT_PREPARATION_PROFILES = ("meeko_only", "autodocktools_only")
PH_SOURCE_USER_ENTERED = "user_entered"
PH_SOURCE_CONFIG_FILE = "config_file"
PH_SOURCE_MODULE_DEFAULT = "module_default"


def _selected_protonation_policy() -> str:
    raw = str(os.environ.get("PDBWIZARD_LIGAND_PREP_PROTONATION_POLICY") or DEFAULT_PROTONATION_POLICY).strip().lower()
    if raw not in PROTONATION_POLICIES:
        raise ProtonationStateError(
            "invalid_protonation_policy",
            f"{raw!r} is not one of: {', '.join(PROTONATION_POLICIES)}",
        )
    return raw


def _selected_ph_with_source() -> Tuple[float, str]:
    """pH for the Open Babel pH model and where it came from.

    The source is ``user_entered`` when the CLI, the interactive prompt or a config file set it
    (PDBWIZARD_LIGAND_PREP_PH_SOURCE). Without any value the module default is used and recorded
    as ``module_default``; the pipeline entry points never reach that branch.
    """
    raw = str(os.environ.get("PDBWIZARD_LIGAND_PREP_PH", "") or "").strip()
    if not raw:
        return float(DEFAULT_PREPARATION_PH), PH_SOURCE_MODULE_DEFAULT
    ok, normalized, error = validate_preparation_ph(raw)
    if not ok:
        raise ProtonationStateError("invalid_ph", error)
    source = str(os.environ.get("PDBWIZARD_LIGAND_PREP_PH_SOURCE") or "").strip()
    if source not in (PH_SOURCE_USER_ENTERED, PH_SOURCE_CONFIG_FILE):
        raise ProtonationStateError(
            "ph_source_missing",
            "PDBWIZARD_LIGAND_PREP_PH_SOURCE must be 'user_entered' or 'config_file' when a pH is set",
        )
    return float(normalized), source


def _selected_state_map_path() -> Optional[Path]:
    raw = str(os.environ.get("PDBWIZARD_LIGAND_PREP_STATE_MAP", "") or "").strip()
    return Path(raw).expanduser().resolve() if raw else None


def load_protonation_state_map(path: Path) -> Dict[str, Dict[str, object]]:
    """Read the per-ligand state map: CSV with a ``ligand`` column and ``smiles`` and/or ``net_charge``."""
    try:
        handle = Path(path).open(newline="", encoding="utf-8")
    except OSError as exc:
        raise ProtonationStateError("state_map_unreadable", str(exc)) from exc
    entries: Dict[str, Dict[str, object]] = {}
    with handle:
        reader = csv.DictReader(handle)
        columns = {str(name or "").strip().lower() for name in (reader.fieldnames or [])}
        if "ligand" not in columns or not ({"smiles", "net_charge"} & columns):
            raise ProtonationStateError(
                "state_map_invalid",
                "the state map needs a 'ligand' column plus 'smiles' and/or 'net_charge' columns",
            )
        for row_number, row in enumerate(reader, start=2):
            fields = {str(key or "").strip().lower(): str(value or "").strip() for key, value in row.items()}
            ligand = fields.get("ligand", "")
            smiles = fields.get("smiles", "")
            net_raw = fields.get("net_charge", "")
            if not ligand:
                raise ProtonationStateError("state_map_invalid", f"row {row_number}: empty ligand")
            if ligand in entries:
                raise ProtonationStateError("state_map_invalid", f"row {row_number}: duplicate ligand {ligand!r}")
            if not smiles and not net_raw:
                raise ProtonationStateError("state_map_invalid", f"row {row_number}: {ligand!r} needs smiles or net_charge")
            smiles_net: Optional[int] = None
            if smiles:
                Chem, _ = _rdkit_or_error()
                parsed = Chem.MolFromSmiles(smiles)
                if parsed is None:
                    raise ProtonationStateError("state_map_invalid", f"row {row_number}: invalid SMILES for {ligand!r}")
                smiles_net = int(Chem.GetFormalCharge(parsed))
            net_charge: Optional[int] = None
            if net_raw:
                try:
                    net_charge = int(net_raw)
                except ValueError as exc:
                    raise ProtonationStateError(
                        "state_map_invalid",
                        f"row {row_number}: net_charge must be an integer, got {net_raw!r}",
                    ) from exc
            if smiles_net is not None and net_charge is not None and smiles_net != net_charge:
                raise ProtonationStateError(
                    "state_map_invalid",
                    f"row {row_number}: {ligand!r} SMILES net charge {smiles_net:+d} disagrees with net_charge {net_charge:+d}",
                )
            entries[ligand] = {
                "ligand": ligand,
                "smiles": smiles,
                "net_charge": smiles_net if smiles_net is not None else net_charge,
            }
    return entries


def _explicit_state_entry(source: Path) -> Dict[str, object]:
    state_map_path = _selected_state_map_path()
    if state_map_path is None:
        raise ProtonationStateError(
            "explicit_state_map_missing",
            "explicit_state requires a state map (--protonation-state-map)",
        )
    entries = load_protonation_state_map(state_map_path)
    for key in (source.name, source.stem):
        if key in entries:
            return entries[key]
    raise ProtonationStateError(
        "explicit_state_missing_entry",
        f"state map {state_map_path.name} has no row for ligand {source.name!r} (ligand = file name or stem)",
    )


def _neutral_heavy_skeleton(mol: Any) -> Any:
    Chem, _ = _rdkit_or_error()
    skeleton = Chem.RemoveHs(mol)
    for atom in skeleton.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
    Chem.SanitizeMol(skeleton)
    return skeleton


def _map_heavy_atoms_by_input_geometry(target: Any, source_heavy: Any) -> Tuple[Optional[List[int]], str, float]:
    """Target index for each input heavy atom, chosen by bond-length strain. Returns (map, reason, strain)."""
    Chem, _ = _rdkit_or_error()
    skeleton_target = _neutral_heavy_skeleton(target)
    skeleton_source = _neutral_heavy_skeleton(source_heavy)
    if Chem.MolToSmiles(skeleton_target) != Chem.MolToSmiles(skeleton_source):
        return None, "input heavy-atom graph differs from the requested SMILES", 0.0
    matches = skeleton_target.GetSubstructMatches(skeleton_source, uniquify=False, maxMatches=2000)
    if not matches:
        return None, "no graph correspondence between input and requested SMILES", 0.0
    conformer = source_heavy.GetConformer()
    positions = [conformer.GetAtomPosition(index) for index in range(source_heavy.GetNumAtoms())]
    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in skeleton_target.GetBonds()]
    best: Optional[Tuple[float, List[int]]] = None
    for match in matches:  # match[s] = target index of input heavy atom s
        inverse = {target_index: source_index for source_index, target_index in enumerate(match)}
        strain = sum(
            abs(positions[inverse[a]].Distance(positions[inverse[b]]) - 1.45)
            for a, b in bonds
        )
        if best is None or strain < best[0]:
            best = (strain, list(match))
    assert best is not None
    return best[1], "input heavy-atom coordinates mapped by bond geometry", float(best[0])


def build_explicit_state_sdf(source: Path, smiles: str, output_sdf: Path) -> Dict[str, object]:
    """Build the requested microspecies with explicit hydrogens (Spec 034 R2c, explicit_state).

    Heavy-atom coordinates come from the input when the input is 3D and its heavy-atom graph is
    the requested SMILES (heavy atoms are copied exactly; hydrogens are placed on them). Otherwise
    coordinates are generated with ETKDG and MMFF/UFF, and the reason is recorded.
    """
    Chem, all_chem = _rdkit_or_error()
    target = Chem.MolFromSmiles(smiles)
    if target is None:
        raise ProtonationStateError("explicit_state_invalid_smiles", f"not a valid SMILES: {smiles}")

    try:
        raw_source = _load_rdkit_molecule(Path(source))
    except Exception as exc:
        raise ProtonationStateError("explicit_state_input_unreadable", f"{Path(source).name}: {exc}") from exc
    if raw_source is None:
        raise ProtonationStateError("explicit_state_input_unreadable", f"RDKit could not read {Path(source).name}")
    source_heavy = Chem.RemoveHs(raw_source)
    # The requested microspecies must be a protonation state of THIS ligand: same heavy-atom graph.
    # A different molecule is refused rather than built from its own SMILES.
    if Chem.MolToSmiles(_neutral_heavy_skeleton(target)) != Chem.MolToSmiles(_neutral_heavy_skeleton(source_heavy)):
        raise ProtonationStateError(
            "explicit_state_mismatch",
            f"the requested SMILES is not a protonation state of {Path(source).name} (different heavy-atom graph)",
        )
    input_is_3d = _has_3d_coordinates(source_heavy)

    mol = Chem.Mol(target)
    coordinate_reason = ""
    strain = None
    if input_is_3d:
        match, coordinate_reason, strain = _map_heavy_atoms_by_input_geometry(target, source_heavy)
    else:
        match = None
        coordinate_reason = "input has no 3D coordinates"

    strain = strain if match is not None else None
    if match is not None:
        conformer = Chem.Conformer(mol.GetNumAtoms())
        source_conformer = source_heavy.GetConformer()
        for source_index, target_index in enumerate(match):
            conformer.SetAtomPosition(target_index, source_conformer.GetAtomPosition(source_index))
        conformer.Set3D(True)
        mol.RemoveAllConformers()
        mol.AddConformer(conformer, assignId=True)
        mol = Chem.AddHs(mol, addCoords=True)
        coordinates_source = "input_heavy_atoms"
        forcefield = ""
    else:
        mol = Chem.AddHs(mol)
        try:
            forcefield = _embed_and_optimize_3d(mol, all_chem)
        except ValueError as exc:
            raise ProtonationStateError("explicit_state_embedding_failed", str(exc)) from exc
        coordinates_source = "generated_rdkit_etkdg"

    mol.SetProp("_Name", Path(source).stem)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_sdf))
    try:
        writer.write(mol)
    finally:
        writer.close()
    return {
        "normalization_backend": "rdkit_explicit_state",
        "source_file": str(source),
        "normalized_sdf": str(output_sdf),
        "state_source": PROTONATION_EXPLICIT_STATE,
        "explicit_smiles": smiles,
        "coordinates_source": coordinates_source,
        "coordinate_reason": coordinate_reason,
        "coordinate_bond_strain_a": round(strain, 4) if strain is not None else None,
        "had_3d_input": input_is_3d,
        "generated_3d": coordinates_source != "input_heavy_atoms",
        "hydrogens_added": True,
        "protonation_step_ran": False,
        "protonation_applied": False,
        "optimized_forcefield": forcefield,
        "output_state": measure_ligand_file(output_sdf),
    }


def prepare_as_input_sdf(source: Path, output_sdf: Path) -> Dict[str, object]:
    """Keep the input's hydrogens and charges. Add hydrogens and 3D only when they are missing."""
    Chem, all_chem = _rdkit_or_error()
    mol = _load_rdkit_molecule(Path(source))
    if mol is None:
        raise ValueError(f"RDKit could not parse ligand file: {source}")
    had_hydrogens = any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())
    had_3d = _has_3d_coordinates(mol)
    if not had_hydrogens:
        mol = Chem.AddHs(mol, addCoords=had_3d)
    forcefield = ""
    generated_3d = False
    if not had_3d:
        forcefield = _embed_and_optimize_3d(mol, all_chem)
        generated_3d = True
    mol.SetProp("_Name", mol.GetProp("_Name") if mol.HasProp("_Name") else Path(source).stem)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_sdf))
    try:
        writer.write(mol)
    finally:
        writer.close()
    return {
        "normalization_backend": "rdkit_as_input",
        "source_file": str(source),
        "normalized_sdf": str(output_sdf),
        "state_source": PROTONATION_AS_INPUT,
        "had_3d_input": had_3d,
        "generated_3d": generated_3d,
        "hydrogens_added": not had_hydrogens,
        "protonation_step_ran": False,
        "protonation_applied": False,
        "optimized_forcefield": forcefield,
        "output_state": measure_ligand_file(output_sdf),
    }


def _normalize_with_openbabel(
    input_path: Path,
    output_sdf: Path,
    *,
    protonation_ph: float = 7.4,
    partial_charge_model: str = "gasteiger",
) -> Dict[str, object]:
    if not shutil.which("obabel"):
        raise FileNotFoundError("Open Babel (obabel) is required for ligand normalization fallback")

    with tempfile.TemporaryDirectory(prefix=f"ligand_normalize_{input_path.stem}_") as tmp_dir:
        temp_dir = Path(tmp_dir)
        temp_source = temp_dir / f"{input_path.stem}_source.sdf"
        temp_ph = temp_dir / f"{input_path.stem}_ph.sdf"
        temp_3d = temp_dir / f"{input_path.stem}_3d.sdf"
        temp_h = temp_dir / f"{input_path.stem}_h.sdf"
        temp_charged = temp_dir / f"{input_path.stem}_charged.sdf"
        temp_min = temp_dir / f"{input_path.stem}_min.sdf"

        source_result = _run_command(["obabel", str(input_path), "-O", str(temp_source)])
        if _is_zero_molecule_obabel_result(source_result) or not _file_has_meaningful_content(temp_source):
            raise ValueError(f"OpenBabel could not parse ligand source into SDF: {input_path}")

        # Spec 034 R2a: the pH model runs on the hydrogen-free molecule. Open Babel does not
        # re-protonate a molecule that already has explicit hydrogens, so -d must come first.
        ph_result = _run_command(
            ["obabel", str(temp_source), "-O", str(temp_ph), "-d", "-p", f"{protonation_ph:.2f}"]
        )
        if _is_zero_molecule_obabel_result(ph_result) or not _file_has_meaningful_content(temp_ph):
            raise ValueError(f"OpenBabel pH model failed for: {input_path}")

        gen3d_result = _run_command(["obabel", str(temp_ph), "-O", str(temp_3d), "--gen3d"])
        if _is_zero_molecule_obabel_result(gen3d_result) or not _file_has_meaningful_content(temp_3d):
            raise ValueError(f"OpenBabel failed to generate 3D conformer for: {input_path}")

        # Explicit hydrogens are added without -p: the pH-model charge state is not re-derived.
        hydrogens_added_after_ph = False
        final_source = temp_3d
        if _sdf_hydrogen_count(temp_3d) == 0:
            addh_result = _run_command(["obabel", str(temp_3d), "-O", str(temp_h), "-h"])
            if _is_zero_molecule_obabel_result(addh_result) or not _file_has_meaningful_content(temp_h):
                raise ValueError(f"OpenBabel failed to add hydrogens for: {input_path}")
            hydrogens_added_after_ph = True
            final_source = temp_h

        partial_charge_applied = False
        if partial_charge_model:
            try:
                charge_result = _run_command(
                    [
                        "obabel",
                        str(final_source),
                        "-O",
                        str(temp_charged),
                        "--partialcharge",
                        partial_charge_model,
                    ]
                )
                if not _is_zero_molecule_obabel_result(charge_result) and _file_has_meaningful_content(temp_charged):
                    final_source = temp_charged
                    partial_charge_applied = True
            except subprocess.CalledProcessError:
                partial_charge_applied = False

        optimized_forcefield = ""
        for forcefield in ("MMFF94", "UFF"):
            try:
                min_result = _run_command(
                    [
                        "obabel",
                        str(final_source),
                        "-O",
                        str(temp_min),
                        "--minimize",
                        "--steps",
                        "250",
                        "--ff",
                        forcefield,
                    ]
                )
                if not _is_zero_molecule_obabel_result(min_result) and _file_has_meaningful_content(temp_min):
                    final_source = temp_min
                    optimized_forcefield = forcefield
                    break
            except subprocess.CalledProcessError:
                continue

        output_sdf.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(final_source, output_sdf)
        if not _file_has_meaningful_content(output_sdf):
            raise ValueError(f"OpenBabel normalization produced an empty SDF for: {input_path}")
        output_state, output_state_error = _try_measure_ligand(output_sdf)

    # Spec 034 R2b: applied only when the pH step ran AND the output was measured.
    return {
        "normalization_backend": "openbabel",
        "source_file": str(input_path),
        "normalized_sdf": str(output_sdf),
        "had_3d_input": False,
        "generated_3d": True,
        "state_source": PROTONATION_PH_MODEL,
        "protonation_ph": protonation_ph,
        "protonation_step_ran": True,
        "protonation_applied": bool(output_state is not None and output_state.get("microspecies_smiles")),
        "protonation_method": "openbabel_ph_model_on_hydrogen_free_input",
        "hydrogens_added_after_ph": hydrogens_added_after_ph,
        "output_state": output_state,
        "output_state_error": output_state_error,
        "partial_charge_model": partial_charge_model,
        "partial_charge_applied": partial_charge_applied,
        "optimized_forcefield": optimized_forcefield,
    }


def normalize_ligand_to_sdf(
    input_path: Path,
    output_sdf: Path,
    *,
    protonation_ph: float = 7.4,
    require_ph_model: bool = False,
) -> Dict[str, object]:
    source = Path(input_path).expanduser().resolve()
    destination = Path(output_sdf).expanduser().resolve()
    errors: list[str] = []

    # Enforce a single OpenBabel-first 3D normalization path across engines.
    try:
        return _normalize_with_openbabel(source, destination, protonation_ph=protonation_ph)
    except Exception as exc:
        errors.append(f"openbabel:{exc}")

    if require_ph_model:
        # The RDKit fallback applies no pH model. Under ph_model the ligand is refused, not left neutral.
        raise ProtonationStateError(
            "ph_model_unavailable",
            "; ".join(errors) + ". The pH model did not run, so the ligand was not prepared.",
        )
    payload = _normalize_with_rdkit(source, destination)
    payload["fallback_reasons"] = errors
    payload["protonation_ph"] = protonation_ph
    payload["state_source"] = "rdkit_fallback_no_ph_model"
    payload["protonation_step_ran"] = False
    payload["protonation_applied"] = False
    payload["output_state"], payload["output_state_error"] = _try_measure_ligand(destination)
    return payload


def _prepare_with_meeko(input_sdf: Path, output_pdbqt: Path) -> None:
    _run_command(["mk_prepare_ligand.py", "-i", str(input_sdf), "-o", str(output_pdbqt)])


def _prepare_with_openbabel_pdbqt(input_sdf: Path, output_pdbqt: Path) -> None:
    _run_command(
        [
            "obabel",
            str(input_sdf),
            "-O",
            str(output_pdbqt),
            "-h",
            "--partialcharge",
            "gasteiger",
        ]
    )


def _effective_profile(profile: str, selected_engines: Tuple[str, ...]) -> str:
    return resolve_effective_ligand_preparation_profile(profile, list(selected_engines))


def prepare_ligand_for_vina_family(input_path: Path, output_pdbqt: Path) -> Dict[str, object]:
    source = Path(input_path).expanduser().resolve()
    destination = Path(output_pdbqt).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)

    selected_profile = _selected_profile()
    selected_engines = _selected_engines()
    compatibility = validate_ligand_preparation_profile(selected_profile, list(selected_engines))
    if not compatibility.is_valid:
        raise RuntimeError("; ".join(compatibility.errors))
    effective_profile = _effective_profile(selected_profile, selected_engines)
    direct_profile = effective_profile in DIRECT_PREPARATION_PROFILES
    policy = _selected_protonation_policy()
    protonation_ph, ph_source = _selected_ph_with_source()
    state_entry = _explicit_state_entry(source) if policy == PROTONATION_EXPLICIT_STATE else None
    explicit_smiles = str(state_entry["smiles"]) if state_entry and state_entry.get("smiles") else ""

    with tempfile.TemporaryDirectory(prefix=f"ligand_prepare_{source.stem}_") as tmp_dir:
        normalized_sdf = Path(tmp_dir) / f"{source.stem}_normalized.sdf"
        input_state, input_state_error = _try_measure_ligand(source)
        normalization: Dict[str, object] = {}

        if policy == PROTONATION_AS_INPUT:
            normalization = prepare_as_input_sdf(source, normalized_sdf)
            state_source = PROTONATION_AS_INPUT
            conversion_input = normalized_sdf
        elif explicit_smiles:
            normalization = build_explicit_state_sdf(source, explicit_smiles, normalized_sdf)
            state_source = "explicit_state_smiles"
            conversion_input = normalized_sdf
        elif policy == PROTONATION_PH_MODEL and direct_profile:
            # Direct profiles convert the input as supplied. The Open Babel pH model is not run for them.
            state_source = "profile_input_unmodified"
            conversion_input = source
        else:
            normalization = normalize_ligand_to_sdf(
                source, normalized_sdf, protonation_ph=protonation_ph, require_ph_model=True
            )
            state_source = PROTONATION_PH_MODEL if policy == PROTONATION_PH_MODEL else "explicit_state_net_charge_checked"
            conversion_input = normalized_sdf

        if conversion_input == source:
            prepared_state, prepared_state_error = input_state, input_state_error
        else:
            prepared_state = normalization.get("output_state") if isinstance(normalization.get("output_state"), dict) else None
            prepared_state_error = str(normalization.get("output_state_error") or "")

        if policy == PROTONATION_EXPLICIT_STATE:
            if prepared_state is None:
                raise ProtonationStateError(
                    "explicit_state_unmeasured",
                    f"the prepared ligand could not be measured: {prepared_state_error}",
                )
            if explicit_smiles:
                verify_microspecies_matches(prepared_state, explicit_smiles, stage="prepared ligand")
            else:
                approved = int(state_entry["net_charge"])
                measured = int(prepared_state["net_formal_charge"])
                if measured != approved:
                    raise ProtonationStateError(
                        "explicit_state_net_charge_mismatch",
                        f"approved net charge {approved:+d} for {source.name}, measured {measured:+d} after the pH model",
                    )

        meeko_error = ""
        obabel_error = ""
        autodocktools_error = ""
        prep_method = ""
        autodocktools_summary: Optional[Dict[str, object]] = None

        if effective_profile == "meeko_only":
            try:
                _prepare_with_meeko(conversion_input, destination)
                prep_method = "meeko_direct"
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                meeko_error = str(exc)
                raise RuntimeError("Ligand preparation failed in Meeko-only mode.") from exc
        elif effective_profile == "autodocktools_only":
            try:
                autodocktools_summary = _prepare_with_autodocktools(conversion_input, destination)
                prep_method = "autodocktools_direct"
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                autodocktools_error = str(exc)
                raise RuntimeError("Ligand preparation failed in AutoDockTools-only mode.") from exc
        elif effective_profile == "openbabel_only":
            try:
                _prepare_with_openbabel_pdbqt(conversion_input, destination)
                prep_method = "openbabel_pdbqt"
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                obabel_error = str(exc)
                raise RuntimeError("Ligand preparation failed in OpenBabel-only mode.") from exc
        elif effective_profile == "openbabel_meeko":
            try:
                _prepare_with_meeko(conversion_input, destination)
                prep_method = "openbabel_then_meeko"
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                meeko_error = str(exc)
                raise RuntimeError("Ligand preparation failed in OpenBabel->Meeko mode.") from exc
        elif effective_profile == "openbabel_autodocktools":
            try:
                autodocktools_summary = _prepare_with_autodocktools(conversion_input, destination)
                prep_method = "openbabel_then_autodocktools"
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                autodocktools_error = str(exc)
                raise RuntimeError("Ligand preparation failed in OpenBabel->AutoDockTools mode.") from exc
        else:
            # openbabel_meeko_autodock
            try:
                _prepare_with_meeko(conversion_input, destination)
                prep_method = "openbabel_then_meeko"
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                meeko_error = str(exc)
                try:
                    autodocktools_summary = _prepare_with_autodocktools(conversion_input, destination)
                    prep_method = "openbabel_then_meeko_fallback_autodocktools"
                except (subprocess.CalledProcessError, FileNotFoundError) as adt_exc:
                    autodocktools_error = str(adt_exc)
                    raise RuntimeError(
                        "Ligand preparation failed in OpenBabel->Meeko->AutoDockTools mode."
                    ) from exc

        # Spec 034 R2b: measure the docking input. A PDBQT that changes the net charge fails here.
        pdbqt_state: Optional[Dict[str, object]] = None
        pdbqt_check = "not_measured_no_remark_smiles"
        canonical_match: Optional[bool] = None
        try:
            pdbqt_state = measure_pdbqt_remark_smiles(destination)
            if pdbqt_state is not None and prepared_state is not None:
                if int(pdbqt_state["net_formal_charge"]) != int(prepared_state["net_formal_charge"]):
                    raise ProtonationStateError(
                        "protonation_state_lost_in_conversion",
                        f"prepared ligand net charge {int(prepared_state['net_formal_charge']):+d}, "
                        f"PDBQT REMARK SMILES net charge {int(pdbqt_state['net_formal_charge']):+d} ({prep_method})",
                    )
                if explicit_smiles:
                    verify_microspecies_matches(pdbqt_state, explicit_smiles, stage="PDBQT REMARK SMILES")
                canonical_match = pdbqt_state["microspecies_smiles"] == prepared_state["microspecies_smiles"]
                pdbqt_check = "passed" if canonical_match else "passed_net_charge_canonical_differs"
            elif pdbqt_state is not None:
                pdbqt_check = "not_compared_prepared_state_unmeasured"
        except ProtonationStateError:
            destination.unlink(missing_ok=True)
            raise

    protonation_applied = bool(normalization.get("protonation_applied", False)) and prepared_state is not None
    protonation_block: Dict[str, object] = {
        "policy": policy,
        "state_source": state_source,
        "protonation_applied": protonation_applied,
        "ph_model_run": bool(normalization.get("protonation_step_ran", False)),
        "ph": protonation_ph if normalization.get("protonation_step_ran", False) else None,
        "ph_source": ph_source if normalization.get("protonation_step_ran", False) else None,
        "human_review_required": policy == PROTONATION_PH_MODEL,
        "profile_ph_step": "not_run_by_profile" if (direct_profile and policy == PROTONATION_PH_MODEL) else None,
        "input_state": input_state,
        "input_state_error": input_state_error,
        "prepared_state": prepared_state,
        "prepared_state_error": prepared_state_error,
        "pdbqt_state": pdbqt_state,
        "pdbqt_state_check": pdbqt_check,
        "pdbqt_canonical_match": canonical_match,
        "coordinates_source": normalization.get("coordinates_source"),
        "hydrogens_added": normalization.get("hydrogens_added", normalization.get("hydrogens_added_after_ph")),
        "explicit_state": (
            {"smiles": explicit_smiles, "net_charge": int(state_entry["net_charge"]), "stereo_compared": False}
            if state_entry is not None else None
        ),
    }
    summary = {
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_file": str(source),
        "output_file": str(destination),
        "preparation_method": prep_method,
        "requested_backend": selected_profile,
        "requested_profile": selected_profile,
        "effective_profile": effective_profile,
        "selected_engines": list(selected_engines),
        "protonation_ph": protonation_ph,
        "protonation_policy": policy,
        "protonation_applied": protonation_applied,
        "protonation": protonation_block,
        "normalization": normalization,
        "compatibility": compatibility.to_dict(),
    }
    if meeko_error:
        summary["meeko_error"] = meeko_error
    if obabel_error:
        summary["obabel_error"] = obabel_error
    if autodocktools_error:
        summary["autodocktools_error"] = autodocktools_error
    if autodocktools_summary:
        summary["autodocktools"] = autodocktools_summary
    report_dir = str(os.environ.get("PDBWIZARD_LIGAND_PREP_REPORT_DIR", "") or "").strip()
    if report_dir:
        report_path = write_ligand_preparation_step_report(summary, report_dir=Path(report_dir))
        if report_path:
            summary["step_report_file"] = str(report_path)
    validation = validate_ligand_preparation_output_contract(summary)
    summary["output_contract_validation"] = validation
    return summary


def _parse_float(raw: object) -> Optional[float]:
    try:
        return float(raw)
    except Exception:
        return None


def _count_line_token(lines: list[str], token: str) -> int:
    return sum(1 for line in lines if line.strip() == token)


def _validate_prepared_pdbqt_contract(output_file: Path) -> Dict[str, object]:
    errors: list[str] = []
    warnings: list[str] = []
    if not output_file.exists():
        errors.append("prepared output file is missing")
        return {"is_valid": False, "errors": errors, "warnings": warnings}
    text = output_file.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()
    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    root_count = _count_line_token(lines, "ROOT")
    torsdof_count = sum(1 for line in lines if line.startswith("TORSDOF"))
    if not atom_lines:
        errors.append("prepared PDBQT has no ATOM/HETATM records")
    if root_count != 1:
        errors.append(f"prepared PDBQT expected exactly one ROOT block, found {root_count}")
    if torsdof_count != 1:
        errors.append(f"prepared PDBQT expected exactly one TORSDOF record, found {torsdof_count}")
    return {"is_valid": not errors, "errors": errors, "warnings": warnings}


def validate_ligand_preparation_output_contract(summary: Dict[str, object]) -> Dict[str, object]:
    """
    Validate mode-specific output contracts for one prepared ligand summary payload.
    """
    requested_profile = normalize_ligand_preparation_profile(str(summary.get("requested_profile") or summary.get("requested_backend") or "engine_aware_full"))
    selected_engines = [str(engine).strip().lower() for engine in (summary.get("selected_engines") or []) if str(engine).strip()]
    compatibility = validate_ligand_preparation_profile(requested_profile, selected_engines)
    effective_profile = str(summary.get("effective_profile") or compatibility.effective_profile or "")
    output_file = Path(str(summary.get("output_file") or "")).expanduser()
    errors: list[str] = []
    warnings: list[str] = []

    if not compatibility.is_valid:
        errors.extend(compatibility.errors)
    warnings.extend(compatibility.warnings)

    expected_methods = {
        "meeko_only": {"meeko_direct"},
        "autodocktools_only": {"autodocktools_direct"},
        "openbabel_only": {"openbabel_pdbqt"},
        "openbabel_meeko": {"openbabel_then_meeko"},
        "openbabel_autodocktools": {"openbabel_then_autodocktools"},
        "openbabel_meeko_autodock": {"openbabel_then_meeko", "openbabel_then_meeko_fallback_autodocktools"},
    }
    preparation_method = str(summary.get("preparation_method") or "")
    allowed_methods = expected_methods.get(effective_profile)
    if allowed_methods and preparation_method not in allowed_methods:
        errors.append(
            f"preparation_method='{preparation_method}' does not match effective_profile='{effective_profile}' "
            f"(expected one of: {', '.join(sorted(allowed_methods))})"
        )

    ph_value = _parse_float(summary.get("protonation_ph"))
    ph_ok, normalized_ph, ph_error = validate_preparation_ph(ph_value if ph_value is not None else DEFAULT_PREPARATION_PH)
    if not ph_ok:
        errors.append(ph_error)
    elif ph_value is not None and abs(ph_value - normalized_ph) > 1e-6:
        warnings.append("protonation_ph was normalized to a valid range")

    normalization = summary.get("normalization") or {}
    if effective_profile.startswith("openbabel_") or effective_profile == "openbabel_only":
        if not isinstance(normalization, dict) or not normalization:
            errors.append("openbabel-based profiles require normalization details in summary")
        else:
            backend = str(normalization.get("normalization_backend") or "").strip().lower()
            if backend not in {"openbabel", "rdkit", "rdkit_as_input", "rdkit_explicit_state"}:
                errors.append(f"unexpected normalization backend: {backend or 'missing'}")
            if backend == "rdkit":
                warnings.append("normalization used RDKit fallback instead of Open Babel")

    # Spec 034 R2b: protonation_applied must come with a pH step and a measured state.
    protonation = summary.get("protonation") if isinstance(summary.get("protonation"), dict) else {}
    if protonation.get("protonation_applied"):
        prepared_state = protonation.get("prepared_state")
        if not protonation.get("ph_model_run"):
            errors.append("protonation_applied is true but the pH model did not run")
        if not isinstance(prepared_state, dict) or not prepared_state.get("microspecies_smiles"):
            errors.append("protonation_applied is true without a measured prepared state")

    output_contract = _validate_prepared_pdbqt_contract(output_file)
    if not output_contract.get("is_valid", False):
        errors.extend([str(issue) for issue in output_contract.get("errors", [])])
    warnings.extend([str(issue) for issue in output_contract.get("warnings", [])])

    return {
        "requested_profile": requested_profile,
        "effective_profile": effective_profile,
        "selected_engines": selected_engines,
        "output_file": str(output_file),
        "protonation_ph": float(ph_value) if ph_value is not None else float(DEFAULT_PREPARATION_PH),
        "ph_bounds": {
            "min": float(MIN_PREPARATION_PH),
            "max": float(MAX_PREPARATION_PH),
            "default": float(DEFAULT_PREPARATION_PH),
        },
        "is_valid": not errors,
        "errors": errors,
        "warnings": warnings,
    }


def write_ligand_preparation_step_report(
    summary: Dict[str, object],
    *,
    report_dir: Optional[Path] = None,
    report_file: Optional[Path] = None,
) -> Optional[Path]:
    """
    Persist a per-ligand preparation step report as JSON.
    """
    target_path: Optional[Path] = None
    if report_file is not None:
        target_path = Path(report_file).expanduser().resolve()
    elif report_dir is not None:
        directory = Path(report_dir).expanduser().resolve()
        input_file = Path(str(summary.get("input_file") or "ligand")).name
        stem = Path(input_file).stem or "ligand"
        safe_stem = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in stem).strip("_") or "ligand"
        target_path = directory / f"{safe_stem}.json"
    if target_path is None:
        return None
    target_path.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(summary)
    payload["output_contract_validation"] = payload.get("output_contract_validation") or validate_ligand_preparation_output_contract(payload)
    target_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return target_path


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize ligands to explicit-H 3D SDF and prepare AutoDock/Vina-family PDBQT output."
    )
    parser.add_argument("--input", required=True, help="Input ligand file (.sdf, .mol, .mol2, .pdb)")
    parser.add_argument("--output", required=True, help="Output PDBQT path")
    parser.add_argument(
        "--backend-profile",
        default="engine_aware_full",
        help=(
            "Ligand preparation profile: "
            "openbabel_only, meeko_only, autodocktools_only, openbabel_meeko, "
            "openbabel_meeko_autodock, openbabel_autodocktools, engine_aware_full"
        ),
    )
    parser.add_argument(
        "--ph",
        type=float,
        default=None,
        help=(
            f"Protonation pH for the Open Babel pH model ({MIN_PREPARATION_PH:.1f}-{MAX_PREPARATION_PH:.1f}). "
            "Required for ph_model and explicit_state; no default is applied."
        ),
    )
    parser.add_argument(
        "--protonation-policy",
        choices=list(PROTONATION_POLICIES),
        default=DEFAULT_PROTONATION_POLICY,
        help="ph_model (default, human review flagged), explicit_state (state map), or as_input (keep input H and charges)",
    )
    parser.add_argument(
        "--protonation-state-map",
        default="",
        help="CSV with columns ligand, smiles and/or net_charge (required for explicit_state)",
    )
    parser.add_argument("--engines", default="", help="Comma-separated selected engines for engine-aware profile decisions")
    parser.add_argument("--summary-file", help="Optional JSON summary output path")
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)
    policy = str(args.protonation_policy)
    if args.ph is not None:
        ph_ok, normalized_ph, ph_error = validate_preparation_ph(args.ph)
        if not ph_ok:
            print(f"❌ Invalid pH value: {ph_error}", file=sys.stderr)
            return 2
        os.environ["PDBWIZARD_LIGAND_PREP_PH"] = str(normalized_ph)
        os.environ["PDBWIZARD_LIGAND_PREP_PH_SOURCE"] = PH_SOURCE_USER_ENTERED
    else:
        needs_ph = policy == PROTONATION_PH_MODEL and not os.environ.get("PDBWIZARD_LIGAND_PREP_PH")
        if policy == PROTONATION_EXPLICIT_STATE and args.protonation_state_map:
            # Only net-charge-only rows run the pH model; SMILES rows do not need a pH.
            try:
                needs_ph = any(
                    not entry["smiles"] for entry in load_protonation_state_map(Path(args.protonation_state_map)).values()
                ) and not os.environ.get("PDBWIZARD_LIGAND_PREP_PH")
            except ProtonationStateError as exc:
                print(json.dumps({"status": "failed", "reason": exc.reason, "details": exc.details}))
                return 2
        if needs_ph:
            print(
                "❌ --ph is required for this protonation policy. No default pH is applied; "
                "7.4 is a common suggestion that must be passed explicitly.",
                file=sys.stderr,
            )
            return 2
    os.environ["PDBWIZARD_LIGAND_PREP_PROTONATION_POLICY"] = policy
    if args.protonation_state_map:
        os.environ["PDBWIZARD_LIGAND_PREP_STATE_MAP"] = str(Path(args.protonation_state_map).expanduser().resolve())
    os.environ["PDBWIZARD_LIGAND_PREP_PROFILE"] = _resolve_profile_alias(args.backend_profile)
    if args.engines:
        os.environ["PDBWIZARD_SELECTED_ENGINES"] = args.engines
    try:
        summary = prepare_ligand_for_vina_family(Path(args.input), Path(args.output))
    except ProtonationStateError as exc:
        print(json.dumps({"status": "failed", "reason": exc.reason, "details": exc.details}))
        return 1
    if args.summary_file:
        write_ligand_preparation_step_report(summary, report_file=Path(args.summary_file))
    else:
        print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
