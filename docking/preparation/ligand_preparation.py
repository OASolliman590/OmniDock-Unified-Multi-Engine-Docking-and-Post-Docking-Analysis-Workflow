from __future__ import annotations

import argparse
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
from .ligand_identity import (
    LigandIdentityError,
    build_identity_ledger,
    ensure_single_component,
    require_finite_3d_coordinates,
    sha256_file,
    validate_mapped_pdbqt,
)


def _resolve_profile_alias(raw: str) -> str:
    normalized = normalize_ligand_preparation_profile(raw)
    validation = validate_ligand_preparation_profile(raw, [])
    if not validation.is_valid:
        raise ValueError("; ".join(validation.errors))
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


def _selected_ph() -> float:
    raw = str(os.environ.get("PDBWIZARD_LIGAND_PREP_PH", str(DEFAULT_PREPARATION_PH)) or "").strip()
    ok, normalized, error = validate_preparation_ph(raw)
    if not ok:
        raise ValueError(f"Invalid ligand preparation pH: {error}")
    return float(normalized)


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
        molecules = list(supplier)
        if len(molecules) != 1 or molecules[0] is None:
            raise ValueError("Each ligand artifact must contain exactly one valid molecule")
        molecule = molecules[0]
    elif suffix == ".mol":
        molecule = Chem.MolFromMolFile(str(input_path), removeHs=False)
    elif suffix == ".mol2":
        molecule = Chem.MolFromMol2File(str(input_path), removeHs=False)
    elif suffix == ".pdb":
        molecule = Chem.MolFromPDBFile(str(input_path), removeHs=False)
    else:
        raise ValueError(f"Unsupported ligand format for RDKit normalization: {input_path.suffix}")
    if molecule is None:
        raise ValueError(f"RDKit could not parse ligand file: {input_path}")
    try:
        ensure_single_component(molecule, context=f"ligand '{input_path.name}'")
    except LigandIdentityError as exc:
        raise ValueError(str(exc)) from exc
    return molecule


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

    source_molecule = _load_rdkit_molecule(input_path)
    if source_molecule is None:
        raise ValueError(f"RDKit could not parse ligand file: {input_path}")

    molecule = Chem.Mol(source_molecule)
    Chem.SanitizeMol(molecule)
    had_3d = _has_3d_coordinates(molecule)
    molecule = Chem.AddHs(molecule, addCoords=had_3d)

    optimized_forcefield = ""
    if not had_3d:
        optimized_forcefield = _embed_and_optimize_3d(molecule, all_chem)

    molecule.SetProp("_Name", molecule.GetProp("_Name") if molecule.HasProp("_Name") else input_path.stem)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    # Keep the staged candidate recognizable by the suffix-dispatch parser.
    # A name ending only in ``.candidate`` is rejected before chemistry can be
    # checked by ``_load_rdkit_molecule``.
    candidate_output = output_sdf.parent / f".{output_sdf.stem}.candidate.sdf"
    try:
        writer = Chem.SDWriter(str(candidate_output))
        try:
            writer.write(molecule)
        finally:
            writer.close()

        normalized_molecule = _load_rdkit_molecule(candidate_output)
        require_finite_3d_coordinates(normalized_molecule, context="RDKit normalized ligand")
        identity_ledger = build_identity_ledger(
            source_molecule,
            normalized_molecule,
            allow_hydrogen_changes=True,
            allow_formal_charge_changes=False,
            preserve_coordinates=True,
        )
        if not identity_ledger["is_valid"]:
            raise ValueError(
                "RDKit normalization changed ligand identity: "
                + "; ".join(str(error) for error in identity_ledger["errors"])
            )
        output_sdf.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(candidate_output, output_sdf)
    finally:
        try:
            candidate_output.unlink()
        except OSError:
            pass

    return {
        "normalization_backend": "rdkit",
        "source_file": str(input_path),
        "normalized_sdf": str(output_sdf),
        "source_sha256": sha256_file(input_path),
        "normalized_sdf_sha256": sha256_file(output_sdf),
        "had_3d_input": had_3d,
        "generated_3d": not had_3d,
        "optimized_forcefield": optimized_forcefield,
        "identity_ledger": identity_ledger,
        "normalized_isomeric_smiles": identity_ledger.get("normalized_isomeric_smiles", ""),
        "coordinate_provenance": "source_3d_preserved" if had_3d else "rdkit_generated_seeded",
    }


def _normalize_with_openbabel(
    input_path: Path,
    output_sdf: Path,
    *,
    protonation_ph: float = 7.4,
    partial_charge_model: str = "gasteiger",
) -> Dict[str, object]:
    ph_ok, normalized_ph, ph_error = validate_preparation_ph(protonation_ph)
    if not ph_ok:
        raise ValueError(f"Invalid ligand preparation pH: {ph_error}")
    protonation_ph = normalized_ph
    if not shutil.which("obabel"):
        raise FileNotFoundError("Open Babel (obabel) is required for ligand normalization fallback")

    with tempfile.TemporaryDirectory(prefix=f"ligand_normalize_{input_path.stem}_") as tmp_dir:
        temp_dir = Path(tmp_dir)
        temp_source = temp_dir / f"{input_path.stem}_source.sdf"
        temp_3d = temp_dir / f"{input_path.stem}_3d.sdf"
        temp_h = temp_dir / f"{input_path.stem}_h.sdf"
        temp_charged = temp_dir / f"{input_path.stem}_charged.sdf"
        temp_min = temp_dir / f"{input_path.stem}_min.sdf"

        molecule = _load_rdkit_molecule(input_path)
        if molecule is None:
            raise ValueError("Ligand chemistry could not be parsed")
        ensure_single_component(molecule, context="source ligand")
        had_3d = _has_3d_coordinates(molecule)
        source_result = _run_command(["obabel", str(input_path), "-O", str(temp_source)])
        if _is_zero_molecule_obabel_result(source_result) or not _file_has_meaningful_content(temp_source):
            raise ValueError(f"OpenBabel could not parse ligand source into SDF: {input_path}")
        if had_3d:
            shutil.copy2(temp_source, temp_3d)
        else:
            gen3d_result = _run_command(["obabel", str(temp_source), "-O", str(temp_3d), "--gen3d"])
            if _is_zero_molecule_obabel_result(gen3d_result) or not _file_has_meaningful_content(temp_3d):
                raise ValueError(f"OpenBabel failed to generate 3D conformer for: {input_path}")
        protonation_applied = False
        try:
            addh_result = _run_command(
                [
                    "obabel",
                    str(temp_3d),
                    "-O",
                    str(temp_h),
                    "-p",
                    f"{protonation_ph:.2f}",
                ]
            )
            protonation_applied = True
        except subprocess.CalledProcessError as exc:
            raise RuntimeError("Required pH-dependent ligand protonation failed; simple hydrogen addition is not an equivalent fallback") from exc
        if _is_zero_molecule_obabel_result(addh_result) or not _file_has_meaningful_content(temp_h):
            raise ValueError(f"OpenBabel failed to add hydrogens for: {input_path}")

        final_source = temp_h
        partial_charge_applied = False
        if partial_charge_model:
            try:
                charge_result = _run_command(
                    [
                        "obabel",
                        str(temp_h),
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
        for forcefield in (() if had_3d else ("MMFF94", "UFF")):
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

        candidate_output = temp_dir / f"{input_path.stem}_validated.sdf"
        shutil.copy2(final_source, candidate_output)
        if not _file_has_meaningful_content(candidate_output):
            raise ValueError(f"OpenBabel normalization produced an empty SDF for: {input_path}")

        normalized_molecule = _load_rdkit_molecule(candidate_output)
        try:
            require_finite_3d_coordinates(normalized_molecule, context="OpenBabel normalized ligand")
        except LigandIdentityError as exc:
            raise ValueError(str(exc)) from exc
        identity_ledger = build_identity_ledger(
            molecule,
            normalized_molecule,
            allow_hydrogen_changes=True,
            allow_formal_charge_changes=protonation_applied,
            preserve_coordinates=True,
        )
        if not identity_ledger["is_valid"]:
            raise ValueError(
                "OpenBabel normalization changed ligand identity: "
                + "; ".join(str(error) for error in identity_ledger["errors"])
            )
        output_sdf.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(candidate_output, output_sdf)

    return {
        "normalization_backend": "openbabel",
        "source_file": str(input_path),
        "normalized_sdf": str(output_sdf),
        "source_sha256": sha256_file(input_path),
        "normalized_sdf_sha256": sha256_file(output_sdf),
        "had_3d_input": had_3d,
        "generated_3d": not had_3d,
        "protonation_ph": protonation_ph,
        "protonation_applied": protonation_applied,
        "partial_charge_model": partial_charge_model,
        "partial_charge_applied": partial_charge_applied,
        "charge_state_status": (
            "normalization_charge_applied" if partial_charge_applied else "normalization_charge_unevaluated"
        ),
        "optimized_forcefield": optimized_forcefield,
        "identity_ledger": identity_ledger,
        "normalized_isomeric_smiles": identity_ledger.get("normalized_isomeric_smiles", ""),
        "coordinate_provenance": (
            "source_3d_preserved" if had_3d else "openbabel_generated_nondeterministic"
        ),
        "normalization_provenance": {
            "procedure": "openbabel_pH_protonation",
            "pH": float(protonation_ph),
            "allowed_changes": ["hydrogens", "formal_charge"],
            "heavy_atom_graph_policy": "exact",
            "source_coordinates": "preserved_and_mapped" if had_3d else "not_present",
        },
    }


def normalize_ligand_to_sdf(
    input_path: Path,
    output_sdf: Path,
    *,
    protonation_ph: float = 7.4,
) -> Dict[str, object]:
    source = Path(input_path).expanduser().resolve()
    destination = Path(output_sdf).expanduser().resolve()
    if source.suffix.lower() == ".pdb":
        raise ValueError("PDB coordinates do not establish ligand chemistry; provide an authoritative SDF with bond orders, charges and stereochemistry")
    ph_ok, normalized_ph, ph_error = validate_preparation_ph(protonation_ph)
    if not ph_ok:
        raise ValueError(f"Invalid ligand preparation pH: {ph_error}")
    # Parse and reject disconnected input before checking or invoking any
    # external normalizer.  This keeps every preparation profile on the same
    # ordinary-ligand policy.
    source_molecule = _load_rdkit_molecule(source)
    ensure_single_component(source_molecule, context="source ligand")
    # Failure is explicit: RDKit AddHs cannot satisfy a requested pH treatment.
    # Keep the public destination staged as well as the backend's own candidate
    # so a malformed 2D/nonfinite result cannot overwrite a prior artifact.
    with tempfile.TemporaryDirectory(prefix=f"ligand_normalized_{source.stem}_") as staging_dir:
        staged_output = Path(staging_dir) / f"{destination.stem}.candidate.sdf"
        result = _normalize_with_openbabel(source, staged_output, protonation_ph=normalized_ph)
        normalized_molecule = _load_rdkit_molecule(staged_output)
        try:
            require_finite_3d_coordinates(normalized_molecule, context="normalized ligand")
        except LigandIdentityError as exc:
            raise ValueError(str(exc)) from exc
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(staged_output, destination)
    result = dict(result)
    result["normalized_sdf"] = str(destination)
    result["normalized_sdf_sha256"] = sha256_file(destination)
    return result


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
    if source.suffix.lower() == ".pdb":
        raise ValueError("Prepare ligands from authoritative SDF chemistry, not inferred PDB connectivity")
    from .asset_identity import sidecar
    chemistry_status = str(sidecar(source).get("chemistry_status", ""))
    if chemistry_status and chemistry_status.startswith("not_evaluated"):
        raise ValueError("Ligand chemistry provenance is not evaluated; retrieve or provide a verified SDF")
    source_molecule = _load_rdkit_molecule(source)
    if source_molecule is None:
        raise ValueError("Ligand chemistry is not evaluable")
    ensure_single_component(source_molecule, context="source ligand")
    Chem, _ = _load_rdkit()
    source_smiles = Chem.MolToSmiles(Chem.RemoveHs(source_molecule), isomericSmiles=True)

    selected_profile = _selected_profile()
    selected_engines = _selected_engines()
    compatibility = validate_ligand_preparation_profile(selected_profile, list(selected_engines))
    if not compatibility.is_valid:
        raise RuntimeError("; ".join(compatibility.errors))
    effective_profile = _effective_profile(selected_profile, selected_engines)
    protonation_ph = _selected_ph()
    input_sha256 = sha256_file(source)

    with tempfile.TemporaryDirectory(prefix=f"ligand_prepare_{source.stem}_") as tmp_dir:
        normalized_sdf = Path(tmp_dir) / f"{source.stem}_normalized.sdf"
        staged_output = Path(tmp_dir) / f"{destination.name}.candidate.pdbqt"
        meeko_error = ""
        obabel_error = ""
        autodocktools_error = ""
        prep_method = ""
        autodocktools_summary: Optional[Dict[str, object]] = None
        normalization: Dict[str, object] = {}
        prepared_input_molecule = source_molecule

        if effective_profile == "meeko_only":
            try:
                require_finite_3d_coordinates(source_molecule, context="source ligand")
            except LigandIdentityError as exc:
                raise ValueError(str(exc)) from exc
            try:
                _prepare_with_meeko(source, staged_output)
                prep_method = "meeko_direct"
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                meeko_error = str(exc)
                raise RuntimeError("Ligand preparation failed in Meeko-only mode.") from exc
        elif effective_profile == "autodocktools_only":
            try:
                require_finite_3d_coordinates(source_molecule, context="source ligand")
            except LigandIdentityError as exc:
                raise ValueError(str(exc)) from exc
            try:
                autodocktools_summary = _prepare_with_autodocktools(source, staged_output)
                prep_method = "autodocktools_direct"
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                autodocktools_error = str(exc)
                raise RuntimeError("Ligand preparation failed in AutoDockTools-only mode.") from exc
        else:
            normalization = normalize_ligand_to_sdf(source, normalized_sdf, protonation_ph=protonation_ph)
            durable_normalized_sdf = destination.with_suffix(".normalized.sdf")
            normalization["source_sha256"] = input_sha256
            prepared_input_molecule = _load_rdkit_molecule(normalized_sdf)
            try:
                require_finite_3d_coordinates(
                    prepared_input_molecule,
                    context="prepared normalized ligand",
                )
            except LigandIdentityError as exc:
                raise ValueError(str(exc)) from exc
            identity_ledger = build_identity_ledger(
                source_molecule,
                prepared_input_molecule,
                allow_hydrogen_changes=True,
                allow_formal_charge_changes=bool(normalization.get("protonation_applied")),
                preserve_coordinates=True,
            )
            if not identity_ledger.get("is_valid"):
                raise ValueError(
                    "Normalized ligand failed the source identity contract: "
                    + "; ".join(str(error) for error in identity_ledger.get("errors", []))
                )
            normalization["identity_ledger"] = identity_ledger
            # Keep the durable normalized SDF staged until backend output has
            # also passed validation.  A failed backend must leave any prior
            # user artifact at the requested normalized path untouched.
            normalization["normalized_sdf"] = str(durable_normalized_sdf)
            normalization["normalized_sdf_sha256"] = sha256_file(normalized_sdf)
            if effective_profile == "openbabel_only":
                try:
                    _prepare_with_openbabel_pdbqt(normalized_sdf, staged_output)
                    prep_method = "openbabel_pdbqt"
                except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                    obabel_error = str(exc)
                    raise RuntimeError("Ligand preparation failed in OpenBabel-only mode.") from exc
            elif effective_profile == "openbabel_meeko":
                try:
                    _prepare_with_meeko(normalized_sdf, staged_output)
                    prep_method = "openbabel_then_meeko"
                except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                    meeko_error = str(exc)
                    raise RuntimeError("Ligand preparation failed in OpenBabel->Meeko mode.") from exc
            elif effective_profile == "openbabel_autodocktools":
                try:
                    autodocktools_summary = _prepare_with_autodocktools(normalized_sdf, staged_output)
                    prep_method = "openbabel_then_autodocktools"
                except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                    autodocktools_error = str(exc)
                    raise RuntimeError("Ligand preparation failed in OpenBabel->AutoDockTools mode.") from exc
            else:
                # openbabel_meeko_autodock
                try:
                    _prepare_with_meeko(normalized_sdf, staged_output)
                    prep_method = "openbabel_then_meeko"
                except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                    meeko_error = str(exc)
                    try:
                        autodocktools_summary = _prepare_with_autodocktools(normalized_sdf, staged_output)
                        prep_method = "openbabel_then_meeko_fallback_autodocktools"
                    except (subprocess.CalledProcessError, FileNotFoundError) as adt_exc:
                        autodocktools_error = str(adt_exc)
                        raise RuntimeError(
                            "Ligand preparation failed in OpenBabel->Meeko->AutoDockTools mode."
                        ) from exc

        staged_contract = _validate_prepared_pdbqt_contract(staged_output)
        if not staged_contract.get("is_valid", False):
            raise ValueError(
                "Prepared ligand failed its chemistry/engine contract before publication: "
                + "; ".join(str(error) for error in staged_contract.get("errors", []))
            )
        output_chemistry_validation = validate_mapped_pdbqt(prepared_input_molecule, staged_output)
        if output_chemistry_validation.get("status") == "invalid":
            raise ValueError(
                "Prepared PDBQT contains malformed or conflicting authoritative chemistry mapping: "
                + "; ".join(str(error) for error in output_chemistry_validation.get("errors", []))
            )
        if "meeko" in prep_method and "fallback_autodocktools" not in prep_method:
            if output_chemistry_validation.get("status") != "exact":
                raise ValueError(
                    "Meeko output could not be certified against the prepared input graph: "
                    + "; ".join(str(error) for error in output_chemistry_validation.get("errors", []))
                )
        elif output_chemistry_validation.get("status") != "exact":
            # Deliberate compatibility for legacy ADT/OpenBabel outputs: the
            # syntax/engine contract is still enforced, but chemistry is
            # limited when the backend does not emit an authoritative map.
            output_chemistry_validation.setdefault(
                "warnings", []
            ).append(
                "The selected OpenBabel/AutoDockTools output has no authoritative atom mapping; "
                "exact output chemistry was not certified."
            )
        # All graph, coordinate, syntax and output-mapping checks have passed;
        # publish the two prepared artifacts only now.
        if normalization:
            durable_normalized_sdf.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(normalized_sdf, durable_normalized_sdf)
        shutil.copy2(staged_output, destination)

    summary = {
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_file": str(source),
        "output_file": str(destination),
        "input_sha256": input_sha256,
        "output_sha256": sha256_file(destination),
        "preparation_method": prep_method,
        "requested_backend": selected_profile,
        "requested_profile": selected_profile,
        "effective_profile": effective_profile,
        "selected_engines": list(selected_engines),
        "protonation_ph": protonation_ph,
        "requested_ph": protonation_ph,
        "protonation_status": "applied" if normalization.get("protonation_applied") else "input_state_not_ph_titrated",
        "source_isomeric_smiles": source_smiles,
        "normalization": normalization,
        "normalization_provenance": normalization.get("normalization_provenance", {}),
        "prepared_input_file": normalization.get("normalized_sdf", str(source)),
        "prepared_input_sha256": (
            normalization.get("normalized_sdf_sha256", input_sha256)
        ),
        "output_chemistry_validation": output_chemistry_validation,
        "chemistry_validation_scope": output_chemistry_validation.get("scope", ""),
        "charge_state_status": (
            "backend_assigned_meeko"
            if "meeko" in prep_method and "fallback_autodocktools" not in prep_method
            else "backend_assigned_autodocktools"
            if "autodocktools" in prep_method
            else "backend_assigned_openbabel_gasteiger"
        ),
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
    validation = validate_ligand_preparation_output_contract(summary)
    summary["output_contract_validation"] = validation
    if not validation["is_valid"]:
        raise ValueError("Prepared ligand failed its chemistry/engine contract: " + "; ".join(validation["errors"]))
    report_dir = str(os.environ.get("PDBWIZARD_LIGAND_PREP_REPORT_DIR", "") or "").strip()
    if report_dir:
        report_path = write_ligand_preparation_step_report(summary, report_dir=Path(report_dir))
        if report_path:
            summary["step_report_file"] = str(report_path)
    from .asset_identity import REFERENCE_FIELDS, sidecar
    source_metadata = sidecar(source)
    for key in REFERENCE_FIELDS:
        summary[key] = source_metadata.get(key, "")
    destination.with_suffix(destination.suffix + ".preparation.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
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
    from .ligand_quality import validate_prepared_ligand_pdbqt
    errors.extend(issue.details for issue in validate_prepared_ligand_pdbqt(output_file, _selected_engines()))
    return {"is_valid": not errors, "errors": errors, "warnings": warnings}


def validate_ligand_preparation_output_contract(summary: Dict[str, object]) -> Dict[str, object]:
    """
    Validate mode-specific output contracts for one prepared ligand summary payload.
    """
    raw_profile = str(
        summary.get("requested_profile")
        or summary.get("requested_backend")
        or "engine_aware_full"
    )
    selected_engines = [str(engine).strip().lower() for engine in (summary.get("selected_engines") or []) if str(engine).strip()]
    compatibility = validate_ligand_preparation_profile(raw_profile, selected_engines)
    requested_profile = compatibility.requested_profile
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

    raw_ph = summary.get("protonation_ph", DEFAULT_PREPARATION_PH)
    if raw_ph is None or raw_ph == "":
        raw_ph = DEFAULT_PREPARATION_PH
    ph_value = _parse_float(raw_ph)
    ph_ok, normalized_ph, ph_error = validate_preparation_ph(raw_ph if ph_value is not None else raw_ph)
    if not ph_ok:
        errors.append(ph_error)
    elif ph_value is not None and abs(ph_value - normalized_ph) > 1e-6:
        warnings.append("protonation_ph was normalized to a valid range")

    normalization = summary.get("normalization") or {}
    uses_normalization = effective_profile.startswith("openbabel_") or effective_profile == "openbabel_only"
    if uses_normalization:
        if not isinstance(normalization, dict) or not normalization:
            errors.append("openbabel-based profiles require normalization details in summary")
        else:
            backend = str(normalization.get("normalization_backend") or "").strip().lower()
            if backend not in {"openbabel", "rdkit"}:
                errors.append(f"unexpected normalization backend: {backend or 'missing'}")
            if backend == "rdkit" or not normalization.get("protonation_applied", False):
                errors.append("Requested pH-dependent normalization was not applied")
            normalized_file_raw = str(normalization.get("normalized_sdf") or "").strip()
            normalized_file = Path(normalized_file_raw).expanduser() if normalized_file_raw else None
            normalized_hash = str(normalization.get("normalized_sdf_sha256") or "").strip().lower()
            if not normalized_file:
                errors.append("normalization must retain a durable normalized SDF path")
            elif not normalized_file.exists():
                errors.append("durable normalized SDF is missing")
            elif not normalized_hash:
                errors.append("normalization must retain the normalized SDF SHA-256 hash")
            elif sha256_file(normalized_file) != normalized_hash:
                errors.append("normalized SDF SHA-256 does not match the durable artifact")
            identity_ledger = normalization.get("identity_ledger")
            if not isinstance(identity_ledger, dict) or not identity_ledger.get("is_valid"):
                errors.append("normalization identity ledger is missing or failed")
            provenance = normalization.get("normalization_provenance")
            if not isinstance(provenance, dict) or not provenance.get("procedure"):
                errors.append("normalization provenance and allowed-change policy are required")

    input_file = Path(str(summary.get("input_file") or "")).expanduser()
    input_hash = str(summary.get("input_sha256") or "").strip().lower()
    if not input_hash:
        errors.append("input SHA-256 provenance is missing")
    elif input_file.exists() and sha256_file(input_file) != input_hash:
        errors.append("input SHA-256 does not match the source artifact")
    output_hash = str(summary.get("output_sha256") or "").strip().lower()
    if not output_hash:
        errors.append("output SHA-256 provenance is missing")
    elif output_file.exists() and sha256_file(output_file) != output_hash:
        errors.append("output SHA-256 does not match the prepared artifact")

    output_chemistry = summary.get("output_chemistry_validation")
    meeko_output = "meeko" in preparation_method and "fallback_autodocktools" not in preparation_method
    if not isinstance(output_chemistry, dict):
        errors.append("output chemistry validation scope is missing")
    elif meeko_output:
        if output_chemistry.get("status") != "exact" or not output_chemistry.get("is_valid"):
            errors.append("Meeko output requires exact authoritative chemistry/mapping validation")
    elif output_chemistry.get("status") not in {"limited", "exact"}:
        errors.append("OpenBabel/AutoDockTools output must declare limited chemistry validation scope")
    elif output_chemistry.get("status") == "exact" and not output_chemistry.get("is_valid"):
        errors.append("Exact output chemistry validation must declare is_valid=true")
    elif output_chemistry.get("status") == "limited":
        warnings.extend(str(item) for item in (output_chemistry.get("warnings") or []))

    if not str(summary.get("charge_state_status") or "").strip():
        errors.append("charge-state provenance is missing")

    output_contract = _validate_prepared_pdbqt_contract(output_file)
    if not output_contract.get("is_valid", False):
        errors.extend([str(issue) for issue in output_contract.get("errors", [])])
    warnings.extend([str(issue) for issue in output_contract.get("warnings", [])])

    return {
        "requested_profile": requested_profile,
        "effective_profile": effective_profile,
        "selected_engines": selected_engines,
        "output_file": str(output_file),
        "input_sha256": input_hash,
        "output_sha256": output_hash,
        "protonation_ph": float(ph_value) if ph_value is not None else float(DEFAULT_PREPARATION_PH),
        "chemistry_validation_scope": str(summary.get("chemistry_validation_scope") or ""),
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
    parser.add_argument("--input", required=True, help="Input ligand with authoritative chemistry (.sdf, .mol, .mol2)")
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
        default=DEFAULT_PREPARATION_PH,
        help=f"Protonation pH for Open Babel normalization stages ({MIN_PREPARATION_PH:.1f}-{MAX_PREPARATION_PH:.1f})",
    )
    parser.add_argument("--engines", default="", help="Comma-separated selected engines for engine-aware profile decisions")
    parser.add_argument("--summary-file", help="Optional JSON summary output path")
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)
    ph_ok, normalized_ph, ph_error = validate_preparation_ph(args.ph)
    if not ph_ok:
        print(f"❌ Invalid pH value: {ph_error}", file=sys.stderr)
        return 2
    profile_validation = validate_ligand_preparation_profile(
        args.backend_profile,
        normalize_engine_names(str(args.engines or "").split(",")),
    )
    if not profile_validation.is_valid:
        print(f"❌ Invalid ligand preparation profile: {'; '.join(profile_validation.errors)}", file=sys.stderr)
        return 2
    os.environ["PDBWIZARD_LIGAND_PREP_PROFILE"] = profile_validation.requested_profile
    os.environ["PDBWIZARD_LIGAND_PREP_PH"] = str(normalized_ph)
    if args.engines:
        os.environ["PDBWIZARD_SELECTED_ENGINES"] = args.engines
    summary = prepare_ligand_for_vina_family(Path(args.input), Path(args.output))
    if args.summary_file:
        write_ligand_preparation_step_report(summary, report_file=Path(args.summary_file))
    else:
        print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
