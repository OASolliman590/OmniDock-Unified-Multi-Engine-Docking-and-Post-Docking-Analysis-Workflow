from __future__ import annotations

import json
import shlex
import shutil
import subprocess
import tempfile
import re
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from ..models import PairlistRow
from .ligand_preparation import prepare_ligand_for_vina_family
from ..project_layout import detect_layout_profile, ensure_project_layout, load_manifest

VALID_AUTODOCK_TYPES = {
    "A",
    "Br",
    "C",
    "Ca",
    "Cl",
    "Cu",
    "F",
    "Fe",
    "H",
    "HD",
    "HS",
    "I",
    "Mg",
    "Mn",
    "N",
    "NA",
    "NS",
    "OA",
    "OS",
    "P",
    "Pd",
    "S",
    "SA",
    "Si",
    "Zn",
}

LIGAND_ADMET_DEFAULTS: Dict[str, object] = {
    "enable_admet_filters": True,
    "max_lipinski_violations": 1,
    "max_molecular_weight": 650.0,
    "max_logp": 6.0,
    "max_tpsa": 180.0,
    "max_rotatable_bonds": 15,
    "max_formal_charge_abs": 2,
    "min_heavy_atom_count": 6,
    "block_pains": True,
    "block_brenk": True,
    "block_reactive": True,
}

REACTIVE_SMARTS_PATTERNS = [
    ("nitro_group", "[N+](=O)[O-]"),
    ("thiol", "[SH]"),
    ("isothiocyanate", "N=C=S"),
    ("acid_halide", "[CX3](=O)[Cl,Br,I]"),
    ("epoxide", "C1OC1"),
]


def _load_rdkit_modules() -> Dict[str, object]:
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

        return {
            "Chem": Chem,
            "Crippen": Crippen,
            "Descriptors": Descriptors,
            "Lipinski": Lipinski,
            "rdMolDescriptors": rdMolDescriptors,
            "FilterCatalog": FilterCatalog,
            "FilterCatalogParams": FilterCatalogParams,
        }
    except Exception:
        return {}


_RDKIT_MODULES = _load_rdkit_modules()
_PAINS_CATALOG = None
_BRENK_CATALOG = None


@dataclass
class LigandValidationIssue:
    ligand: str
    prepared_file: str
    reason: str
    details: str = ""
    raw_pdb_file: str = ""
    root_count: int = 0
    torsdof_count: int = 0
    affected_pairs: int = 0
    affected_tags: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


def _normalize_ligand_admet_thresholds(overrides: Optional[Dict[str, object]] = None) -> Dict[str, object]:
    merged = dict(LIGAND_ADMET_DEFAULTS)
    for key, value in (overrides or {}).items():
        if key in merged and value is not None:
            merged[key] = value

    normalized: Dict[str, object] = {}
    normalized["enable_admet_filters"] = bool(merged["enable_admet_filters"])
    normalized["max_lipinski_violations"] = max(0, int(merged["max_lipinski_violations"]))
    normalized["max_molecular_weight"] = float(merged["max_molecular_weight"])
    normalized["max_logp"] = float(merged["max_logp"])
    normalized["max_tpsa"] = float(merged["max_tpsa"])
    normalized["max_rotatable_bonds"] = max(0, int(merged["max_rotatable_bonds"]))
    normalized["max_formal_charge_abs"] = max(0, int(merged["max_formal_charge_abs"]))
    normalized["min_heavy_atom_count"] = max(0, int(merged["min_heavy_atom_count"]))
    normalized["block_pains"] = bool(merged["block_pains"])
    normalized["block_brenk"] = bool(merged["block_brenk"])
    normalized["block_reactive"] = bool(merged["block_reactive"])
    return normalized


def _get_filter_catalog(kind: str):
    global _PAINS_CATALOG, _BRENK_CATALOG
    if not _RDKIT_MODULES:
        return None
    FilterCatalog = _RDKIT_MODULES["FilterCatalog"]
    FilterCatalogParams = _RDKIT_MODULES["FilterCatalogParams"]

    if kind == "pains":
        if _PAINS_CATALOG is None:
            params = FilterCatalogParams()
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
            _PAINS_CATALOG = FilterCatalog(params)
        return _PAINS_CATALOG

    if kind == "brenk":
        if _BRENK_CATALOG is None:
            params = FilterCatalogParams()
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
            _BRENK_CATALOG = FilterCatalog(params)
        return _BRENK_CATALOG
    return None


def _count_line_token(lines: Sequence[str], token: str) -> int:
    return sum(1 for line in lines if line.strip() == token)


def _extract_atom_type(line: str) -> str:
    parts = line.split()
    return parts[-1].strip() if parts else ""


def _has_invalid_autodock_type(atom_type: str) -> bool:
    if not atom_type:
        return True
    if atom_type in VALID_AUTODOCK_TYPES:
        return False
    if any(char.isdigit() for char in atom_type):
        return True
    return atom_type not in VALID_AUTODOCK_TYPES


def _replace_atom_type(line: str, atom_type: str) -> str:
    return re.sub(r"(\s+)(\S+)\s*$", lambda match: f"{match.group(1)}{atom_type}", line)


def _infer_fixed_atom_type(atom_name: str, atom_type: str) -> Optional[str]:
    raw = str(atom_type).strip()
    if raw == "G0":
        return None
    if raw == "CG0":
        return "C"

    letters = "".join(ch for ch in raw if ch.isalpha())
    if letters in VALID_AUTODOCK_TYPES:
        return letters

    atom_letters = "".join(ch for ch in str(atom_name).strip() if ch.isalpha()).upper()
    if atom_letters.startswith("CL"):
        return "Cl"
    if atom_letters.startswith("BR"):
        return "Br"
    if atom_letters.startswith("ZN"):
        return "Zn"
    if atom_letters.startswith("FE"):
        return "Fe"
    if atom_letters.startswith("MG"):
        return "Mg"
    if atom_letters.startswith("MN"):
        return "Mn"
    if atom_letters.startswith("CA"):
        return "Ca"
    if atom_letters.startswith("CU"):
        return "Cu"
    if atom_letters.startswith("PD"):
        return "Pd"
    if atom_letters.startswith("OA"):
        return "OA"
    if atom_letters.startswith("OS"):
        return "OS"
    if atom_letters.startswith("NA"):
        return "NA"
    if atom_letters.startswith("NS"):
        return "NS"
    if atom_letters.startswith("SA"):
        return "SA"
    if atom_letters.startswith("SI"):
        return "Si"
    if atom_letters.startswith("HD"):
        return "HD"
    if atom_letters.startswith("HS"):
        return "HS"
    if atom_letters:
        first = atom_letters[0]
        if first in {"C", "N", "O", "S", "P", "H", "F", "I", "A"}:
            return "OA" if first == "O" else first
    return None


def sanitize_prepared_ligand_pdbqt(source_file: Path, destination_file: Path) -> Dict[str, object]:
    source = Path(source_file).expanduser().resolve()
    destination = Path(destination_file).expanduser()
    lines = source.read_text(encoding="utf-8", errors="ignore").splitlines()
    cleaned_lines: List[str] = []
    dropped_atoms = 0
    replaced_atoms = 0
    invalid_before: List[str] = []

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            cleaned_lines.append(line)
            continue
        atom_type = _extract_atom_type(line)
        if not _has_invalid_autodock_type(atom_type):
            cleaned_lines.append(line)
            continue

        invalid_before.append(atom_type)
        atom_name = line[12:16].strip()
        fixed_type = _infer_fixed_atom_type(atom_name, atom_type)
        if fixed_type is None:
            dropped_atoms += 1
            continue
        replaced_atoms += 1
        cleaned_lines.append(_replace_atom_type(line, fixed_type))

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(cleaned_lines) + "\n", encoding="utf-8")
    return {
        "source_file": str(source),
        "destination_file": str(destination),
        "replaced_atoms": replaced_atoms,
        "dropped_atoms": dropped_atoms,
        "invalid_atom_types_before": sorted(set(invalid_before)),
    }


def _parse_float(segment: str, default: float) -> float:
    try:
        return float(segment.strip())
    except ValueError:
        return default


def _infer_element(atom_name: str, element_field: str) -> str:
    element = "".join(char for char in str(element_field).strip() if char.isalpha())
    if element:
        return element.title()
    letters = "".join(char for char in str(atom_name).strip() if char.isalpha())
    if not letters:
        return "X"
    if len(letters) >= 2 and letters[:2].title() in {"Cl", "Br", "Na", "Ca", "Mg", "Zn", "Fe", "Mn", "Cu", "Pd"}:
        return letters[:2].title()
    return letters[0].upper()


def _safe_atom_name(element: str, count: int) -> str:
    prefix = "".join(char for char in element.upper() if char.isalpha())[:2] or "X"
    return f"{prefix}{count}"[:4]


def _format_atom_name(name: str, element: str) -> str:
    if len(element.strip()) == 1:
        return f"{name:>4s}"[-4:]
    return f"{name:<4s}"[:4]


def _format_pdb_atom_line(
    record: str,
    serial: int,
    atom_name: str,
    resname: str,
    chain_id: str,
    resseq: str,
    icode: str,
    x: float,
    y: float,
    z: float,
    occupancy: float,
    bfactor: float,
    element: str,
) -> str:
    atom_field = _format_atom_name(atom_name, element)
    return (
        f"{record:<6}{serial:>5} {atom_field}{' ':1}"
        f"{resname:>3} {chain_id[:1] or 'A'}{resseq:>4}{icode[:1] or ' '}"
        f"   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{bfactor:>6.2f}"
        f"          {element:>2}"
    )


def sanitize_ligand_pdb_file(source_file: Path, destination_file: Path, preferred_altloc: str = "A") -> Dict[str, object]:
    source = Path(source_file).expanduser().resolve()
    destination = Path(destination_file).expanduser()
    lines = source.read_text(encoding="utf-8", errors="ignore").splitlines()

    grouped: Dict[tuple[str, str, str, str, str, str], List[Dict[str, object]]] = defaultdict(list)
    order: List[tuple[str, str, str, str, str, str]] = []

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_name = line[12:16]
        key = (
            line[0:6].strip() or "HETATM",
            atom_name.strip(),
            line[17:20].strip() or "LIG",
            line[21].strip() or "A",
            line[22:26].strip() or "1",
            line[26].strip(),
        )
        if key not in grouped:
            order.append(key)
        grouped[key].append(
            {
                "record": line[0:6].strip() or "HETATM",
                "atom_name": atom_name.strip() or "X",
                "altloc": line[16].strip(),
                "resname": line[17:20].strip() or "LIG",
                "chain_id": line[21].strip() or "A",
                "resseq": line[22:26].strip() or "1",
                "icode": line[26].strip(),
                "x": _parse_float(line[30:38], 0.0),
                "y": _parse_float(line[38:46], 0.0),
                "z": _parse_float(line[46:54], 0.0),
                "occupancy": _parse_float(line[54:60], 1.0),
                "bfactor": _parse_float(line[60:66], 0.0),
                "element": _infer_element(atom_name, line[76:78]),
            }
        )

    selected_atoms: List[Dict[str, object]] = []
    removed_altlocs = 0
    for key in order:
        candidates = grouped[key]
        if len(candidates) == 1:
            selected_atoms.append(candidates[0])
            continue
        blank = [candidate for candidate in candidates if not candidate["altloc"]]
        if blank:
            selected = blank[0]
        else:
            preferred = [candidate for candidate in candidates if str(candidate["altloc"]).upper() == preferred_altloc.upper()]
            if preferred:
                selected = preferred[0]
            else:
                selected = max(candidates, key=lambda candidate: float(candidate["occupancy"]))
        removed_altlocs += len(candidates) - 1
        selected_atoms.append(selected)

    element_counts: Dict[str, int] = defaultdict(int)
    output_lines: List[str] = []
    for serial, atom in enumerate(selected_atoms, start=1):
        element = str(atom["element"] or "X")
        element_counts[element] += 1
        output_lines.append(
            _format_pdb_atom_line(
                record=str(atom["record"]),
                serial=serial,
                atom_name=_safe_atom_name(element, element_counts[element]),
                resname=str(atom["resname"]),
                chain_id=str(atom["chain_id"]),
                resseq=str(atom["resseq"]),
                icode=str(atom["icode"]),
                x=float(atom["x"]),
                y=float(atom["y"]),
                z=float(atom["z"]),
                occupancy=1.0,
                bfactor=float(atom["bfactor"]),
                element=element,
            )
        )
    output_lines.append(f"TER   {len(selected_atoms) + 1:>5}      {selected_atoms[-1]['resname']:>3} {selected_atoms[-1]['chain_id'][:1]}{selected_atoms[-1]['resseq']:>4}")
    output_lines.append("END")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    return {
        "source_file": str(source),
        "destination_file": str(destination),
        "atom_count": len(selected_atoms),
        "removed_altloc_atoms": removed_altlocs,
    }


def validate_prepared_ligand_pdbqt(ligand_file: Path) -> List[LigandValidationIssue]:
    path = Path(ligand_file).expanduser()
    if not path.exists():
        return [LigandValidationIssue(ligand=path.name, prepared_file=str(path), reason="missing_prepared_ligand", details="Prepared ligand file does not exist.")]

    text = path.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()
    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    invalid_atom_types = sorted({_extract_atom_type(line) for line in atom_lines if _has_invalid_autodock_type(_extract_atom_type(line))})
    root_count = _count_line_token(lines, "ROOT")
    torsdof_count = sum(1 for line in lines if line.startswith("TORSDOF"))
    issues: List[LigandValidationIssue] = []

    if not atom_lines:
        issues.append(
            LigandValidationIssue(
                ligand=path.name,
                prepared_file=str(path),
                reason="missing_atom_records",
                details="No ATOM/HETATM records were found in the prepared ligand file.",
            )
        )
    if root_count != 1:
        issues.append(
            LigandValidationIssue(
                ligand=path.name,
                prepared_file=str(path),
                reason="invalid_root_count",
                details="Prepared ligand should contain exactly one ROOT block.",
                root_count=root_count,
                torsdof_count=torsdof_count,
            )
        )
    if torsdof_count != 1:
        issues.append(
            LigandValidationIssue(
                ligand=path.name,
                prepared_file=str(path),
                reason="invalid_torsdof_count",
                details="Prepared ligand should contain exactly one TORSDOF record.",
                root_count=root_count,
                torsdof_count=torsdof_count,
            )
        )
    if invalid_atom_types:
        issues.append(
            LigandValidationIssue(
                ligand=path.name,
                prepared_file=str(path),
                reason="invalid_atom_types",
                details=f"Prepared ligand contains invalid AutoDock atom types: {', '.join(invalid_atom_types)}",
                root_count=root_count,
                torsdof_count=torsdof_count,
            )
        )
    return issues


def audit_project_ligands(
    project_root: Path,
    rows: Iterable[PairlistRow],
    enable_admet_filters: bool = True,
    admet_thresholds: Optional[Dict[str, object]] = None,
    report_path: Optional[Path] = None,
) -> Dict[str, object]:
    root = Path(project_root).expanduser().resolve()
    manifest = load_manifest(root)
    ligands_dir = Path(manifest.get("ligands_dir") or "")
    raw_ligands_dir = Path(manifest.get("raw_ligands_dir") or "") if manifest.get("raw_ligands_dir") else None
    raw_ligands_sdf_dir = Path(manifest.get("raw_ligands_sdf_dir") or "") if manifest.get("raw_ligands_sdf_dir") else None

    pair_tags_by_ligand: Dict[str, List[str]] = defaultdict(list)
    for row in rows:
        pair_tags_by_ligand[row.ligand].append(row.tag)

    issues: List[LigandValidationIssue] = []
    warnings: List[str] = []
    admet_rows: List[Dict[str, object]] = []
    normalized_admet = _normalize_ligand_admet_thresholds(admet_thresholds)
    normalized_admet["enable_admet_filters"] = bool(
        enable_admet_filters and normalized_admet.get("enable_admet_filters", True)
    )

    for ligand_name in sorted(pair_tags_by_ligand):
        prepared_file = ligands_dir / ligand_name
        ligand_issues = validate_prepared_ligand_pdbqt(prepared_file)
        for issue in ligand_issues:
            issue.affected_tags = pair_tags_by_ligand[ligand_name]
            issue.affected_pairs = len(issue.affected_tags)
            if raw_ligands_dir:
                raw_pdb = raw_ligands_dir / f"{Path(ligand_name).stem}.pdb"
                if raw_pdb.exists():
                    issue.raw_pdb_file = str(raw_pdb)
            issues.append(issue)

        if normalized_admet["enable_admet_filters"]:
            raw_sdf = raw_ligands_sdf_dir / f"{Path(ligand_name).stem}.sdf" if raw_ligands_sdf_dir else None
            raw_pdb = raw_ligands_dir / f"{Path(ligand_name).stem}.pdb" if raw_ligands_dir else None
            mol, source, warning = _load_ligand_mol_from_sources(
                ligand_name=ligand_name,
                prepared_file=prepared_file,
                raw_sdf=raw_sdf,
                raw_pdb=raw_pdb,
            )
            if warning:
                warnings.append(warning)
                continue
            if mol is None:
                warnings.append(f"ADMET filtering skipped for {ligand_name}: no parseable chemistry source.")
                continue

            admet_issues, metrics = _assess_ligand_admet_issues(
                ligand_name=ligand_name,
                prepared_file=prepared_file,
                mol=mol,
                thresholds=normalized_admet,
            )
            for issue in admet_issues:
                issue.affected_tags = pair_tags_by_ligand[ligand_name]
                issue.affected_pairs = len(issue.affected_tags)
                if raw_ligands_dir:
                    raw_pdb_for_issue = raw_ligands_dir / f"{Path(ligand_name).stem}.pdb"
                    if raw_pdb_for_issue.exists():
                        issue.raw_pdb_file = str(raw_pdb_for_issue)
                issues.append(issue)
            admet_rows.append(
                {
                    "ligand": ligand_name,
                    "prepared_file": str(prepared_file),
                    "source": source,
                    **metrics,
                    "admet_issue_count": len(admet_issues),
                }
            )

    payload = {
        "project_root": str(root),
        "ligands_dir": str(ligands_dir),
        "issue_count": len(issues),
        "affected_ligands": sorted({issue.ligand for issue in issues}),
        "issues": [issue.to_dict() for issue in issues],
        "warning_count": len(warnings),
        "warnings": warnings,
        "admet_filters_enabled": bool(normalized_admet["enable_admet_filters"]),
        "admet_thresholds": normalized_admet,
        "admet_metrics": admet_rows,
    }
    if report_path:
        report = Path(report_path).expanduser()
        try:
            report.parent.mkdir(parents=True, exist_ok=True)
            report.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            payload["report_path"] = str(report)
        except OSError:
            payload["report_path"] = str(report)
            payload["report_write_error"] = True
    return payload


def _run_command(command: List[str]) -> None:
    subprocess.run(command, check=True, capture_output=True, text=True)


def _prepare_ligand_pdbqt(input_sdf: Path, output_pdbqt: Path) -> str:
    summary = prepare_ligand_for_vina_family(input_sdf, output_pdbqt)
    return str(summary.get("preparation_method", ""))


def _extract_smiles_from_prepared_ligand(prepared_file: Path) -> str:
    path = Path(prepared_file).expanduser()
    if not path.exists():
        return ""
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = re.match(r"^REMARK\s+SMILES\s+(.+?)\s*$", line.strip())
        if match:
            smiles = match.group(1).strip()
            if smiles and " IDX " not in smiles and " H PARENT " not in smiles:
                return smiles
    return ""


def _load_ligand_mol_from_sources(
    ligand_name: str,
    prepared_file: Path,
    raw_sdf: Optional[Path],
    raw_pdb: Optional[Path],
) -> tuple[object, str, str]:
    if not _RDKIT_MODULES:
        return None, "", "RDKit is not available for ligand ADMET filtering."

    Chem = _RDKIT_MODULES["Chem"]

    if raw_sdf and raw_sdf.exists():
        try:
            supplier = Chem.SDMolSupplier(str(raw_sdf), removeHs=False)
            if supplier and len(supplier) > 0:
                for mol in supplier:
                    if mol is not None:
                        return mol, "raw_sdf", ""
        except Exception:
            pass

    if raw_pdb and raw_pdb.exists():
        try:
            mol = Chem.MolFromPDBFile(str(raw_pdb), removeHs=False, sanitize=True)
            if mol is not None:
                return mol, "raw_pdb", ""
        except Exception:
            pass

    smiles = _extract_smiles_from_prepared_ligand(prepared_file)
    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                return Chem.AddHs(mol), "prepared_smiles", ""
        except Exception:
            pass

    return None, "", f"Could not resolve ADMET chemistry source for ligand {ligand_name}."


def _collect_reactive_alerts(mol) -> List[str]:
    if not _RDKIT_MODULES:
        return []
    Chem = _RDKIT_MODULES["Chem"]
    alerts: List[str] = []
    for label, pattern in REACTIVE_SMARTS_PATTERNS:
        try:
            smarts = Chem.MolFromSmarts(pattern)
            if smarts is not None and mol.HasSubstructMatch(smarts):
                alerts.append(label)
        except Exception:
            continue
    return sorted(set(alerts))


def _assess_ligand_admet_issues(
    ligand_name: str,
    prepared_file: Path,
    mol,
    thresholds: Dict[str, object],
) -> tuple[List[LigandValidationIssue], Dict[str, object]]:
    Crippen = _RDKIT_MODULES["Crippen"]
    Descriptors = _RDKIT_MODULES["Descriptors"]
    Lipinski = _RDKIT_MODULES["Lipinski"]
    rdMolDescriptors = _RDKIT_MODULES["rdMolDescriptors"]

    issues: List[LigandValidationIssue] = []

    molecular_weight = float(Descriptors.MolWt(mol))
    logp = float(Crippen.MolLogP(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    rotatable_bonds = int(Lipinski.NumRotatableBonds(mol))
    heavy_atoms = int(rdMolDescriptors.CalcNumHeavyAtoms(mol))
    h_donors = int(rdMolDescriptors.CalcNumHBD(mol))
    h_acceptors = int(rdMolDescriptors.CalcNumHBA(mol))
    formal_charge = int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))

    lipinski_violations = 0
    lipinski_violations += int(molecular_weight > 500.0)
    lipinski_violations += int(logp > 5.0)
    lipinski_violations += int(h_donors > 5)
    lipinski_violations += int(h_acceptors > 10)

    pains_alerts: List[str] = []
    brenk_alerts: List[str] = []
    if bool(thresholds.get("block_pains", True)):
        catalog = _get_filter_catalog("pains")
        if catalog is not None:
            try:
                pains_alerts = [str(match.GetDescription() or "PAINS") for match in catalog.GetMatches(mol)]
            except Exception:
                pains_alerts = []

    if bool(thresholds.get("block_brenk", True)):
        catalog = _get_filter_catalog("brenk")
        if catalog is not None:
            try:
                brenk_alerts = [str(match.GetDescription() or "BRENK") for match in catalog.GetMatches(mol)]
            except Exception:
                brenk_alerts = []

    reactive_alerts = _collect_reactive_alerts(mol) if bool(thresholds.get("block_reactive", True)) else []

    if lipinski_violations > int(thresholds["max_lipinski_violations"]):
        issues.append(
            LigandValidationIssue(
                ligand=ligand_name,
                prepared_file=str(prepared_file),
                reason="admet_lipinski_violations",
                details=(
                    f"Lipinski violations={lipinski_violations} exceed "
                    f"max_lipinski_violations={int(thresholds['max_lipinski_violations'])}."
                ),
            )
        )
    if molecular_weight > float(thresholds["max_molecular_weight"]):
        issues.append(
            LigandValidationIssue(
                ligand=ligand_name,
                prepared_file=str(prepared_file),
                reason="admet_excess_molecular_weight",
                details=(
                    f"molecular_weight={molecular_weight:.2f} exceeds "
                    f"max_molecular_weight={float(thresholds['max_molecular_weight']):.2f}."
                ),
            )
        )
    if logp > float(thresholds["max_logp"]):
        issues.append(
            LigandValidationIssue(
                ligand=ligand_name,
                prepared_file=str(prepared_file),
                reason="admet_excess_logp",
                details=f"logp={logp:.3f} exceeds max_logp={float(thresholds['max_logp']):.3f}.",
            )
        )
    if tpsa > float(thresholds["max_tpsa"]):
        issues.append(
            LigandValidationIssue(
                ligand=ligand_name,
                prepared_file=str(prepared_file),
                reason="admet_excess_tpsa",
                details=f"tpsa={tpsa:.3f} exceeds max_tpsa={float(thresholds['max_tpsa']):.3f}.",
            )
        )
    if rotatable_bonds > int(thresholds["max_rotatable_bonds"]):
        issues.append(
            LigandValidationIssue(
                ligand=ligand_name,
                prepared_file=str(prepared_file),
                reason="admet_excess_rotatable_bonds",
                details=(
                    f"rotatable_bonds={rotatable_bonds} exceeds "
                    f"max_rotatable_bonds={int(thresholds['max_rotatable_bonds'])}."
                ),
            )
        )
    if abs(formal_charge) > int(thresholds["max_formal_charge_abs"]):
        issues.append(
            LigandValidationIssue(
                ligand=ligand_name,
                prepared_file=str(prepared_file),
                reason="admet_excess_formal_charge",
                details=(
                    f"abs(formal_charge)={abs(formal_charge)} exceeds "
                    f"max_formal_charge_abs={int(thresholds['max_formal_charge_abs'])}."
                ),
            )
        )
    if heavy_atoms < int(thresholds["min_heavy_atom_count"]):
        issues.append(
            LigandValidationIssue(
                ligand=ligand_name,
                prepared_file=str(prepared_file),
                reason="admet_low_heavy_atom_count",
                details=(
                    f"heavy_atom_count={heavy_atoms} below "
                    f"min_heavy_atom_count={int(thresholds['min_heavy_atom_count'])}."
                ),
            )
        )
    if bool(thresholds.get("block_pains", True)) and pains_alerts:
        issues.append(
            LigandValidationIssue(
                ligand=ligand_name,
                prepared_file=str(prepared_file),
                reason="admet_pains_alert",
                details=f"PAINS alerts: {', '.join(sorted(set(pains_alerts))[:5])}",
            )
        )
    if bool(thresholds.get("block_brenk", True)) and brenk_alerts:
        issues.append(
            LigandValidationIssue(
                ligand=ligand_name,
                prepared_file=str(prepared_file),
                reason="admet_brenk_alert",
                details=f"Brenk alerts: {', '.join(sorted(set(brenk_alerts))[:5])}",
            )
        )
    if bool(thresholds.get("block_reactive", True)) and reactive_alerts:
        issues.append(
            LigandValidationIssue(
                ligand=ligand_name,
                prepared_file=str(prepared_file),
                reason="admet_reactive_alert",
                details=f"Reactive alerts: {', '.join(reactive_alerts)}",
            )
        )

    metrics = {
        "molecular_weight": molecular_weight,
        "logp": logp,
        "tpsa": tpsa,
        "rotatable_bonds": rotatable_bonds,
        "heavy_atom_count": heavy_atoms,
        "h_donors": h_donors,
        "h_acceptors": h_acceptors,
        "formal_charge": formal_charge,
        "lipinski_violations": lipinski_violations,
        "pains_alert_count": len(pains_alerts),
        "brenk_alert_count": len(brenk_alerts),
        "reactive_alert_count": len(reactive_alerts),
    }
    return issues, metrics


def _write_sdf_from_smiles(smiles: str, output_sdf: Path) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse ligand SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    embed_status = AllChem.EmbedMolecule(mol, params)
    if embed_status != 0:
        raise ValueError("RDKit could not embed a 3D conformer from ligand SMILES")
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception:
            pass
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_sdf))
    try:
        writer.write(mol)
    finally:
        writer.close()


def repair_project_ligands(
    project_root: Path,
    ligand_names: Sequence[str],
    preferred_altloc: str = "A",
    backup_dir: Optional[Path] = None,
) -> Dict[str, object]:
    root = Path(project_root).expanduser().resolve()
    manifest = load_manifest(root)
    layout_profile = detect_layout_profile(root)
    layout = ensure_project_layout(root, layout_profile)

    raw_ligands_dir = Path(manifest.get("raw_ligands_dir") or layout["raw_ligands"]).expanduser().resolve()
    raw_ligands_sdf_dir = Path(manifest.get("raw_ligands_sdf_dir") or layout["raw_ligands_sdf"]).expanduser().resolve()
    prepared_ligands_dir = Path(manifest.get("prepared_ligands_dir") or layout["prepared_ligands"]).expanduser().resolve()
    docking_ligands_dir = Path(manifest.get("ligands_dir") or layout["ligands"]).expanduser().resolve()
    metadata_dir = Path(manifest.get("metadata_dir") or layout["metadata"]).expanduser().resolve()

    backup_root = Path(backup_dir).expanduser().resolve() if backup_dir else metadata_dir / "ligand_repair_backups"
    repaired: List[Dict[str, object]] = []

    for ligand_name in ligand_names:
        stem = Path(ligand_name).stem
        raw_pdb = raw_ligands_dir / f"{stem}.pdb"
        raw_sdf = raw_ligands_sdf_dir / f"{stem}.sdf"
        prepared_pdbqt = prepared_ligands_dir / f"{stem}.pdbqt"
        docking_pdbqt = docking_ligands_dir / f"{stem}.pdbqt"

        stamp_dir = backup_root / stem
        stamp_dir.mkdir(parents=True, exist_ok=True)
        for original in (raw_pdb, raw_sdf, prepared_pdbqt, docking_pdbqt):
            if original.exists():
                backup_path = stamp_dir / original.name
                if not backup_path.exists():
                    shutil.copy2(original, backup_path)

        with tempfile.TemporaryDirectory(prefix=f"ligand_repair_{stem}_") as tmp_dir:
            temp_dir = Path(tmp_dir)
            sanitized_pdb = temp_dir / f"{stem}.pdb"
            sanitized_sdf = temp_dir / f"{stem}.sdf"
            sanitized_pdbqt = temp_dir / f"{stem}.pdbqt"
            source_mode = "raw_pdb"
            if raw_pdb.exists():
                sanitize_summary = sanitize_ligand_pdb_file(raw_pdb, sanitized_pdb, preferred_altloc=preferred_altloc)
                _run_command(["obabel", str(sanitized_pdb), "-O", str(sanitized_sdf)])
            elif raw_sdf.exists():
                source_mode = "raw_sdf"
                sanitize_summary = {
                    "source_file": str(raw_sdf),
                    "destination_file": str(sanitized_sdf),
                    "atom_count": 0,
                    "removed_altloc_atoms": 0,
                }
                shutil.copy2(raw_sdf, sanitized_sdf)
                _run_command(["obabel", str(sanitized_sdf), "-O", str(sanitized_pdb)])
            else:
                smiles = _extract_smiles_from_prepared_ligand(prepared_pdbqt if prepared_pdbqt.exists() else docking_pdbqt)
                if not smiles:
                    raise FileNotFoundError(
                        f"Could not repair {ligand_name}: no raw PDB, no raw SDF, and no usable REMARK SMILES were found"
                    )
                source_mode = "prepared_smiles"
                _write_sdf_from_smiles(smiles, sanitized_sdf)
                _run_command(["obabel", str(sanitized_sdf), "-O", str(sanitized_pdb)])
                sanitize_summary = {
                    "source_file": str(prepared_pdbqt if prepared_pdbqt.exists() else docking_pdbqt),
                    "destination_file": str(sanitized_sdf),
                    "atom_count": 0,
                    "removed_altloc_atoms": 0,
                    "smiles": smiles,
                }
            prep_method = _prepare_ligand_pdbqt(sanitized_sdf, sanitized_pdbqt)
            atom_type_repair = sanitize_prepared_ligand_pdbqt(sanitized_pdbqt, sanitized_pdbqt)
            issues = validate_prepared_ligand_pdbqt(sanitized_pdbqt)
            if issues:
                raise ValueError(f"Repaired ligand {ligand_name} is still invalid: {[issue.reason for issue in issues]}")

            shutil.copy2(sanitized_pdb, raw_pdb)
            shutil.copy2(sanitized_sdf, raw_sdf)
            shutil.copy2(sanitized_pdbqt, prepared_pdbqt)
            if docking_pdbqt.exists() or docking_pdbqt.is_symlink():
                docking_pdbqt.unlink()
            shutil.copy2(sanitized_pdbqt, docking_pdbqt)
            repaired.append(
                {
                    "ligand": ligand_name,
                    "source_mode": source_mode,
                    "raw_pdb": str(raw_pdb),
                    "raw_sdf": str(raw_sdf),
                    "prepared_pdbqt": str(prepared_pdbqt),
                    "docking_pdbqt": str(docking_pdbqt),
                    "prep_method": prep_method,
                    "sanitize_summary": sanitize_summary,
                    "atom_type_repair": atom_type_repair,
                }
            )

    report_path = metadata_dir / "ligand_repair_report.json"
    payload = {"project_root": str(root), "repaired_ligands": repaired, "backup_dir": str(backup_root)}
    report_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    payload["report_path"] = str(report_path)
    return payload
