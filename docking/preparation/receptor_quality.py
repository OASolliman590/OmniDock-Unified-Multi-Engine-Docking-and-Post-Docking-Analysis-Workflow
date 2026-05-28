from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ..models import PairlistRow
from ..project_layout import load_manifest


@dataclass
class ReceptorValidationIssue:
    receptor: str
    receptor_file: str
    reason: str
    details: str = ""
    atom_count: int = 0
    heavy_atom_count: int = 0
    chain_count: int = 0
    max_coordinate_span: float = 0.0
    affected_pairs: int = 0
    affected_tags: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


def _safe_float(token: str) -> Optional[float]:
    try:
        return float(str(token).strip())
    except Exception:
        return None


def _parse_xyz(line: str) -> Optional[Tuple[float, float, float]]:
    if len(line) >= 54:
        x = _safe_float(line[30:38])
        y = _safe_float(line[38:46])
        z = _safe_float(line[46:54])
        if x is not None and y is not None and z is not None:
            return x, y, z
    tokens = line.split()
    float_values: List[float] = []
    for token in tokens:
        value = _safe_float(token)
        if value is not None:
            float_values.append(value)
        if len(float_values) >= 3:
            return float_values[0], float_values[1], float_values[2]
    return None


def _infer_element(line: str) -> str:
    element_field = str(line[76:78]).strip() if len(line) >= 78 else ""
    if element_field:
        return "".join(char for char in element_field if char.isalpha()).upper()
    atom_name = str(line[12:16]).strip() if len(line) >= 16 else ""
    atom_letters = "".join(char for char in atom_name if char.isalpha()).upper()
    if not atom_letters:
        return ""
    if len(atom_letters) >= 2 and atom_letters[:2] in {"CL", "BR", "ZN", "FE", "MG", "MN", "CU", "NA", "CA", "PD"}:
        return atom_letters[:2]
    return atom_letters[0]


def _is_heavy_atom(line: str) -> bool:
    element = _infer_element(line)
    return bool(element and element != "H")


def validate_receptor_file(
    receptor_file: Path,
    *,
    min_atom_count: int = 100,
    min_heavy_atom_count: int = 60,
    min_chain_count: int = 1,
    max_coordinate_span: float = 500.0,
) -> List[ReceptorValidationIssue]:
    path = Path(receptor_file).expanduser()
    receptor_name = path.name
    if not path.exists():
        return [
            ReceptorValidationIssue(
                receptor=receptor_name,
                receptor_file=str(path),
                reason="missing_receptor_file",
                details="Prepared receptor file does not exist.",
            )
        ]

    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    atom_count = len(atom_lines)
    heavy_atom_count = sum(1 for line in atom_lines if _is_heavy_atom(line))

    chains = set()
    for line in atom_lines:
        chain = str(line[21]).strip() if len(line) >= 22 else ""
        chains.add(chain or "A")
    chain_count = len(chains)

    coords = [_parse_xyz(line) for line in atom_lines]
    valid_coords = [coord for coord in coords if coord is not None]
    if valid_coords:
        xs = [coord[0] for coord in valid_coords]
        ys = [coord[1] for coord in valid_coords]
        zs = [coord[2] for coord in valid_coords]
        span = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
    else:
        span = 0.0

    issues: List[ReceptorValidationIssue] = []
    if atom_count == 0:
        issues.append(
            ReceptorValidationIssue(
                receptor=receptor_name,
                receptor_file=str(path),
                reason="missing_atom_records",
                details="No ATOM/HETATM records were found in receptor file.",
                atom_count=atom_count,
                heavy_atom_count=heavy_atom_count,
                chain_count=chain_count,
                max_coordinate_span=float(span),
            )
        )
        return issues

    if atom_count < int(min_atom_count):
        issues.append(
            ReceptorValidationIssue(
                receptor=receptor_name,
                receptor_file=str(path),
                reason="low_atom_count",
                details=f"Receptor atom_count={atom_count} is below min_atom_count={int(min_atom_count)}.",
                atom_count=atom_count,
                heavy_atom_count=heavy_atom_count,
                chain_count=chain_count,
                max_coordinate_span=float(span),
            )
        )
    if heavy_atom_count < int(min_heavy_atom_count):
        issues.append(
            ReceptorValidationIssue(
                receptor=receptor_name,
                receptor_file=str(path),
                reason="low_heavy_atom_count",
                details=f"Receptor heavy_atom_count={heavy_atom_count} is below min_heavy_atom_count={int(min_heavy_atom_count)}.",
                atom_count=atom_count,
                heavy_atom_count=heavy_atom_count,
                chain_count=chain_count,
                max_coordinate_span=float(span),
            )
        )
    if chain_count < int(min_chain_count):
        issues.append(
            ReceptorValidationIssue(
                receptor=receptor_name,
                receptor_file=str(path),
                reason="low_chain_count",
                details=f"Receptor chain_count={chain_count} is below min_chain_count={int(min_chain_count)}.",
                atom_count=atom_count,
                heavy_atom_count=heavy_atom_count,
                chain_count=chain_count,
                max_coordinate_span=float(span),
            )
        )
    if float(span) > float(max_coordinate_span):
        issues.append(
            ReceptorValidationIssue(
                receptor=receptor_name,
                receptor_file=str(path),
                reason="high_coordinate_span",
                details=(
                    f"Receptor max_coordinate_span={float(span):.3f} exceeds "
                    f"max_coordinate_span={float(max_coordinate_span):.3f}."
                ),
                atom_count=atom_count,
                heavy_atom_count=heavy_atom_count,
                chain_count=chain_count,
                max_coordinate_span=float(span),
            )
        )
    return issues


def audit_project_receptors(
    project_root: Path,
    rows: Iterable[PairlistRow],
    *,
    min_atom_count: int = 100,
    min_heavy_atom_count: int = 60,
    min_chain_count: int = 1,
    max_coordinate_span: float = 500.0,
    report_path: Optional[Path] = None,
) -> Dict[str, object]:
    root = Path(project_root).expanduser().resolve()
    manifest = load_manifest(root)
    receptors_dir = Path(manifest.get("receptors_dir") or "")

    pair_tags_by_receptor: Dict[str, List[str]] = {}
    for row in rows:
        pair_tags_by_receptor.setdefault(row.receptor, []).append(row.tag)

    issues: List[ReceptorValidationIssue] = []
    for receptor_name in sorted(pair_tags_by_receptor):
        receptor_file = receptors_dir / receptor_name
        receptor_issues = validate_receptor_file(
            receptor_file,
            min_atom_count=min_atom_count,
            min_heavy_atom_count=min_heavy_atom_count,
            min_chain_count=min_chain_count,
            max_coordinate_span=max_coordinate_span,
        )
        for issue in receptor_issues:
            issue.affected_tags = pair_tags_by_receptor[receptor_name]
            issue.affected_pairs = len(issue.affected_tags)
            issues.append(issue)

    payload = {
        "project_root": str(root),
        "receptors_dir": str(receptors_dir),
        "issue_count": len(issues),
        "affected_receptors": sorted({issue.receptor for issue in issues}),
        "issues": [issue.to_dict() for issue in issues],
        "thresholds": {
            "min_atom_count": int(min_atom_count),
            "min_heavy_atom_count": int(min_heavy_atom_count),
            "min_chain_count": int(min_chain_count),
            "max_coordinate_span": float(max_coordinate_span),
        },
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
