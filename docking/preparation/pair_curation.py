from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

from ..project_layout import (
    LAYOUT_DOCKING_LEGACY,
    detect_layout_profile,
    ensure_project_layout,
    pair_curation_state_path,
    pair_intent_path,
    pairlist_path,
)
from .pairlist_builder import (
    PAIRLIST_COLUMNS,
    PAIR_INTENT_COLUMNS,
    _write_csv,
    generate_pairlists,
)


DEFAULT_ROUND_ID = "round_001"


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _default_state(project_root: Path, layout_profile: str) -> Dict[str, object]:
    timestamp = _now_iso()
    return {
        "version": 1,
        "project_root": str(project_root),
        "layout_profile": layout_profile,
        "created_at": timestamp,
        "updated_at": timestamp,
        "current_round": "",
        "rounds": {},
    }


def load_pair_curation_state(project_root: Path, layout_profile: Optional[str] = None) -> Dict[str, object]:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    path = pair_curation_state_path(root, profile)
    if not path.exists():
        return _default_state(root, profile)
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def save_pair_curation_state(
    project_root: Path,
    payload: Dict[str, object],
    layout_profile: Optional[str] = None,
) -> Path:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    ensure_project_layout(root, profile)
    payload["updated_at"] = _now_iso()
    payload.setdefault("project_root", str(root))
    payload.setdefault("layout_profile", profile)
    path = pair_curation_state_path(root, profile)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    return path


def next_round_id(project_root: Path, layout_profile: Optional[str] = None) -> str:
    state = load_pair_curation_state(project_root, layout_profile)
    existing = []
    for round_id in state.get("rounds", {}):
        if str(round_id).startswith("round_"):
            try:
                existing.append(int(str(round_id).split("_", 1)[1]))
            except ValueError:
                continue
    next_index = max(existing, default=0) + 1
    return f"round_{next_index:03d}"


def upsert_pair_curation_round(
    project_root: Path,
    prepared_proteins: Path,
    prepared_ligands: Path,
    excel_path: Path,
    mode: str,
    default_site_id: str = "site_1",
    default_box_size: float = 20.0,
    curated_receptors: Optional[List[str]] = None,
    curated_ligands: Optional[List[str]] = None,
    curated_mapping: Optional[Dict[str, List[str]]] = None,
    layout_profile: str = LAYOUT_DOCKING_LEGACY,
    round_id: Optional[str] = None,
    freeze: bool = False,
) -> Dict[str, object]:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    summary = generate_pairlists(
        project_root=root,
        prepared_proteins=prepared_proteins,
        prepared_ligands=prepared_ligands,
        excel_path=excel_path,
        mode=mode,
        default_site_id=default_site_id,
        default_box_size=default_box_size,
        curated_receptors=curated_receptors,
        curated_ligands=curated_ligands,
        curated_mapping=curated_mapping,
        layout_profile=profile,
    )
    state = load_pair_curation_state(root, profile)
    resolved_round = round_id or str(state.get("current_round") or next_round_id(root, profile))
    rounds = dict(state.get("rounds", {}))
    existing = dict(rounds.get(resolved_round, {}))
    rounds[resolved_round] = {
        **existing,
        "round_id": resolved_round,
        "mode": mode,
        "default_site_id": default_site_id,
        "default_box_size": default_box_size,
        "prepared_proteins": str(Path(prepared_proteins).expanduser().resolve()),
        "prepared_ligands": str(Path(prepared_ligands).expanduser().resolve()),
        "excel_path": str(Path(excel_path).expanduser().resolve()),
        "curated_receptors": curated_receptors or [],
        "curated_ligands": curated_ligands or [],
        "curated_mapping": curated_mapping or {},
        "pair_count": summary["pair_count"],
        "receptor_count": summary["receptor_count"],
        "ligand_count": summary["ligand_count"],
        "warning_count": summary["warning_count"],
        "warnings": list(summary.get("warnings", [])),
        "has_cocrystal_benchmark_rows": bool(summary.get("has_cocrystal_benchmark_rows", False)),
        "pairlist_rows": summary["pairlist_rows"],
        "pair_intent_rows": summary["pair_intent_rows"],
        "updated_at": _now_iso(),
    }
    state["rounds"] = rounds
    state["current_round"] = resolved_round
    state_path = save_pair_curation_state(root, state, profile)

    result = {
        **summary,
        "pair_curation_state_file": str(state_path),
        "round_id": resolved_round,
        "current_round": resolved_round,
    }
    if freeze:
        freeze_summary = materialize_pair_curation_round(root, resolved_round, layout_profile=profile)
        result.update(freeze_summary)
    return result


def materialize_pair_curation_round(
    project_root: Path,
    round_id: Optional[str] = None,
    layout_profile: Optional[str] = None,
) -> Dict[str, object]:
    root = Path(project_root).expanduser().resolve()
    profile = detect_layout_profile(root, layout_profile)
    state = load_pair_curation_state(root, profile)
    resolved_round = round_id or str(state.get("current_round") or "")
    if not resolved_round:
        raise ValueError("No pair-curation round is available to freeze")

    rounds = state.get("rounds", {})
    if resolved_round not in rounds:
        raise ValueError(f"Pair-curation round '{resolved_round}' was not found")

    round_payload = dict(rounds[resolved_round])
    pairlist_rows = list(round_payload.get("pairlist_rows", []))
    pair_intent_rows = list(round_payload.get("pair_intent_rows", []))
    if not pairlist_rows:
        raise ValueError(f"Pair-curation round '{resolved_round}' has no pairlist rows")

    pairlist_file = pairlist_path(root, profile)
    intent_file = pair_intent_path(root, profile)
    _write_csv(intent_file, pair_intent_rows, PAIR_INTENT_COLUMNS)
    _write_csv(pairlist_file, pairlist_rows, PAIRLIST_COLUMNS)

    round_payload["frozen_at"] = _now_iso()
    round_payload["pairlist_file"] = str(pairlist_file)
    round_payload["pair_intent_file"] = str(intent_file)
    rounds[resolved_round] = round_payload
    state["rounds"] = rounds
    state["current_round"] = resolved_round
    state_path = save_pair_curation_state(root, state, profile)

    return {
        "project_root": str(root),
        "layout_profile": profile,
        "round_id": resolved_round,
        "pairlist_file": str(pairlist_file),
        "pair_intent_file": str(intent_file),
        "pair_curation_state_file": str(state_path),
        "pair_count": len(pairlist_rows),
        "receptor_count": len({row["receptor"] for row in pairlist_rows}),
        "ligand_count": len({row["ligand"] for row in pairlist_rows}),
        "warning_count": int(round_payload.get("warning_count", 0)),
        "warnings": list(round_payload.get("warnings", [])),
        "pair_mode": str(round_payload.get("mode", "")),
        "has_cocrystal_benchmark_rows": bool(round_payload.get("has_cocrystal_benchmark_rows", False)),
    }
