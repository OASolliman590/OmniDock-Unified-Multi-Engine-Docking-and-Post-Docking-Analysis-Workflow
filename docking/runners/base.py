from __future__ import annotations

import json
import hashlib
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from ..models import EnsembleReceptorConfig, EngineJobResult, PairlistRow, normalize_ensemble_aggregation_strategy
from ..project_layout import ensure_engine_layout, pairlist_path, shared_ligands_dir, shared_receptors_dir
from ..parameter_schema import validate_pairlist_geometry
from .job_contract import completion_valid, execute_job, file_hash, validate_pose_file, write_json


class DockingEngineRunner(ABC):
    name: str = ""
    binary_name: str = ""
    pose_extension: str = ""

    def __init__(self, project_root: Path, runtime: Optional[Dict[str, object]] = None):
        self.project_root = Path(project_root).expanduser().resolve()
        self.runtime = runtime or {}
        self.layout = ensure_engine_layout(self.project_root, self.name)
        self.shared_receptors = shared_receptors_dir(self.project_root)
        self.shared_ligands = shared_ligands_dir(self.project_root)
        self.pairlist_file = pairlist_path(self.project_root)
        self.ensemble_config = self._resolve_ensemble_config(self.runtime)

    @staticmethod
    def _resolve_ensemble_config(runtime: Optional[Dict[str, object]]) -> EnsembleReceptorConfig:
        payload = dict((runtime or {}).get("ensemble_receptors") or {})
        return EnsembleReceptorConfig(
            enabled=bool(payload.get("enabled", False)),
            strategy=normalize_ensemble_aggregation_strategy(payload.get("strategy")),
            conformer_paths=[str(path) for path in payload.get("conformer_paths", []) if str(path).strip()],
            max_conformers=max(int(payload.get("max_conformers", 0) or 0), 0),
            selection_note=str(payload.get("selection_note", "") or ""),
        )

    def supports_ensemble_receptors(self) -> bool:
        """
        Extension-point hook for future ensemble receptor execution.
        """
        return False

    def build_ensemble_plan(self, pairlist_rows: List[PairlistRow]) -> Dict[str, object]:
        """
        Return planning metadata for future multi-conformation runs.
        """
        if not self.ensemble_config.enabled:
            return {"enabled": False, "supported": self.supports_ensemble_receptors(), "reason": "disabled"}
        if not self.supports_ensemble_receptors():
            return {
                "enabled": True,
                "supported": False,
                "reason": f"{self.name} runner does not yet implement ensemble receptor execution",
                "config": self.ensemble_config.to_dict(),
                "pair_count": int(len(pairlist_rows)),
            }
        return {
            "enabled": True,
            "supported": True,
            "reason": "runner supports ensemble receptor planning",
            "config": self.ensemble_config.to_dict(),
            "pair_count": int(len(pairlist_rows)),
        }

    def run(
        self,
        pairlist_rows: List[PairlistRow],
        dry_run: bool = False,
        skip_completed: bool = False,
    ) -> Dict[str, object]:
        jobs = self.plan_jobs(pairlist_rows, skip_completed=skip_completed)
        if dry_run:
            for job in jobs:
                if job.status == "planned":
                    job.status = "dry_run"
            payload = {
                "engine": self.name,
                "runtime": self.runtime,
                "jobs": [job.to_dict() for job in jobs],
            }
            manifest_path = self.layout["root"] / "run_manifest.json"
            with open(manifest_path, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2)
            return payload

        manifest_path = self.layout["root"] / "run_manifest.json"
        # An interrupted invocation must not leave yesterday's accepted table
        # looking like today's results. Valid skipped jobs are restored below.
        pd.DataFrame(columns=["engine", "tag", "protein", "ligand", "site_id", "pose", "affinity_kcal_mol"]).to_csv(self.layout["scores"] / "normalized_scores.csv", index=False)
        if self.name == "gnina":
            pd.DataFrame(columns=["tag", "mode", "vina_affinity", "cnn_affinity", "cnn_score"]).to_csv(self.layout["scores"] / "all_scores.csv", index=False)
        for job in jobs:
            if job.status != "planned":
                continue
            cwd_value = Path(str(self.runtime.get("execution_workdir") or self.project_root))
            result = execute_job(job.to_dict(), workdir=cwd_value, skip_completed=skip_completed)
            job.status, job.returncode, job.error = result["status"], result["returncode"], result["error"]
            write_json(manifest_path, {"engine": self.name, "runtime": self.runtime, "jobs": [item.to_dict() for item in jobs]})

        payload = {
            "engine": self.name,
            "runtime": self.runtime,
            "jobs": [job.to_dict() for job in jobs],
        }
        manifest_path = self.layout["root"] / "run_manifest.json"
        with open(manifest_path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)

        accepted = {job.tag for job in jobs if job.status in {"completed", "skipped"}}
        normalized = self.collect_normalized_scores([row for row in pairlist_rows if row.tag in accepted])
        if normalized is None or normalized.empty:
            normalized = pd.DataFrame(columns=["engine", "tag", "protein", "ligand", "site_id", "pose", "affinity_kcal_mol", "score_name_primary", "score_primary", "score_name_secondary", "score_secondary", "rmsd_lb", "rmsd_ub", "pose_file", "log_file"])
        normalized.to_csv(self.layout["scores"] / "normalized_scores.csv", index=False)

        return payload

    def plan_jobs(
        self,
        pairlist_rows: List[PairlistRow],
        skip_completed: bool = False,
    ) -> List[EngineJobResult]:
        errors = validate_pairlist_geometry(pairlist_rows)
        if errors:
            raise ValueError("\n".join(errors))
        jobs: List[EngineJobResult] = []
        for row in pairlist_rows:
            pose_file = self.layout["poses"] / f"{row.tag}{self.pose_extension}"
            log_file = self.layout["logs"] / f"{row.tag}.log"
            command = self.build_command(row, pose_file=pose_file, log_file=log_file)
            inputs = self.input_files(row)
            # Ignore transport/round labels; scientific and executable settings
            # remain part of the protocol. Bytes, not file names, identify inputs.
            protocol = {key: value for key, value in self.runtime.items() if key not in {"execution_workdir", "deployment_preamble", "round_id", "rerun_manifest_file", "job_python"}}
            fingerprint = hashlib.sha256(json.dumps({
                "contract_version": 1, "engine": self.name, "pair": row.to_dict(), "runtime": protocol,
                "command": [token.replace(str(self.project_root), "<project>") for token in command],
                "inputs": [{"name": Path(item["path"]).name, "sha256": item["sha256"]} for item in inputs],
            }, sort_keys=True, allow_nan=False).encode()).hexdigest()
            job = EngineJobResult(
                    tag=row.tag,
                    command=command,
                    pose_file=str(pose_file),
                    log_file=str(log_file),
                    status="planned",
                    engine=self.name,
                    fingerprint=fingerprint,
                    completion_file=str(pose_file) + ".completion.json",
                    input_files=inputs,
                    executables=self.executables(),
                    skip_completed=skip_completed,
                )
            if skip_completed and completion_valid(job.to_dict()):
                job.status = "skipped"
            jobs.append(job)
        return jobs

    def input_files(self, row: PairlistRow) -> List[Dict[str, str]]:
        receptor = self._resolve_receptor_path(row.receptor) if hasattr(self, "_resolve_receptor_path") else self.shared_receptors / row.receptor
        ligand = self._resolve_ligand_path(row.ligand) if hasattr(self, "_resolve_ligand_path") else self.shared_ligands / row.ligand
        paths = [receptor, ligand]
        for key in ("image", "parameter_file_path", "parameter_file", "autodocktools_prepare_gpf4", "autodocktools_prepare_dpf4"):
            if self.runtime.get(key):
                path = Path(str(self.runtime[key])).expanduser()
                if not path.is_absolute():
                    path = self.project_root / path
                if path.is_file():
                    paths.append(path)
        return [{"path": str(path), "sha256": file_hash(path)} for path in paths]

    def executables(self) -> List[str]:
        binary = str(self.runtime.get("binary") or (self.runtime.get("autodock_binary") if self.name == "autodock4" else None) or self.binary_name)
        # A bare engine name inside `conda run` resolves in that environment,
        # not on the caller's PATH. Without an absolute engine path, provenance
        # is unknown and completion_valid deliberately declines reuse.
        if self.runtime.get("conda_env") and not Path(binary).is_absolute():
            return []
        binaries = [binary]
        if self.name == "autodock4":
            binaries.append(str(self.runtime.get("autogrid_binary") or "autogrid4"))
        if self.runtime.get("image"):
            image_path = Path(str(self.runtime["image"])).expanduser()
            if not image_path.is_absolute():
                image_path = self.project_root / image_path
            if not image_path.is_file():
                return []
            binaries = ["apptainer"]
        return binaries

    def accepted_pose_files(self, pairlist_rows: List[PairlistRow]):
        for row in pairlist_rows:
            path = self.layout["poses"] / f"{row.tag}{self.pose_extension}"
            try:
                validate_pose_file(path, self.name)
            except (OSError, ValueError, KeyError, IndexError, StopIteration):
                continue
            yield path

    @abstractmethod
    def build_command(self, row: PairlistRow, pose_file: Path, log_file: Path) -> List[str]:
        raise NotImplementedError

    @abstractmethod
    def collect_normalized_scores(self, pairlist_rows: List[PairlistRow]) -> pd.DataFrame:
        raise NotImplementedError
