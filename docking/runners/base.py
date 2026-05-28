from __future__ import annotations

import json
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from ..models import EnsembleReceptorConfig, EngineJobResult, PairlistRow, normalize_ensemble_aggregation_strategy
from ..project_layout import ensure_engine_layout, pairlist_path, shared_ligands_dir, shared_receptors_dir


class DockingEngineRunner(ABC):
    name: str = ""
    binary_name: str = ""
    pose_extension: str = ""

    def __init__(self, project_root: Path, runtime: Optional[Dict[str, object]] = None):
        self.project_root = Path(project_root)
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

        for job in jobs:
            if job.status != "planned":
                continue
            runner_log = self.layout["logs"] / f"{job.tag}.runner.log"
            pair_log = Path(job.log_file)
            cwd_value = Path(str(self.runtime.get("execution_workdir") or self.project_root))
            completed = subprocess.run(
                job.command,
                cwd=cwd_value,
                capture_output=True,
                text=True,
            )
            if not pair_log.exists() or pair_log.stat().st_size == 0:
                pair_log.write_text(
                    f"COMMAND: {' '.join(job.command)}\n\nSTDOUT:\n{completed.stdout}\n\nSTDERR:\n{completed.stderr}",
                    encoding="utf-8",
                )
            runner_log.write_text(
                f"COMMAND: {' '.join(job.command)}\n\nSTDOUT:\n{completed.stdout}\n\nSTDERR:\n{completed.stderr}",
                encoding="utf-8",
            )
            job.status = "completed" if completed.returncode == 0 else "failed"
            job.returncode = completed.returncode
            job.error = completed.stderr.strip()

        payload = {
            "engine": self.name,
            "runtime": self.runtime,
            "jobs": [job.to_dict() for job in jobs],
        }
        manifest_path = self.layout["root"] / "run_manifest.json"
        with open(manifest_path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)

        if not dry_run:
            normalized = self.collect_normalized_scores(pairlist_rows)
            if normalized is not None and not normalized.empty:
                normalized.to_csv(self.layout["scores"] / "normalized_scores.csv", index=False)

        return payload

    def plan_jobs(
        self,
        pairlist_rows: List[PairlistRow],
        skip_completed: bool = False,
    ) -> List[EngineJobResult]:
        jobs: List[EngineJobResult] = []
        for row in pairlist_rows:
            pose_file = self.layout["poses"] / f"{row.tag}{self.pose_extension}"
            log_file = self.layout["logs"] / f"{row.tag}.log"
            command = self.build_command(row, pose_file=pose_file, log_file=log_file)
            status = "planned"
            if skip_completed and pose_file.exists():
                status = "skipped"
            jobs.append(
                EngineJobResult(
                    tag=row.tag,
                    command=command,
                    pose_file=str(pose_file),
                    log_file=str(log_file),
                    status=status,
                )
            )
        return jobs

    @abstractmethod
    def build_command(self, row: PairlistRow, pose_file: Path, log_file: Path) -> List[str]:
        raise NotImplementedError

    @abstractmethod
    def collect_normalized_scores(self, pairlist_rows: List[PairlistRow]) -> pd.DataFrame:
        raise NotImplementedError
