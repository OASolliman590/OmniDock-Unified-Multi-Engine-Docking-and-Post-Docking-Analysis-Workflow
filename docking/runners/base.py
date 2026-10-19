from __future__ import annotations

import hashlib
import json
import time
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from ..models import EnsembleReceptorConfig, EngineJobResult, PairlistRow, normalize_ensemble_aggregation_strategy
from ..project_layout import ensure_engine_layout, pairlist_path, shared_ligands_dir, shared_receptors_dir
from .job_contract import completion_valid, execute_job, file_hash, validate_pose_file, write_json

# Spec 036 R5b naming helpers live in post_docking_analysis.replicate_names (shared with post-docking analysis).
from post_docking_analysis.replicate_names import (  # noqa: F401  (re-exported for runner modules and tests)
    REPLICATE_STEM_PATTERN,
    REPLICATE_TOKEN,
    replicate_job_tag,
    split_replicate_stem,
)

# Engines whose runner takes a seed and runs independent replicates. AutoDock4 is not in
# this set: its replicate execution is outside Spec 036 R5b and is recorded as replicates: 1.
REPLICATE_ENGINES = {"vina", "smina", "gnina"}
VINA_FAMILY_ENGINES = {"vina", "smina"}
# Execution-only keys: they change where or how a job runs, not the scientific protocol.
TRANSPORT_RUNTIME_KEYS = {"execution_workdir", "deployment_preamble", "round_id", "rerun_manifest_file", "job_python"}
NORMALIZED_SCORE_COLUMNS = [
    "engine", "tag", "protein", "ligand", "site_id", "pose", "affinity_kcal_mol", "score_name_primary",
    "score_primary", "score_name_secondary", "score_secondary", "rmsd_lb", "rmsd_ub", "pose_file", "log_file",
]


def build_pair_index(pairlist_rows: List[PairlistRow]) -> Dict[str, PairlistRow]:
    """Index pairs by current tag and by the pre-Spec 036 legacy tag, so old pose files still resolve."""
    index: Dict[str, PairlistRow] = {}
    for row in pairlist_rows:
        index.setdefault(row.tag, row)
        index.setdefault(row.legacy_tag, row)
    return index


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

    # ----- Spec 036 R5b: replicates and seeds -------------------------------------------
    @property
    def effective_replicates(self) -> int:
        if self.name not in REPLICATE_ENGINES:
            return 1
        requested = self.runtime.get("replicates")
        if requested in (None, ""):
            return 1
        count = int(requested)
        if count < 1:
            raise ValueError("`replicates` must be >= 1")
        return count

    def replicate_seeds(self) -> List[Optional[int]]:
        base_seed = self.runtime.get("seed")
        count = self.effective_replicates
        if base_seed in (None, ""):
            if count > 1:
                raise ValueError("Replicates require a base seed (`--seed`); none was given.")
            return [None] * count
        base = int(base_seed)
        return [base + offset for offset in range(count)]

    def seed_for_replicate(self, replicate_id: Optional[int]) -> Optional[int]:
        seeds = self.replicate_seeds()
        if replicate_id is None or not seeds or seeds[0] is None:
            return None
        index = int(replicate_id) - 1
        if index < 0 or index >= len(seeds):
            return None
        return seeds[index]

    # ----- Spec 036 R5c/R5d: effective engine parameters --------------------------------
    @property
    def effective_exhaustiveness(self) -> Optional[int]:
        value = self.runtime.get("exhaustiveness")
        return None if value in (None, "") else int(value)

    @property
    def effective_num_modes(self) -> Optional[int]:
        value = self.runtime.get("num_modes")
        return None if value in (None, "") else int(value)

    @property
    def effective_energy_range(self) -> Optional[float]:
        if self.name not in VINA_FAMILY_ENGINES:
            return None
        value = self.runtime.get("energy_range", 3.0)
        return None if value in (None, "") else float(value)

    def effective_settings(self) -> Dict[str, object]:
        seeds = self.replicate_seeds()
        return {
            "engine": self.name,
            "exhaustiveness": self.effective_exhaustiveness,
            "num_modes": self.effective_num_modes,
            "energy_range": self.effective_energy_range,
            "replicates": self.effective_replicates,
            "seeds": seeds,
            # Spec 037 R3a: where exhaustiveness came from (user_explicit, basic_preset:<name>, advanced_default).
            "exhaustiveness_source": str(self.runtime.get("exhaustiveness_source") or ""),
            # Spec 036 R5c/R5d: where num_modes came from (user_explicit, basic_preset:<name>, advanced_default).
            "num_modes_source": str(self.runtime.get("num_modes_source") or ""),
            "replicate_contract": "specs/036-scientific-consistency/replicate_contract.md",
        }

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
                "effective": self.effective_settings(),
                "jobs": [job.to_dict() for job in jobs],
            }
            manifest_path = self.layout["root"] / "run_manifest.json"
            with open(manifest_path, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2)
            return payload

        manifest_path = self.layout["root"] / "run_manifest.json"
        # An interrupted invocation must not leave yesterday's accepted table
        # looking like today's results. Valid skipped jobs are restored below.
        pd.DataFrame(columns=NORMALIZED_SCORE_COLUMNS).to_csv(self.layout["scores"] / "normalized_scores.csv", index=False)
        for job in jobs:
            if job.status != "planned":
                continue
            cwd_value = Path(str(self.runtime.get("execution_workdir") or self.project_root))
            started = time.monotonic()
            result = execute_job(job.to_dict(), workdir=cwd_value, skip_completed=skip_completed)
            elapsed = round(time.monotonic() - started, 3)
            job.wall_time_s = elapsed
            job.duration_s = elapsed
            job.status, job.returncode, job.error = result["status"], result["returncode"], result["error"]
            if job.status == "completed":
                self._record_returned_poses(job)
            write_json(manifest_path, self._run_payload(jobs))

        payload = self._run_payload(jobs)
        write_json(manifest_path, payload)

        normalized = self.collect_normalized_scores(pairlist_rows)
        if normalized is None or normalized.empty:
            normalized = pd.DataFrame(columns=NORMALIZED_SCORE_COLUMNS)
        normalized.to_csv(self.layout["scores"] / "normalized_scores.csv", index=False)

        return payload

    def _run_payload(self, jobs: List[EngineJobResult]) -> Dict[str, object]:
        return {
            "engine": self.name,
            "runtime": self.runtime,
            "effective": self.effective_settings(),
            "jobs": [job.to_dict() for job in jobs],
        }

    def _record_returned_poses(self, job: EngineJobResult) -> None:
        """Spec 036 R5d: compare poses returned with num_modes and warn when fewer come back."""
        pose_path = Path(job.pose_file)
        if not pose_path.exists():
            job.poses_returned = 0
            job.warnings.append("pose file missing after a successful run")
            return
        returned = self.count_returned_poses(pose_path)
        job.poses_returned = returned
        if returned is None or job.num_modes is None:
            return
        if returned < int(job.num_modes):
            job.warnings.append(
                f"poses_returned={returned} is fewer than num_modes={job.num_modes}"
                + (
                    f" (energy_range={job.energy_range} kcal/mol may exclude the remaining modes)"
                    if job.energy_range is not None
                    else ""
                )
            )

    def count_returned_poses(self, pose_file: Path) -> Optional[int]:
        """Engine-specific pose count. None when the engine has no defined count."""
        return None

    def plan_jobs(
        self,
        pairlist_rows: List[PairlistRow],
        skip_completed: bool = False,
    ) -> List[EngineJobResult]:
        jobs: List[EngineJobResult] = []
        seeds = self.replicate_seeds()
        replicate_count = len(seeds)
        for row in pairlist_rows:
            pair_tag = row.tag
            siblings = []
            for replicate_id in range(1, replicate_count + 1):
                job_tag = replicate_job_tag(pair_tag, replicate_id)
                pose_file = self.layout["poses"] / f"{job_tag}{self.pose_extension}"
                siblings.append(
                    {
                        "replicate_id": replicate_id,
                        "seed": seeds[replicate_id - 1],
                        "output_path": str(pose_file),
                    }
                )
            for replicate_id in range(1, replicate_count + 1):
                seed = seeds[replicate_id - 1]
                job_tag = replicate_job_tag(pair_tag, replicate_id)
                pose_file = self.layout["poses"] / f"{job_tag}{self.pose_extension}"
                log_file = self.layout["logs"] / f"{job_tag}.log"
                command = self.build_command(row, pose_file=pose_file, log_file=log_file, seed=seed)
                inputs = self.input_files(row)
                # Ignore transport/round labels; scientific and executable settings
                # remain part of the protocol. Bytes, not file names, identify inputs.
                protocol = {key: value for key, value in self.runtime.items() if key not in TRANSPORT_RUNTIME_KEYS}
                fingerprint = hashlib.sha256(json.dumps({
                    "contract_version": 1, "engine": self.name, "pair": row.to_dict(), "runtime": protocol,
                    "command": [token.replace(str(self.project_root), "<project>") for token in command],
                    "inputs": [{"name": Path(item["path"]).name, "sha256": item["sha256"]} for item in inputs],
                }, sort_keys=True, allow_nan=False).encode()).hexdigest()
                job = EngineJobResult(
                    tag=job_tag,
                    command=command,
                    pose_file=str(pose_file),
                    log_file=str(log_file),
                    status="planned",
                    pair_tag=pair_tag,
                    replicate_id=replicate_id,
                    seed=seed,
                    replicates=[dict(item) for item in siblings],
                    exhaustiveness=self.effective_exhaustiveness,
                    num_modes=self.effective_num_modes,
                    energy_range=self.effective_energy_range,
                    box_method=row.box_method,
                    ligand_rg_angstrom=row.ligand_rg_angstrom,
                    edge_angstrom=row.edge_angstrom if row.edge_angstrom is not None else row.size_x,
                    warnings=[item for item in str(row.box_warnings or "").split(" | ") if item],
                    engine=self.name,
                    fingerprint=fingerprint,
                    completion_file=str(pose_file) + ".completion.json",
                    input_files=inputs,
                    executables=self.executables(),
                    skip_completed=skip_completed,
                )
                # Reuse requires a valid completion record for this exact protocol and inputs,
                # not merely an existing pose file.
                if skip_completed and completion_valid(job.to_dict()):
                    job.status = "skipped"
                jobs.append(job)
        return jobs

    def input_files(self, row: PairlistRow) -> List[Dict[str, object]]:
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
        # A missing input gets an empty digest; execution then fails it rather than planning crashing.
        return [{"path": str(path), "sha256": file_hash(path) if path.is_file() else ""} for path in paths]

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
        """Expected replicate pose files that pass the job-contract validation (partial output is never scored)."""
        for row in pairlist_rows:
            for replicate_id in range(1, self.effective_replicates + 1):
                path = self.layout["poses"] / f"{replicate_job_tag(row.tag, replicate_id)}{self.pose_extension}"
                try:
                    validate_pose_file(path, self.name)
                except (OSError, ValueError, KeyError, IndexError, StopIteration):
                    continue
                yield path

    @abstractmethod
    def build_command(
        self,
        row: PairlistRow,
        pose_file: Path,
        log_file: Path,
        seed: Optional[int] = None,
    ) -> List[str]:
        raise NotImplementedError

    @abstractmethod
    def collect_normalized_scores(self, pairlist_rows: List[PairlistRow]) -> pd.DataFrame:
        raise NotImplementedError
