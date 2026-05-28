from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple


ExecutionEnvironment = Literal["local_cpu", "local_gpu", "conda_env", "container", "remote_hpc"]
SUPPORTED_EXECUTION_ENVIRONMENTS = {"local_cpu", "local_gpu", "conda_env", "container", "remote_hpc"}


@dataclass
class ExecutionEnvironmentConfig:
    environment: str = "local_cpu"
    workdir: str = ""
    prerequisites_dir: str = ""
    required_files_dir: str = ""
    shared_conda_env: str = ""
    container_image: str = ""
    remote_profile: str = ""

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


def normalize_environment(raw_value: str) -> str:
    token = str(raw_value or "").strip().lower()
    alias_map = {
        "": "local_cpu",
        "cpu": "local_cpu",
        "gpu": "local_gpu",
        "conda": "conda_env",
        "containerized": "container",
        "hpc": "remote_hpc",
        "remote": "remote_hpc",
    }
    resolved = alias_map.get(token, token)
    return resolved if resolved in SUPPORTED_EXECUTION_ENVIRONMENTS else "local_cpu"


def validate_environment_selection(
    engines: List[str],
    config: ExecutionEnvironmentConfig,
    runtime_by_engine: Optional[Dict[str, Dict[str, object]]] = None,
) -> Tuple[List[str], List[str]]:
    errors: List[str] = []
    warnings: List[str] = []

    selected_environment = normalize_environment(config.environment)
    if selected_environment == "local_gpu" and "gnina" not in engines:
        warnings.append("`local_gpu` selected, but GNINA is not in selected engines.")

    if selected_environment == "container":
        runtime = runtime_by_engine or {}
        gnina_image = str((runtime.get("gnina", {}) or {}).get("image") or "").strip()
        container_image = str(config.container_image or "").strip()
        if "gnina" in engines and not (gnina_image or container_image):
            errors.append(
                "Container environment selected with GNINA, but no container image was provided. "
                "Set --gnina-image or --container-image."
            )

    if selected_environment == "conda_env":
        runtime = runtime_by_engine or {}
        requires_conda = any(engine in {"vina", "smina"} for engine in engines)
        has_vina_conda = bool(str((runtime.get("vina", {}) or {}).get("conda_env") or "").strip())
        has_smina_conda = bool(str((runtime.get("smina", {}) or {}).get("conda_env") or "").strip())
        shared_conda = bool(str(config.shared_conda_env or "").strip())
        if requires_conda and not (has_vina_conda or has_smina_conda or shared_conda):
            errors.append(
                "Conda environment mode selected, but no conda env was configured for Vina/Smina. "
                "Set --shared-conda-env or engine-specific conda env flags."
            )

    if selected_environment == "remote_hpc":
        warnings.append(
            "`remote_hpc` in `dock run` performs local planning/execution only. "
            "Use `dock deploy`, `dock sync`, and `dock submit` for HPC submission."
        )

    for label, raw_path in [
        ("execution workdir", config.workdir),
        ("prerequisites dir", config.prerequisites_dir),
        ("required files dir", config.required_files_dir),
    ]:
        path_value = str(raw_path or "").strip()
        if not path_value:
            continue
        path = Path(path_value).expanduser()
        if not path.exists():
            errors.append(f"Configured {label} does not exist: {path}")
        elif not path.is_dir():
            errors.append(f"Configured {label} is not a directory: {path}")

    return errors, warnings


def apply_environment_runtime(
    engines: List[str],
    runtime_by_engine: Dict[str, Dict[str, object]],
    config: ExecutionEnvironmentConfig,
) -> Dict[str, Dict[str, object]]:
    selected_environment = normalize_environment(config.environment)
    merged: Dict[str, Dict[str, object]] = {}
    for engine, runtime in runtime_by_engine.items():
        payload = dict(runtime or {})
        payload["execution_environment"] = selected_environment
        if config.workdir:
            payload["execution_workdir"] = str(Path(config.workdir).expanduser().resolve())
        if config.prerequisites_dir:
            payload["prerequisites_dir"] = str(Path(config.prerequisites_dir).expanduser().resolve())
        if config.required_files_dir:
            payload["required_files_dir"] = str(Path(config.required_files_dir).expanduser().resolve())
        merged[engine] = payload

    if selected_environment == "local_cpu":
        if "gnina" in engines and "gnina" in merged:
            gnina_runtime = dict(merged["gnina"])
            gnina_runtime["use_gpu"] = False
            merged["gnina"] = gnina_runtime

    if selected_environment == "local_gpu":
        if "gnina" in engines and "gnina" in merged:
            gnina_runtime = dict(merged["gnina"])
            if gnina_runtime.get("use_gpu") in (None, ""):
                gnina_runtime["use_gpu"] = True
            merged["gnina"] = gnina_runtime

    if selected_environment == "conda_env":
        shared_env = str(config.shared_conda_env or "").strip()
        if shared_env:
            for engine in ("vina", "smina"):
                if engine in merged:
                    payload = dict(merged[engine])
                    if not str(payload.get("conda_env") or "").strip():
                        payload["conda_env"] = shared_env
                    merged[engine] = payload

    if selected_environment == "container":
        container_image = str(config.container_image or "").strip()
        if container_image and "gnina" in merged:
            payload = dict(merged["gnina"])
            if not str(payload.get("image") or "").strip():
                payload["image"] = container_image
            merged["gnina"] = payload

    return merged
