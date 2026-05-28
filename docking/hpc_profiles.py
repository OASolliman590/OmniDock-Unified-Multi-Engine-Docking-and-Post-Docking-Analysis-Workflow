from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from typing import Dict, List, Optional

import yaml


BUILTIN_HPC_PROFILES: Dict[str, Dict[str, object]] = {
    "nmrbox-condor": {
        "name": "nmrbox-condor",
        "scheduler": "condor",
        "description": (
            "Safe public template for NMRbox HTCondor cluster. "
            "Copy to .workflow/hpc_profiles/nmrbox-condor.json and fill in local values."
        ),
        "shell_preamble": [
            "if [ -f ~/.bashrc ]; then source ~/.bashrc; fi",
        ],
        "condor": {
            "job_name_prefix": "dock",
        },
        "remote": {
            "ssh_target": "<user>@sulfur.nmrbox.org",
            "project_root_base": "/home/<user>/docking",
            "sync_excludes": [
                ".DS_Store",
                "__pycache__/",
                ".workflow/state.json",
            ],
        },
        "engines": {
            "gnina": {
                "runtime": {
                    "binary": "/home/<user>/gnina",
                    "image": "/home/<user>/gnina.sif",
                    "use_gpu": True,
                },
                "preamble": [],
                "condor_gpu": {
                    "submit_mode": "condor_array",
                    "cpus": 16,
                    "memory": "48GB",
                    "disk": "20GB",
                    "gpus": 1,
                    "parallelism": 1,
                },
                "condor_cpu": {
                    "submit_mode": "single_job",
                    "cpus": 16,
                    "memory": "48GB",
                    "disk": "20GB",
                    "gpus": 0,
                    "parallelism": 0,
                },
            },
            "vina": {
                "runtime": {
                    "binary": "~/.local/bin/vina",
                },
                "condor": {
                    "submit_mode": "condor_array",
                    "cpus": 8,
                    "memory": "24GB",
                    "disk": "10GB",
                    "gpus": 0,
                    "parallelism": 4,
                },
            },
            "smina": {
                "runtime": {
                    "binary": "~/.local/bin/smina",
                },
                "condor": {
                    "submit_mode": "condor_array",
                    "cpus": 8,
                    "memory": "24GB",
                    "disk": "10GB",
                    "gpus": 0,
                    "parallelism": 4,
                },
            },
        },
    },
    "generic-apptainer-conda": {
        "name": "generic-apptainer-conda",
        "description": (
            "Safe generic template for clusters that run GNINA through Apptainer "
            "and CPU-bound engines through Conda-provided CLIs."
        ),
        "shell_preamble": [
            "if [ -f ~/.bashrc ]; then source ~/.bashrc; fi",
        ],
        "slurm": {},
        "remote": {},
        "engines": {
            "gnina": {
                "preamble": [
                    "if command -v module >/dev/null 2>&1; then module load apptainer >/dev/null 2>&1 || true; fi",
                ],
                "slurm_gpu": {
                    "time": "72:00:00",
                    "mem": "48G",
                    "cpus_per_task": 16,
                    "gpus": 1,
                    "array_parallelism": 1,
                },
                "slurm_cpu": {
                    "time": "72:00:00",
                    "mem": "48G",
                    "cpus_per_task": 16,
                    "gpus": 0,
                    "array_parallelism": 2,
                },
            },
            "vina": {
                "slurm": {
                    "time": "24:00:00",
                    "mem": "24G",
                    "cpus_per_task": 16,
                    "gpus": 0,
                    "array_parallelism": 4,
                }
            },
            "smina": {
                "slurm": {
                    "time": "24:00:00",
                    "mem": "24G",
                    "cpus_per_task": 16,
                    "gpus": 0,
                    "array_parallelism": 4,
                }
            },
            "autodock4": {
                "slurm": {
                    "time": "24:00:00",
                    "mem": "24G",
                    "cpus_per_task": 16,
                    "gpus": 0,
                    "array_parallelism": 4,
                }
            },
        },
    },
    "bibalex-apptainer-conda": {
        "name": "bibalex-apptainer-conda",
        "description": (
            "Safe public template for the BibAlex-style deployment pattern. "
            "Real account names and home-directory paths belong in a local ignored profile file."
        ),
        "shell_preamble": [
            "if [ -f ~/.bashrc ]; then source ~/.bashrc; fi",
        ],
        "slurm": {},
        "remote": {
            "ssh_target": "<user>@<hpc-host>",
            "project_root_base": "/cluster/users/<user>/docking",
            "sync_excludes": [
                ".DS_Store",
                "__pycache__/",
                ".workflow/state.json",
            ],
        },
        "engines": {
            "gnina": {
                "preamble": [
                    "if command -v module >/dev/null 2>&1; then module load apptainer >/dev/null 2>&1 || true; fi",
                ],
                "runtime": {
                    "binary": "<set-gnina-binary>",
                    "image": "<set-gnina-image>",
                    "use_gpu": False,
                },
                "slurm_gpu": {
                    "submit_mode": "single_job",
                    "partition": "gpu",
                    "time": "72:00:00",
                    "mem": "48G",
                    "cpus_per_task": 16,
                    "gpus": 1,
                    "array_parallelism": 1,
                },
                "slurm_cpu": {
                    "submit_mode": "single_job",
                    "partition": "cpu",
                    "time": "72:00:00",
                    "mem": "48G",
                    "cpus_per_task": 16,
                    "gpus": 0,
                    "array_parallelism": 2,
                },
            },
            "vina": {
                "runtime": {
                    "conda_env": "<set-cpu-conda-env>",
                    "binary": "vina",
                },
                "slurm": {
                    "submit_mode": "single_job",
                    "partition": "cpu",
                    "time": "24:00:00",
                    "mem": "24G",
                    "cpus_per_task": 16,
                    "gpus": 0,
                    "array_parallelism": 4,
                },
            },
            "smina": {
                "runtime": {
                    "conda_env": "<set-cpu-conda-env>",
                    "binary": "smina",
                },
                "slurm": {
                    "submit_mode": "single_job",
                    "partition": "cpu",
                    "time": "24:00:00",
                    "mem": "24G",
                    "cpus_per_task": 16,
                    "gpus": 0,
                    "array_parallelism": 4,
                },
            },
            "autodock4": {
                "runtime": {
                    "binary": "autodock4",
                    "autogrid_binary": "autogrid4",
                    "parameter_file": "AD4.1_bound.dat",
                    "parameter_file_path": "<set-absolute-AD4.1_bound.dat-path>",
                    "autodocktools_python": "<set-mgltools-pythonsh>",
                    "autodocktools_prepare_gpf4": "<set-prepare_gpf4.py>",
                    "autodocktools_prepare_dpf4": "<set-prepare_dpf4.py>",
                },
                "slurm": {
                    "submit_mode": "single_job",
                    "partition": "cpu",
                    "time": "24:00:00",
                    "mem": "24G",
                    "cpus_per_task": 16,
                    "gpus": 0,
                    "array_parallelism": 4,
                },
            },
        },
    },
}


def _coerce_lines(raw_value: object) -> List[str]:
    if raw_value in (None, ""):
        return []
    if isinstance(raw_value, list):
        return [str(item).strip() for item in raw_value if str(item).strip()]
    return [str(raw_value).strip()]


def _coerce_mapping(raw_value: object) -> Dict[str, object]:
    if isinstance(raw_value, dict):
        return dict(raw_value)
    return {}


def _gnina_uses_gpu(runtime: Dict[str, object], engine_profile: Optional[Dict[str, object]] = None) -> bool:
    explicit = runtime.get("use_gpu")
    if explicit not in (None, ""):
        return bool(explicit)
    if engine_profile:
        base_runtime = _coerce_mapping(engine_profile.get("runtime"))
        if base_runtime.get("use_gpu") not in (None, ""):
            return bool(base_runtime.get("use_gpu"))
    return True


def _local_profile_candidates(project_root: Path, profile_name: str) -> List[Path]:
    root = Path(project_root).expanduser().resolve()
    def _profile_paths(base: Path) -> List[Path]:
        return [
            base / ".workflow" / "hpc_profiles" / f"{profile_name}.json",
            base / ".workflow" / "hpc_profiles" / f"{profile_name}.yaml",
            base / ".workflow" / "hpc_profiles" / f"{profile_name}.yml",
            base / ".codex" / "hpc_profiles" / f"{profile_name}.json",
            base / ".codex" / "hpc_profiles" / f"{profile_name}.yaml",
            base / ".codex" / "hpc_profiles" / f"{profile_name}.yml",
        ]

    candidates = _profile_paths(root)
    cwd_root = Path.cwd().resolve()
    if cwd_root not in {root}:
        candidates.extend(_profile_paths(cwd_root))
    candidates.extend(
        [
            Path.home() / ".config" / "pdb-prepare-wizard" / "hpc_profiles" / f"{profile_name}.json",
            Path.home() / ".config" / "pdb-prepare-wizard" / "hpc_profiles" / f"{profile_name}.yaml",
            Path.home() / ".config" / "pdb-prepare-wizard" / "hpc_profiles" / f"{profile_name}.yml",
        ]
    )
    return candidates


def _load_profile_file(path: Path) -> Dict[str, object]:
    suffix = str(path.suffix or "").strip().lower()
    with open(path, "r", encoding="utf-8") as handle:
        if suffix in {".yaml", ".yml"}:
            payload = yaml.safe_load(handle)
        else:
            payload = json.load(handle)
    if payload is None:
        payload = {}
    if not isinstance(payload, dict):
        raise ValueError(f"HPC profile must be a JSON/YAML object: {path}")
    return payload


def load_hpc_profile(
    project_root: Path,
    profile_name: str = "",
    profile_file: str = "",
) -> Dict[str, object]:
    if profile_file:
        path = Path(profile_file).expanduser().resolve()
        payload = _load_profile_file(path)
        payload.setdefault("name", path.stem)
        payload["source"] = "profile-file"
        return normalize_hpc_profile(payload)

    normalized_name = str(profile_name or "").strip()
    if not normalized_name:
        return {}

    for path in _local_profile_candidates(project_root, normalized_name):
        if path.exists():
            payload = _load_profile_file(path)
            payload.setdefault("name", normalized_name)
            payload["source"] = "local-file"
            return normalize_hpc_profile(payload)

    if normalized_name in BUILTIN_HPC_PROFILES:
        payload = deepcopy(BUILTIN_HPC_PROFILES[normalized_name])
        payload["source"] = f"builtin:{normalized_name}"
        return normalize_hpc_profile(payload)

    locations = ", ".join(str(path) for path in _local_profile_candidates(project_root, normalized_name))
    raise ValueError(
        f"Unknown HPC profile '{normalized_name}'. Expected a builtin profile or a YAML/JSON file at one of: {locations}"
    )


def normalize_hpc_profile(payload: Dict[str, object]) -> Dict[str, object]:
    profile = deepcopy(payload)
    profile.setdefault("name", "custom")
    profile.setdefault("source", "")
    raw_scheduler = str(profile.get("scheduler") or "slurm").strip().lower()
    profile["scheduler"] = raw_scheduler if raw_scheduler in {"slurm", "condor"} else "slurm"
    profile["shell_preamble"] = _coerce_lines(profile.get("shell_preamble"))
    profile["slurm"] = _coerce_mapping(profile.get("slurm"))
    profile["condor"] = _coerce_mapping(profile.get("condor"))
    remote = _coerce_mapping(profile.get("remote"))
    profile["remote"] = {
        "ssh_target": str(remote.get("ssh_target", "") or "").strip(),
        "project_root": str(remote.get("project_root", "") or "").strip(),
        "project_root_base": str(remote.get("project_root_base", "") or "").strip(),
        "sync_excludes": _coerce_lines(remote.get("sync_excludes")),
    }

    normalized_engines: Dict[str, Dict[str, object]] = {}
    for engine_name, engine_payload in _coerce_mapping(profile.get("engines")).items():
        engine_cfg = _coerce_mapping(engine_payload)
        normalized_engines[str(engine_name)] = {
            "runtime": _coerce_mapping(engine_cfg.get("runtime")),
            "runtime_cpu": _coerce_mapping(engine_cfg.get("runtime_cpu")),
            "runtime_gpu": _coerce_mapping(engine_cfg.get("runtime_gpu")),
            "slurm": _coerce_mapping(engine_cfg.get("slurm")),
            "slurm_cpu": _coerce_mapping(engine_cfg.get("slurm_cpu")),
            "slurm_gpu": _coerce_mapping(engine_cfg.get("slurm_gpu")),
            "condor": _coerce_mapping(engine_cfg.get("condor")),
            "condor_cpu": _coerce_mapping(engine_cfg.get("condor_cpu")),
            "condor_gpu": _coerce_mapping(engine_cfg.get("condor_gpu")),
            "preamble": _coerce_lines(engine_cfg.get("preamble")),
        }
    profile["engines"] = normalized_engines
    return profile


def resolve_remote_project_dir(
    local_project_root: Path,
    profile: Dict[str, object],
    override: str = "",
) -> str:
    explicit = str(override or "").strip()
    if explicit:
        return explicit
    remote = _coerce_mapping(profile.get("remote"))
    project_root = str(remote.get("project_root", "") or "").strip()
    if project_root:
        return project_root
    project_root_base = str(remote.get("project_root_base", "") or "").strip()
    if project_root_base:
        return str(Path(project_root_base) / Path(local_project_root).expanduser().resolve().name)
    return ""


def resolve_ssh_target(profile: Dict[str, object], override: str = "") -> str:
    explicit = str(override or "").strip()
    if explicit:
        return explicit
    remote = _coerce_mapping(profile.get("remote"))
    return str(remote.get("ssh_target", "") or "").strip()


def merge_runtime_with_profile(
    runtime_by_engine: Dict[str, Dict[str, object]],
    profile: Dict[str, object],
) -> Dict[str, Dict[str, object]]:
    if not profile:
        return runtime_by_engine

    merged_runtime: Dict[str, Dict[str, object]] = {}
    global_preamble = _coerce_lines(profile.get("shell_preamble"))
    engine_profiles = _coerce_mapping(profile.get("engines"))

    for engine, runtime in runtime_by_engine.items():
        engine_profile = _coerce_mapping(engine_profiles.get(engine))
        engine_runtime = dict(_coerce_mapping(engine_profile.get("runtime")))
        variant_key = "runtime_gpu" if engine == "gnina" and _gnina_uses_gpu(runtime, engine_profile) else "runtime_cpu"
        engine_runtime.update(_coerce_mapping(engine_profile.get(variant_key)))
        for key, value in runtime.items():
            if value not in (None, ""):
                engine_runtime[key] = value
            else:
                engine_runtime.setdefault(key, value)
        combined_preamble = global_preamble + _coerce_lines(engine_profile.get("preamble"))
        if combined_preamble:
            existing = _coerce_lines(engine_runtime.get("deployment_preamble"))
            engine_runtime["deployment_preamble"] = combined_preamble + existing
        merged_runtime[engine] = engine_runtime
    return merged_runtime


def merge_slurm_with_profile(
    engine: str,
    runtime: Dict[str, object],
    cli_slurm: Dict[str, object],
    profile: Dict[str, object],
) -> Dict[str, object]:
    merged = {}
    if profile:
        merged.update(_coerce_mapping(profile.get("slurm")))
        engine_profile = _coerce_mapping(_coerce_mapping(profile.get("engines")).get(engine))
        merged.update(_coerce_mapping(engine_profile.get("slurm")))
        variant_key = "slurm_gpu" if engine == "gnina" and _gnina_uses_gpu(runtime, engine_profile) else "slurm_cpu"
        merged.update(_coerce_mapping(engine_profile.get(variant_key)))

        extra_args = []
        extra_args.extend(_coerce_lines(_coerce_mapping(profile.get("slurm")).get("extra_args")))
        extra_args.extend(_coerce_lines(_coerce_mapping(engine_profile.get("slurm")).get("extra_args")))
        extra_args.extend(_coerce_lines(_coerce_mapping(engine_profile.get(variant_key)).get("extra_args")))
        if extra_args:
            merged["extra_args"] = extra_args

    for key, value in cli_slurm.items():
        if value in (None, ""):
            continue
        if key == "extra_args":
            merged["extra_args"] = _coerce_lines(merged.get("extra_args")) + _coerce_lines(value)
        else:
            merged[key] = value
    return merged


def merge_condor_with_profile(
    engine: str,
    runtime: Dict[str, object],
    cli_condor: Dict[str, object],
    profile: Dict[str, object],
) -> Dict[str, object]:
    merged = {}
    if profile:
        merged.update(_coerce_mapping(profile.get("condor")))
        engine_profile = _coerce_mapping(_coerce_mapping(profile.get("engines")).get(engine))
        merged.update(_coerce_mapping(engine_profile.get("condor")))
        variant_key = "condor_gpu" if engine == "gnina" and _gnina_uses_gpu(runtime, engine_profile) else "condor_cpu"
        merged.update(_coerce_mapping(engine_profile.get(variant_key)))

    for key, value in cli_condor.items():
        if value in (None, ""):
            continue
        merged[key] = value
    return merged
