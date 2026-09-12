"""Portable, standard-library job execution and output validation.

Deployment bundles copy this module verbatim so local and scheduled jobs use
the same acceptance rules without requiring the analysis stack on workers.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import uuid


def file_hash(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Keep temporary names short: repeating a long pose name and UUID can
    # exceed Windows MAX_PATH even when the published file itself fits.
    descriptor, temporary = tempfile.mkstemp(prefix=".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, allow_nan=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _finite_numbers(tokens):
    return bool(tokens) and all(math.isfinite(float(value)) for value in tokens)


def parse_sdf_scores(path):
    records = []
    text = Path(path).read_text(encoding="utf-8", errors="strict")
    if not text.rstrip().endswith("$$$$"):
        raise ValueError("SDF is incomplete (missing final record delimiter)")
    blocks = text.split("$$$$")[:-1]
    for number, block in enumerate(blocks, 1):
        # The delimiter is followed by a newline before the next title line.
        if number > 1:
            block = block.removeprefix("\r\n").removeprefix("\n")
        lines = block.splitlines()
        counts_index = next((i for i, line in enumerate(lines[:5]) if "V2000" in line or "V3000" in line), None)
        if counts_index is None or not any(line.strip() == "M  END" for line in lines):
            raise ValueError(f"SDF record {number} has no valid structure header/end")
        if "V3000" in lines[counts_index]:
            count_line = next(line for line in lines if line.startswith("M  V30 COUNTS "))
            atom_count = int(count_line.split()[3])
            start = lines.index("M  V30 BEGIN ATOM") + 1
            stop = lines.index("M  V30 END ATOM")
            atoms = lines[start:stop]
            coordinates = [line.split()[4:7] for line in atoms]
        else:
            atom_count = int(lines[counts_index][:3])
            atoms = lines[counts_index + 1:counts_index + 1 + atom_count]
            coordinates = [[line[:10], line[10:20], line[20:30]] for line in atoms]
        if atom_count < 1 or len(atoms) != atom_count or any(len(xyz) != 3 or not _finite_numbers(xyz) for xyz in coordinates):
            raise ValueError(f"SDF record {number} has invalid atoms/coordinates")
        properties = {}
        for index, line in enumerate(lines):
            match = re.match(r">\s*<([^>]+)>.*", line)
            if match and index + 1 < len(lines):
                properties[match.group(1)] = lines[index + 1].strip()
        affinity = float(properties["minimizedAffinity"])
        if not math.isfinite(affinity):
            raise ValueError("Nonfinite SDF affinity")
        record = {"pose": number, "affinity": affinity}
        for prop, name in [("CNNaffinity", "cnn_affinity"), ("CNNscore", "cnn_score")]:
            value = float(properties[prop]) if prop in properties else None
            if value is not None and not math.isfinite(value):
                raise ValueError(f"Nonfinite SDF {prop}")
            record[name] = value
        records.append(record)
    if not records:
        raise ValueError("No SDF poses")
    return records


def validate_pose_file(path, engine):
    """Reject missing, partial, unscored and coordinate-invalid results."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("Pose file is missing or empty")
    if engine == "gnina":
        return len(parse_sdf_scores(path))
    text = path.read_text(encoding="utf-8", errors="strict")
    if engine == "autodock4":
        if "Successful Completion" not in text:
            raise ValueError("AutoDock4 did not report Successful Completion")
        lines = [line[len("DOCKED: "):] for line in text.splitlines() if line.startswith("DOCKED: ")]
    else:
        lines = text.splitlines()
    models = []
    current = []
    in_model = False
    for line in lines:
        if line.startswith("MODEL "):
            if in_model:
                raise ValueError("Nested/incomplete pose MODEL")
            current = []
            in_model = True
        elif line.startswith("ENDMDL"):
            if not in_model:
                raise ValueError("ENDMDL without MODEL")
            models.append(current)
            current = []
            in_model = False
        else:
            current.append(line)
    if in_model:
        raise ValueError("Incomplete pose MODEL")
    if not models:
        models = [current]
    for model in models:
        atoms = [line for line in model if line.startswith(("ATOM  ", "HETATM"))]
        if not atoms or any(not _finite_numbers([line[30:38], line[38:46], line[46:54]]) for line in atoms):
            raise ValueError("Pose has no atoms or invalid coordinates")
        scores = []
        for line in model:
            if line.startswith("REMARK VINA RESULT:"):
                scores.append(float(line.split()[3]))
            elif line.startswith("REMARK minimizedAffinity"):
                scores.append(float(line.split()[2]))
            elif engine == "autodock4" and "Estimated Free Energy of Binding" in line:
                scores.append(float(line.split("=", 1)[1].split()[0]))
        if not scores or not all(math.isfinite(value) for value in scores):
            raise ValueError("Pose has no finite affinity")
    return len(models)


def execution_identity(job):
    """Capture actual engine executables, rather than a conda/bash wrapper alone."""
    identities = {}
    for executable in job.get("executables", []):
        path = shutil.which(executable)
        if path is None and Path(executable).is_file():
            path = executable
        if path is None:
            # A partial identity (e.g. AutoDock found but AutoGrid missing on
            # this PATH) cannot certify the full scientific computation.
            return {}
        identities[executable] = file_hash(path)
    return identities


def completion_valid(job):
    try:
        record = json.loads(Path(job["completion_file"]).read_text(encoding="utf-8"))
        identity = execution_identity(job)
        return bool(identity) and (
            record.get("status") == "completed"
            and record.get("fingerprint") == job["fingerprint"]
            and record.get("execution_identity") == identity
            and record.get("pose_sha256") == file_hash(job["pose_file"])
            and all(file_hash(item["path"]) == item["sha256"] for item in job.get("input_files", []))
            and validate_pose_file(job["pose_file"], job["engine"]) > 0
        )
    except (OSError, ValueError, KeyError, IndexError, StopIteration):
        return False


def _execute_job(job, workdir=None, skip_completed=False):
    job = dict(job)
    if skip_completed and completion_valid(job):
        job.update(status="skipped", returncode=0, error="")
        return job
    pose = Path(job["pose_file"])
    pose.parent.mkdir(parents=True, exist_ok=True)
    attempt = pose.parent / ".attempts" / uuid.uuid4().hex
    attempt.mkdir(parents=True)
    # Retire old published artifacts before starting so a failure cannot expose
    # a previous pose as the result of this attempt. Originals remain archived.
    for raw in [job["pose_file"], job["log_file"], job["completion_file"]]:
        old = Path(raw)
        if old.exists():
            shutil.move(str(old), str(attempt / old.name))
    job.update(status="running", returncode=None, error="", execution_identity=execution_identity(job))
    write_json(job["completion_file"], job)
    try:
        for item in job.get("input_files", []):
            if file_hash(item["path"]) != item["sha256"]:
                raise ValueError(f"Input changed after planning: {item['path']}")
        completed = subprocess.run(job["command"], cwd=workdir, capture_output=True, text=True)
        log_text = f"COMMAND: {job['command']!r}\n\nSTDOUT:\n{completed.stdout}\n\nSTDERR:\n{completed.stderr}"
        pair_log = Path(job["log_file"])
        pair_log.parent.mkdir(parents=True, exist_ok=True)
        if not pair_log.exists() or pair_log.stat().st_size == 0:
            pair_log.write_text(log_text, encoding="utf-8")
        (attempt / "runner.log").write_text(log_text, encoding="utf-8")
        job["returncode"] = completed.returncode
        if completed.returncode != 0:
            raise ValueError(completed.stderr.strip() or f"Engine returned {completed.returncode}")
        for item in job.get("input_files", []):
            if file_hash(item["path"]) != item["sha256"]:
                raise ValueError(f"Input changed during execution: {item['path']}")
        if execution_identity(job) != job["execution_identity"]:
            raise ValueError("Engine executable changed during execution")
        job["pose_count"] = validate_pose_file(pose, job["engine"])
        job["pose_sha256"] = file_hash(pose)
        job["status"] = "completed"
    except (OSError, ValueError, KeyError, IndexError, StopIteration) as exc:
        job.update(status="failed", error=str(exc))
        if pose.exists():
            shutil.move(str(pose), str(attempt / ("invalid_" + pose.name)))
    write_json(job["completion_file"], job)
    write_json(attempt / "result.json", job)
    return job


def execute_job(job, workdir=None, skip_completed=False):
    """Serialize attempts for one output, including independent scheduler jobs."""
    lock = Path(job["completion_file"] + ".lock")
    lock.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except OSError as exc:
        return {**job, "status": "failed", "returncode": None, "error": f"Cannot acquire job lock {lock}: {exc}. Check whether another attempt is active before clearing a stale lock."}
    try:
        with os.fdopen(descriptor, "w") as handle:
            json.dump({"pid": os.getpid(), "host": __import__("socket").gethostname()}, handle)
        return _execute_job(job, workdir, skip_completed)
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest")
    parser.add_argument("index", type=int)
    args = parser.parse_args()
    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    job = manifest["jobs"][args.index]
    result = execute_job(job, workdir=manifest["execution_project_root"], skip_completed=job.get("skip_completed", False))
    if result["status"] not in {"completed", "skipped"}:
        print(result["error"], file=__import__("sys").stderr)
        raise SystemExit(1)
