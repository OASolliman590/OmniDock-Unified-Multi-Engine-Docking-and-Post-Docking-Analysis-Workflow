"""SQLite storage backend for DockForge comparative/post-docking outputs.

This module is intentionally feature-flag friendly:
- CSV remains the source-of-truth default.
- SQLite dual-write can be enabled to persist tabular outputs for query workflows.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

import pandas as pd


SQLITE_PRAGMAS = (
    "PRAGMA journal_mode=WAL;",
    "PRAGMA synchronous=NORMAL;",
    "PRAGMA foreign_keys=ON;",
)


@dataclass
class SQLiteWriteSummary:
    database_path: str
    run_id: str
    datasets_written: int
    row_counts: Dict[str, int]
    started_at: str
    completed_at: str

    def to_dict(self) -> Dict[str, object]:
        return {
            "database_path": self.database_path,
            "run_id": self.run_id,
            "datasets_written": int(self.datasets_written),
            "row_counts": {str(key): int(value) for key, value in self.row_counts.items()},
            "started_at": self.started_at,
            "completed_at": self.completed_at,
        }


@dataclass
class ParityResult:
    dataset_key: str
    csv_path: str
    csv_exists: bool
    csv_rows: int
    sqlite_rows: int
    column_match: bool
    status: str
    note: str

    def to_dict(self) -> Dict[str, object]:
        return {
            "dataset_key": self.dataset_key,
            "csv_path": self.csv_path,
            "csv_exists": bool(self.csv_exists),
            "csv_rows": int(self.csv_rows),
            "sqlite_rows": int(self.sqlite_rows),
            "column_match": bool(self.column_match),
            "status": self.status,
            "note": self.note,
        }


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sanitize_identifier(value: str) -> str:
    token = "".join(ch if (ch.isalnum() or ch == "_") else "_" for ch in str(value or "").strip().lower())
    token = token.strip("_")
    if not token:
        token = "unnamed"
    if token[0].isdigit():
        token = f"x_{token}"
    return token


def _json_dumps(value: object) -> str:
    try:
        return json.dumps(value, sort_keys=True)
    except Exception:
        return json.dumps(str(value))


def _coerce_sql_value(value: object) -> Optional[str]:
    if value is None:
        return None
    if pd.isna(value):
        return None
    if isinstance(value, (dict, list, tuple, set)):
        return _json_dumps(value)
    if isinstance(value, (datetime,)):
        return value.isoformat()
    return str(value)


def _connect(db_path: Path) -> sqlite3.Connection:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    for pragma in SQLITE_PRAGMAS:
        conn.execute(pragma)
    return conn


def initialize_schema(db_path: Path) -> Path:
    target = Path(db_path).expanduser().resolve()
    with _connect(target) as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS runs (
                run_id TEXT PRIMARY KEY,
                generated_at TEXT NOT NULL,
                project_dir TEXT,
                output_dir TEXT,
                analysis_mode TEXT,
                analysis_scope TEXT,
                consensus_mode TEXT,
                normalization_method TEXT,
                metadata_json TEXT,
                created_at TEXT NOT NULL
            )
            """
        )
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS dataset_registry (
                run_id TEXT NOT NULL,
                dataset_key TEXT NOT NULL,
                table_name TEXT NOT NULL,
                row_count INTEGER NOT NULL,
                columns_json TEXT NOT NULL,
                source_csv TEXT,
                updated_at TEXT NOT NULL,
                PRIMARY KEY (run_id, dataset_key),
                FOREIGN KEY (run_id) REFERENCES runs(run_id) ON DELETE CASCADE
            )
            """
        )
        conn.commit()
    return target


def _ensure_dataset_table(conn: sqlite3.Connection, table_name: str, columns: Iterable[str]) -> List[str]:
    normalized = [_sanitize_identifier(column) for column in columns]
    normalized = [column if column != "run_id" else "csv_run_id" for column in normalized]

    conn.execute(
        f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            run_id TEXT NOT NULL,
            row_index INTEGER NOT NULL,
            PRIMARY KEY (run_id, row_index)
        )
        """
    )

    existing_cols = {
        str(row[1])
        for row in conn.execute(f"PRAGMA table_info({table_name})").fetchall()
        if len(row) >= 2
    }
    for column in normalized:
        if column in {"run_id", "row_index"}:
            continue
        if column not in existing_cols:
            conn.execute(f"ALTER TABLE {table_name} ADD COLUMN {column} TEXT")
    return normalized


def _insert_dataset_rows(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    dataset_key: str,
    table_name: str,
    frame: pd.DataFrame,
    source_csv: str = "",
) -> int:
    safe_frame = frame.copy()
    normalized_columns = _ensure_dataset_table(conn, table_name, list(safe_frame.columns))

    conn.execute(f"DELETE FROM {table_name} WHERE run_id = ?", (run_id,))

    insert_cols = ["run_id", "row_index"] + [
        column for column in normalized_columns if column not in {"run_id", "row_index"}
    ]
    placeholders = ", ".join(["?"] * len(insert_cols))
    insert_sql = f"INSERT INTO {table_name} ({', '.join(insert_cols)}) VALUES ({placeholders})"

    rows: List[Tuple[object, ...]] = []
    for row_idx, row in enumerate(safe_frame.itertuples(index=False), start=1):
        values = [run_id, int(row_idx)]
        for original_idx, original_column in enumerate(safe_frame.columns):
            normalized_column = normalized_columns[original_idx]
            if normalized_column in {"run_id", "row_index"}:
                continue
            values.append(_coerce_sql_value(row[original_idx]))
        rows.append(tuple(values))

    if rows:
        conn.executemany(insert_sql, rows)

    registry_payload = {
        "run_id": run_id,
        "dataset_key": dataset_key,
        "table_name": table_name,
        "row_count": int(len(safe_frame)),
        "columns_json": _json_dumps([str(column) for column in safe_frame.columns]),
        "source_csv": str(source_csv or ""),
        "updated_at": _now_iso(),
    }
    conn.execute(
        """
        INSERT INTO dataset_registry (
            run_id, dataset_key, table_name, row_count, columns_json, source_csv, updated_at
        ) VALUES (
            :run_id, :dataset_key, :table_name, :row_count, :columns_json, :source_csv, :updated_at
        )
        ON CONFLICT(run_id, dataset_key)
        DO UPDATE SET
            table_name=excluded.table_name,
            row_count=excluded.row_count,
            columns_json=excluded.columns_json,
            source_csv=excluded.source_csv,
            updated_at=excluded.updated_at
        """,
        registry_payload,
    )
    return int(len(safe_frame))


def write_comparative_bundle(
    db_path: Path,
    *,
    run_id: str,
    run_metadata: Optional[Mapping[str, object]] = None,
    dataset_frames: Mapping[str, pd.DataFrame],
    source_csv_map: Optional[Mapping[str, str]] = None,
) -> SQLiteWriteSummary:
    """Persist comparative outputs into SQLite tables keyed by run_id."""

    started_at = _now_iso()
    target = initialize_schema(Path(db_path))
    row_counts: Dict[str, int] = {}

    run_meta = dict(run_metadata or {})
    with _connect(target) as conn:
        conn.execute(
            """
            INSERT INTO runs (
                run_id, generated_at, project_dir, output_dir, analysis_mode, analysis_scope,
                consensus_mode, normalization_method, metadata_json, created_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(run_id)
            DO UPDATE SET
                generated_at=excluded.generated_at,
                project_dir=excluded.project_dir,
                output_dir=excluded.output_dir,
                analysis_mode=excluded.analysis_mode,
                analysis_scope=excluded.analysis_scope,
                consensus_mode=excluded.consensus_mode,
                normalization_method=excluded.normalization_method,
                metadata_json=excluded.metadata_json
            """,
            (
                str(run_id),
                str(run_meta.get("generated_at") or _now_iso()),
                str(run_meta.get("project_dir") or ""),
                str(run_meta.get("output_dir") or ""),
                str(run_meta.get("analysis_mode") or ""),
                str(run_meta.get("analysis_scope") or ""),
                str(run_meta.get("consensus_mode") or ""),
                str(run_meta.get("normalization_method") or ""),
                _json_dumps(run_meta),
                _now_iso(),
            ),
        )

        for dataset_key, frame in dataset_frames.items():
            table_name = f"dataset_{_sanitize_identifier(dataset_key)}"
            csv_source = str((source_csv_map or {}).get(dataset_key, ""))
            row_count = _insert_dataset_rows(
                conn,
                run_id=str(run_id),
                dataset_key=str(dataset_key),
                table_name=table_name,
                frame=frame if isinstance(frame, pd.DataFrame) else pd.DataFrame(),
                source_csv=csv_source,
            )
            row_counts[str(dataset_key)] = int(row_count)

        conn.commit()

    completed_at = _now_iso()
    return SQLiteWriteSummary(
        database_path=str(target),
        run_id=str(run_id),
        datasets_written=len(dataset_frames),
        row_counts=row_counts,
        started_at=started_at,
        completed_at=completed_at,
    )


def read_dataset_registry(db_path: Path, run_id: str) -> pd.DataFrame:
    target = Path(db_path).expanduser().resolve()
    if not target.exists():
        return pd.DataFrame(
            columns=["run_id", "dataset_key", "table_name", "row_count", "columns_json", "source_csv", "updated_at"]
        )
    with _connect(target) as conn:
        return pd.read_sql_query(
            "SELECT run_id, dataset_key, table_name, row_count, columns_json, source_csv, updated_at "
            "FROM dataset_registry WHERE run_id = ? ORDER BY dataset_key",
            conn,
            params=(str(run_id),),
        )


def _sha256_text(payload: str) -> str:
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def validate_csv_sqlite_parity(
    db_path: Path,
    *,
    run_id: str,
    dataset_csv_map: Mapping[str, Path],
) -> pd.DataFrame:
    """Validate row-count and column-parity between CSV outputs and SQLite registry."""

    registry = read_dataset_registry(db_path, run_id)
    registry_index = {
        str(row["dataset_key"]): row for _, row in registry.iterrows()
    }

    rows: List[Dict[str, object]] = []
    for dataset_key, csv_path_raw in dataset_csv_map.items():
        csv_path = Path(csv_path_raw).expanduser().resolve()
        reg = registry_index.get(str(dataset_key))

        csv_exists = csv_path.exists() and csv_path.is_file()
        if csv_exists:
            csv_df = pd.read_csv(csv_path)
            csv_rows = int(len(csv_df))
            csv_columns = [str(column) for column in csv_df.columns]
            csv_columns_hash = _sha256_text(json.dumps(csv_columns))
        else:
            csv_rows = 0
            csv_columns = []
            csv_columns_hash = ""

        sqlite_rows = int(reg["row_count"]) if reg is not None else 0
        sqlite_columns = []
        if reg is not None:
            try:
                sqlite_columns = json.loads(str(reg.get("columns_json") or "[]"))
            except Exception:
                sqlite_columns = []
        sqlite_columns_hash = _sha256_text(json.dumps(sqlite_columns)) if sqlite_columns else ""

        column_match = bool(csv_columns_hash and sqlite_columns_hash and csv_columns_hash == sqlite_columns_hash)
        if not csv_exists:
            status = "missing_csv"
            note = "CSV file not found"
        elif reg is None:
            status = "missing_sqlite"
            note = "Dataset missing in SQLite registry"
        elif csv_rows != sqlite_rows:
            status = "row_mismatch"
            note = f"csv_rows={csv_rows} sqlite_rows={sqlite_rows}"
        elif not column_match:
            status = "column_mismatch"
            note = "CSV and SQLite column signatures differ"
        else:
            status = "ok"
            note = "row and column parity match"

        rows.append(
            ParityResult(
                dataset_key=str(dataset_key),
                csv_path=str(csv_path),
                csv_exists=csv_exists,
                csv_rows=csv_rows,
                sqlite_rows=sqlite_rows,
                column_match=column_match,
                status=status,
                note=note,
            ).to_dict()
        )

    return pd.DataFrame(rows)
