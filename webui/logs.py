"""Incremental reads over job log files.

Docking runs produce large logs over many hours, so the log view tails by
byte offset rather than loading the file. Nothing here ever reads the whole
file into memory (FR-015).
"""

from __future__ import annotations

from pathlib import Path

# Per-request read ceiling. Large enough to catch up quickly after a
# reconnect, small enough that one request cannot exhaust memory.
MAX_CHUNK_BYTES = 256 * 1024

# Initial window shown when a page loads.
DEFAULT_TAIL_BYTES = 64 * 1024


def size_of(log_path: str | Path) -> int:
    try:
        return Path(log_path).stat().st_size
    except OSError:
        return 0


def read_from(
    log_path: str | Path,
    offset: int = 0,
    max_bytes: int = MAX_CHUNK_BYTES,
) -> tuple[str, int, int]:
    """Read up to ``max_bytes`` starting at ``offset``.

    Returns ``(text, new_offset, total_size)``. Decoding is lenient because
    a chunk boundary can split a multi-byte character.
    """
    path = Path(log_path)
    total = size_of(path)

    if total == 0 or offset >= total:
        return "", max(0, min(offset, total)), total

    # A truncated or replaced file (offset past EOF) restarts from 0.
    if offset < 0:
        offset = 0

    try:
        with open(path, "rb") as handle:
            handle.seek(offset)
            data = handle.read(max_bytes)
    except OSError:
        return "", offset, total

    return data.decode("utf-8", errors="replace"), offset + len(data), total


def tail(log_path: str | Path, max_bytes: int = DEFAULT_TAIL_BYTES) -> tuple[str, int, int]:
    """Read the last ``max_bytes`` of the log, for the initial page load."""
    total = size_of(log_path)
    start = max(0, total - max_bytes)
    text, new_offset, total = read_from(log_path, start, max_bytes)
    if start > 0:
        text = f"… showing the last {len(text)} bytes of {total} …\n{text}"
    return text, new_offset, total
