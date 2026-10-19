"""Reentrant file locks shared by threads and independent local processes."""
from __future__ import annotations

import os
import threading
import time
from contextlib import contextmanager
from pathlib import Path

_locks_guard = threading.Lock()
_locks = {}
_held = threading.local()


@contextmanager
def _thread_lock(lock, deadline, path):
    if not lock.acquire(timeout=max(0.0, deadline - time.monotonic())):
        raise TimeoutError(f"Timed out waiting for workflow lock: {path}")
    try:
        yield
    finally:
        lock.release()


@contextmanager
def file_lock(path: Path, timeout: float = 60.0):
    """Hold an OS lock on a stable sidecar, never on a replaced data file.

    Locks are per path, so a graph worker may update unrelated workflow state.
    The same thread may nest transactions without releasing the outer lock.
    """
    path = Path(path).resolve()
    key = os.path.normcase(str(path))
    with _locks_guard:
        local_lock = _locks.setdefault(key, threading.RLock())
    deadline = time.monotonic() + timeout
    with _thread_lock(local_lock, deadline, path):
        held = getattr(_held, "paths", None)
        if held is None:
            held = _held.paths = set()
        if key in held:
            yield
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a+b") as handle:
            if handle.seek(0, os.SEEK_END) == 0:
                handle.write(b"\0")
                handle.flush()
            while True:
                try:
                    handle.seek(0)
                    if os.name == "nt":
                        import msvcrt
                        msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                    else:
                        import fcntl
                        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except (OSError, BlockingIOError):
                    if time.monotonic() >= deadline:
                        raise TimeoutError(f"Timed out waiting for workflow lock: {path}")
                    time.sleep(0.025)
            held.add(key)
            try:
                yield
            finally:
                held.remove(key)
                handle.seek(0)
                if os.name == "nt":
                    import msvcrt
                    msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
                else:
                    import fcntl
                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
