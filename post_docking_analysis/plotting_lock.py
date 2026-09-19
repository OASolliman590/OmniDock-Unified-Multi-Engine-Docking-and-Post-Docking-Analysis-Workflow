"""Process-wide serialization for Matplotlib's pyplot state."""

from __future__ import annotations

from contextlib import contextmanager
from functools import wraps
from threading import RLock
from typing import Callable, Iterator, ParamSpec, TypeVar


_P = ParamSpec("_P")
_R = TypeVar("_R")

_MATPLOTLIB_LOCK = RLock()


@contextmanager
def matplotlib_critical_section() -> Iterator[None]:
    """Hold the shared reentrant lock while touching pyplot global state."""

    with _MATPLOTLIB_LOCK:
        yield


def serialized_matplotlib(function: Callable[_P, _R]) -> Callable[_P, _R]:
    """Run one complete pyplot operation under the process-wide lock."""

    @wraps(function)
    def wrapped(*args: _P.args, **kwargs: _P.kwargs) -> _R:
        with matplotlib_critical_section():
            return function(*args, **kwargs)

    return wrapped


__all__ = ["matplotlib_critical_section", "serialized_matplotlib"]
