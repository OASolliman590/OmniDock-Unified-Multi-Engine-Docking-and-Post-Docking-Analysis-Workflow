#!/usr/bin/env python3
"""
Omni-DockForge - Unified entry point.
"""

from __future__ import annotations

import sys

from workflow.cli import main as workflow_main


LEGACY_TOP_LEVEL = {"interactive", "cli", "batch", "prepare-docking", "dock"}
DOCK_SUBCOMMANDS = {"run", "dry-run", "engine", "deploy", "sync", "submit"}


def _ensure_console_encoding() -> None:
    """Make stdout/stderr able to carry the non-ASCII output this CLI emits.

    A Windows console defaults to cp1252, which cannot encode the banner or
    the status glyphs used throughout. Without this, every command --
    including `--help` -- dies with UnicodeEncodeError before doing any
    work. Reconfiguring beats requiring users to set PYTHONIOENCODING.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            encoding = (getattr(stream, "encoding", "") or "").lower()
            if encoding.replace("-", "") != "utf8" and hasattr(stream, "reconfigure"):
                stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            # Never let console setup stop the program from running.
            pass


def _has_tty() -> bool:
    return sys.stdin.isatty() and sys.stdout.isatty()


def _normalize_legacy_args(argv: list[str]) -> list[str]:
    if not argv:
        return ["workflow", "interactive"] if _has_tty() else ["--help"]

    first = argv[0]
    if first in {"-h", "--help"}:
        return argv

    if first == "interactive":
        return ["workflow", "interactive", *argv[1:]]

    if first == "cli":
        return ["pdb", "run", *argv[1:]]

    if first == "batch":
        return ["pdb", "batch", *argv[1:]]

    if first == "prepare-docking":
        return ["prep", "project", *argv[1:]]

    if first == "dock":
        if len(argv) == 1:
            return ["dock", "run", *argv[1:]]
        if argv[1] in {"-h", "--help"}:
            return ["dock", *argv[1:]]
        if argv[1].startswith("-"):
            return ["dock", "run", *argv[1:]]
        if argv[1] not in DOCK_SUBCOMMANDS:
            return ["dock", "run", *argv[1:]]
        return argv

    return argv


def main() -> int:
    _ensure_console_encoding()
    argv = _normalize_legacy_args(sys.argv[1:])
    print("🔬 Omni-DockForge: End-to-End Docking, Consensus by Design. v4.0.0-dev")
    print("=" * 40)
    return workflow_main(argv)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\n\n⚠️  Pipeline interrupted by user")
        raise SystemExit(1)
    except SystemExit:
        raise
    except Exception as exc:
        print(f"\n❌ Unexpected error: {exc}")
        raise SystemExit(1)
