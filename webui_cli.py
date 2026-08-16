#!/usr/bin/env python3
"""Entry point for the OmniDock web UI.

    python webui_cli.py
    python main.py webui

Binds to localhost by default. There is no authentication in this feature,
so exposing it on another interface requires the explicit --allow-remote
flag (FR-034).
"""

from __future__ import annotations

import argparse
import sys

from webui.app import create_app
from webui.config import HOST, PORT


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="omnidock-webui",
        description="Run the OmniDock web UI (localhost only by default).",
    )
    parser.add_argument("--host", default=HOST, help=f"Bind address (default: {HOST})")
    parser.add_argument("--port", type=int, default=PORT, help=f"Bind port (default: {PORT})")
    parser.add_argument(
        "--allow-remote",
        action="store_true",
        help="Permit binding to a non-localhost address. The UI has no authentication; "
        "only use this on a trusted network.",
    )
    parser.add_argument("--debug", action="store_true", help="Run Flask in debug mode")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    is_local = args.host in ("127.0.0.1", "localhost", "::1")
    if not is_local and not args.allow_remote:
        print(
            f"error: refusing to bind to {args.host} without --allow-remote.\n"
            "The web UI has no authentication; anyone who can reach this address "
            "could launch jobs and read project files.",
            file=sys.stderr,
        )
        return 2

    app = create_app()
    print(f"OmniDock web UI: http://{args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=args.debug, threaded=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
