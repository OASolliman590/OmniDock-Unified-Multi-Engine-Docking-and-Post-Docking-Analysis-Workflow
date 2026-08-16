"""Flask application factory and route registration.

Routes only -- all behavior lives in the adapter modules. This is the one
module in the package permitted to import Flask.
"""

from __future__ import annotations

from flask import Flask, jsonify

from . import __version__
from .config import ensure_app_dirs


def create_app() -> Flask:
    app = Flask(__name__)
    ensure_app_dirs()

    @app.get("/api/health")
    def api_health():
        return jsonify({"ok": True, "version": __version__})

    return app
