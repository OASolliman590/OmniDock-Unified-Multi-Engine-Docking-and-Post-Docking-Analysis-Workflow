"""OmniDock web UI: a browser surface over the existing workflow and CLI.

This package is a thin adapter. It reads project state through
``workflow.state`` accessors, sources its vocabulary from
``workflow.interactive`` constants, and executes work by launching the
existing CLI as a subprocess. No pipeline stage is reimplemented here.

See ``specs/030-web-ui-platform/`` for the governing specification.
"""

__version__ = "0.1.0"
