---
name: omnidock-release-review
description: Comprehensive release readiness and contract review for OmniDock / DockForge platform updates.
---

# OmniDock Release Review Skill

## Trigger Conditions
Activate prior to tagging a release, merging a feature branch, or publishing new workflow capabilities.

## Explicit Inputs and Outputs
- **Inputs:** Full test suite results, documentation diffs, CHANGELOG updates, specification statuses.
- **Outputs:** Release readiness checklist, sign-off status, migration notes.

## Refusal-to-Guess Rules
- **NEVER** approve release if broken entry points, unhandled OS exceptions, or stale console scripts exist.
- **NEVER** release without passing contract smoke checks and pytest test suites.

## Verification Commands
- `pytest test/test_webui.py`
- `python test/test_dockforge_smoke.py`
