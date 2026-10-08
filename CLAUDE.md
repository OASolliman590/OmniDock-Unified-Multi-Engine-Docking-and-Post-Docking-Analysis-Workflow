# OmniDock Developer & Agent Instructions

See [AGENTS.md](./AGENTS.md) for the authoritative **Agentic Bioinformatics Operating Contract**, **Karpathy Coding Guidelines**, scientific interpretation safeguards, and output topology contracts.

## Quick Test Commands
- Pytest WebUI tests: pytest test/test_webui.py
- Contract Smoke tests: python test/test_dockforge_smoke.py --skip-all-engines --skip-prep-matrix
- Full dry-run smoke test: python test/test_dockforge_smoke.py


## Omar progress reporting

For Omar-facing reports and Notion updates, read docs/OMAR_PROGRESS_CONTRACT.md. Keep a plain-language briefing and separate software evidence from scientific evidence. First adoption is orientation only. Preserve the existing scientific spec, task IDs, approval gates and orchestration limits. Only the coordinator publishes shared Notion progress. Never claim a sync without read-back verification.
