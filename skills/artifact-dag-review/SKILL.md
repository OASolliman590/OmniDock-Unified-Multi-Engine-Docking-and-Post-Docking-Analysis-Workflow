---
name: artifact-dag-review
description: Audit the artifact dependency graph, node execution states, cache invalidation rules, and execution plans.
---

# Artifact DAG Review Skill

## Trigger Conditions
Activate when updating DAG execution nodes, cache key generation, or selective stage reruns.

## Explicit Inputs and Outputs
- **Inputs:** `post_docking_analysis/artifact_graph.py`, project directory, requested targets.
- **Outputs:** Execution plan, cache hit/miss report, DAG status dictionary.

## Refusal-to-Guess Rules
- **NEVER** allow a cache hit if input files, parameters, or code revisions have changed.
- **NEVER** execute downstream nodes if their mandatory upstream prerequisites failed.

## Verification Commands
- `pytest test/test_webui.py -k dag`
