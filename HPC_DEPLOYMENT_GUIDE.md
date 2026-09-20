# HPC deployment

OmniDock generates a portable deployment bundle locally, then optionally syncs
and submits it to a remote scheduler. Generation is not execution, and a
successful submission is not evidence that engine outputs passed validation.

## Profiles

Start from a public-safe template:

- `examples/hpc_profiles/bibalex-apptainer-conda.template.json`: generic Slurm
  with engine/container/Conda placeholders;
- `examples/hpc_profiles/nmrbox-condor.template.json`: NMRBox-style HTCondor
  placeholders;
- `examples/hpc_profiles/ssh_systems.template.yaml`: example SSH system routing.

Copy a template to a private local file, fill in values, and pass
`--hpc-profile-file`. Do not commit usernames, hostnames, account IDs, absolute
personal paths, tokens, keys, or filled connection profiles. The templates are
examples, not verified allocations or live connection details.

## Slurm flow

```bash
python main.py dock deploy \
  --project-dir docking_project --engines gnina,vina \
  --mode screen --hpc-profile-file local-slurm-profile.json

python main.py dock sync \
  --project-dir docking_project \
  --hpc-profile-file local-slurm-profile.json --dry-run

python main.py dock submit \
  --project-dir docking_project --mode screen \
  --hpc-profile-file local-slurm-profile.json --dry-run
```

Deployment resolves engine runtimes and Slurm resources, writes commands/jobs
and scheduler manifests below the project deployment area, and keeps paths
portable. Review the pair count, engines, boxes, parameters, runtime identities,
resource requests, remote root, sync exclusions, and generated manifest before
removing either dry-run flag.

## NMRBox / HTCondor flow

Use the public `nmrbox-condor.template.json` only as a schema. Set private SSH,
remote-root, binary/image, CPU/GPU, memory, disk, and parallelism values in an
untracked copy. Then use the same `deploy`, `sync`, and `submit` verbs with that
file. OmniDock selects the scheduler from the profile. This guide does not claim
that NMRBox credentials, queues, binaries, or GPU access are available.

## Screen and exhaustive safeguards

`--mode screen` is the normal first deployment. Exhaustive deployment should
consume a reviewed rerun manifest produced by comparative analysis:

```bash
python main.py analyze comparative \
  --project-dir docking_project \
  --promote-exhaustive --rerun-engine gnina --top-per-protein 3

python main.py dock deploy \
  --project-dir docking_project --mode exhaustive \
  --from-rerun-manifest path/to/rerun_manifest.csv \
  --hpc-profile-file local-profile.json
```

The explicit `--allow-full-exhaustive` option bypasses this safeguard. Use it
only after deliberate review of scope and cost. A rerun manifest is a selection
record, not a scientific endorsement.

## Sync, submit, and recovery

`dock sync` transfers the prepared project/deployment using the profile or
explicit SSH arguments. `--delete` can remove remote files absent locally; do
not use it without reviewing the dry run and destination. `dock submit` submits
a previously generated/synced deployment and supports `--dry-run`.

After execution, preserve scheduler logs and sync results back before analysis.
Resume requires matching inputs, executable identity, protocol, and validated
outputs. A pose filename alone is insufficient. Inspect failed/empty attempts
and locks; never relabel old outputs as the current attempt.

## Security and environment review

- use SSH keys/agents and host verification; never place secrets in profiles;
- restrict profile and scheduler-log permissions;
- pin/record containers, binaries, modules, Conda environments, and parameter
  files;
- validate Python 3.10+ and every engine on worker nodes;
- test a small screen bundle before scaling;
- report excluded/failed jobs and coverage.

See [dependencies](DEPENDENCIES.md), [CLI recipes](USAGE.md), and
[testing/evidence](docs/testing.md).
