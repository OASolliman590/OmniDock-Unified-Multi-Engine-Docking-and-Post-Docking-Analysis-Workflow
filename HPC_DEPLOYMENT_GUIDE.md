# HPC Deployment Guide

This pipeline supports secure HPC deployment by separating public templates from local site-specific profiles.

## Security model

- Keep public templates in `examples/hpc_profiles/`.
- Keep real cluster settings in `.workflow/hpc_profiles/`.
- Do not commit Slurm accounts, home-directory paths, or scheduler conventions into tracked files.
- Generated docking projects and deployment bundles already live under ignored runtime roots.

## Profiles

- Public-safe template: `examples/hpc_profiles/bibalex-apptainer-conda.template.json`
- Local BibAlex profile: `.workflow/hpc_profiles/alex-bibalex.json`

Generate deployment assets with the local profile:

```bash
python main.py dock deploy \
  --project-dir /path/to/project \
  --engines gnina,vina,smina \
  --hpc-profile alex-bibalex
```

Use the public-safe builtin template instead:

```bash
python main.py dock deploy \
  --project-dir /path/to/project \
  --engines gnina,vina,smina \
  --hpc-profile bibalex-apptainer-conda
```

## Controller environment

Use the project Conda environment for project orchestration, preparation, and post-docking analysis:

```bash
conda env create -f environment.yml
conda activate pdb-prepare-wizard
pip install -e .
```

## Engine bootstrap on HPC

### GNINA

The pipeline expects a GNINA binary and, for the current BibAlex setup, an Apptainer image:

```bash
ls -l "$HOME/gnina"
ls -l "$HOME/cuda12.3.2-cudnn9-runtime-ubuntu22.04.sif"
apptainer exec --nv "$HOME/cuda12.3.2-cudnn9-runtime-ubuntu22.04.sif" "$HOME/gnina" --version
```

### AutoDock Vina CLI

The Python `vina` package is not enough for this pipeline. Install the CLI executable separately:

```bash
mkdir -p "$HOME/.local/bin"
curl -L \
  https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64 \
  -o "$HOME/.local/bin/vina"
chmod +x "$HOME/.local/bin/vina"
"$HOME/.local/bin/vina" --help | head
```

If you want to keep the existing `vina_dock` environment for helper Python packages, retain it, but the docking pipeline still needs the CLI binary above available on `PATH` or referenced explicitly.

For the local BibAlex deployment profile, the simplest contract is to point the
profile directly at `$HOME/.local/bin/vina` rather than relying on `conda run`.

### Smina CLI

Install the static binary and expose it on `PATH`:

```bash
mkdir -p "$HOME/.local/bin"
curl -L https://sourceforge.net/projects/smina/files/smina.static/download -o "$HOME/.local/bin/smina"
chmod +x "$HOME/.local/bin/smina"
"$HOME/.local/bin/smina" --help | head
```

### PATH wiring

Add user binaries to the shell startup once:

```bash
grep -q 'export PATH="$HOME/.local/bin:$PATH"' ~/.bashrc || \
  echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

For the local BibAlex deployment profile, point the Smina runtime to
`$HOME/.local/bin/smina`.

## Verification

Check the final runtime surface on the cluster:

```bash
command -v apptainer
command -v conda
command -v vina
command -v smina
ls -l "$HOME/gnina"
ls -l "$HOME/cuda12.3.2-cudnn9-runtime-ubuntu22.04.sif"
```

Then generate Slurm bundles:

```bash
python main.py dock deploy \
  --project-dir /path/to/project \
  --engines gnina,vina,smina \
  --hpc-profile alex-bibalex
```

Generated scripts will include:

- shell bootstrap via `~/.bashrc`
- Apptainer module-loading fallback for GNINA
- engine-specific Slurm defaults from the selected profile
- one Slurm array job per engine, not one `sbatch` call per docking pair
- optional array parallelism caps from the selected profile or `--slurm-array-parallelism`

## Screening vs exhaustive

- `screen` is the normal first pass for the full `pairlist.csv`.
- `exhaustive` is the follow-up stage for a smaller rerun manifest produced by
  comparative analysis.
- The CLI now refuses `dock deploy --mode exhaustive` unless you pass
  `--from-rerun-manifest` or an explicit `--allow-full-exhaustive` override.
- The CLI also refuses `dock submit --mode exhaustive` if the local deployment
  bundle does not record a rerun manifest, unless you pass the same override.
- The interactive workflow mirrors those checks and will warn before allowing a
  deliberate full-pair exhaustive campaign.

Generate a proper exhaustive rerun manifest after screening:

```bash
python main.py analyze comparative \
  --project-dir /path/to/project \
  --promote-exhaustive \
  --rerun-engine gnina
```

Then deploy from that manifest:

```bash
python main.py dock deploy \
  --project-dir /path/to/project \
  --mode exhaustive \
  --from-rerun-manifest /path/to/exhaustive_rerun_manifest.csv \
  --hpc-profile alex-bibalex
```

## Upstream engine contract verification

The current wrappers were checked against the upstream GitHub documentation on
March 9, 2026:

- `gnina/gnina` README documents both `--out` and `--log`, so the GNINA runner
  still writes a dedicated log file.
- `ccsb-scripps/AutoDock-Vina` documents `--out` for single-ligand docking and
  `--dir` for batch mode. No documented `--log` flag was found in the official
  repo docs, so the Vina runner no longer passes `--log`.
- The available GitHub mirror README for `smina` documents both `--out` and
  `--log`, so the Smina runner remains compatible with that contract. The
  original project homepage is still SourceForge.

## Current HPC observations

On `<hpc-login-host>` as of March 9, 2026:

- `apptainer` and `singularity` are available.
- `conda` is available.
- `$HOME/gnina` exists.
- `$HOME/cuda12.3.2-cudnn9-runtime-ubuntu22.04.sif` exists.
- Conda env `vina_dock` exists, but it currently contains the Python `vina` package rather than a `vina` CLI binary.
- No `smina` executable was found in the home directory or shared software paths that were inspected.
