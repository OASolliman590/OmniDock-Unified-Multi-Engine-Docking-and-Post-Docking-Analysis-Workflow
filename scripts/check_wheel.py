"""Check installable entrypoints and required runtime resources in a built wheel."""
from pathlib import Path
from zipfile import ZipFile
import configparser
import json
import subprocess
import sys
import tempfile

wheels = sorted(Path("dist").glob("*.whl"), key=lambda path: path.stat().st_mtime)
if not wheels:
    raise SystemExit("Build a wheel before running this check")
with ZipFile(wheels[-1]) as archive:
    names = set(archive.namelist())
    required = {
        "main.py", "cli_pipeline.py", "interactive_pipeline.py", "core_pipeline.py",
        "batch_pdb_preparation.py", "autodock_preparation.py", "prep_autodock_enhanced.sh",
        "workflow/cli.py", "workflow/state.py", "post_docking_analysis/config/schema.yaml",
    }
    missing = required - names
    if missing:
        raise SystemExit(f"Wheel omits runtime files: {sorted(missing)}")
    entries = next(name for name in names if name.endswith(".dist-info/entry_points.txt"))
    config = configparser.ConfigParser()
    config.read_string(archive.read(entries).decode())
    for line in archive.read(entries).decode().splitlines():
        if "=" not in line:
            continue
        module = line.split("=", 1)[1].split(":", 1)[0].strip().replace(".", "/")
        if module + ".py" not in names and module + "/__init__.py" not in names:
            raise SystemExit(f"Console entrypoint module absent: {module}")
    with tempfile.TemporaryDirectory(prefix="omnidock-wheel-") as directory:
        archive.extractall(directory)
        code = '''
import importlib, json, pathlib, sys
root = pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root))
for entry in json.loads(sys.argv[2]):
    name, attribute = entry.split(":", 1)
    module = importlib.import_module(name)
    assert root in pathlib.Path(module.__file__).resolve().parents, (name, module.__file__)
    assert callable(getattr(module, attribute)), entry
print("Imported every console entrypoint from the extracted wheel outside the checkout")
'''
        subprocess.run([sys.executable, "-I", "-c", code, directory,
                        json.dumps(list(config["console_scripts"].values()))],
                       cwd=directory, check=True)
print(f"Verified runtime files and console entrypoints: {wheels[-1].name}")
