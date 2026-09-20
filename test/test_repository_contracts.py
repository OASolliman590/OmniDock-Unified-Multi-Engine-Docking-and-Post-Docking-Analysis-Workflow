"""Repository integrity contracts for tracked links, artifacts, and symlinks."""
from dataclasses import FrozenInstanceError
from pathlib import Path
import os
import subprocess
import sys

import pytest

from scripts.check_repository import Finding, check_repository


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_repository.py"


def run_git(root, *args):
    result = subprocess.run(["git", *args], cwd=root, capture_output=True, text=True)
    if result.returncode != 0:
        raise AssertionError(f"git {args} failed: {result.stderr or result.stdout}")
    return result


def init_repo(tmp_path):
    root = tmp_path / "repo"
    root.mkdir(parents=True)
    run_git(root, "init")
    run_git(root, "config", "user.email", "contract-test@example.com")
    run_git(root, "config", "user.name", "Contract Test")
    run_git(root, "config", "core.autocrlf", "false")
    return root


def track(root, relpath, content=""):
    path = root / relpath
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(content, bytes):
        path.write_bytes(content)
    else:
        path.write_text(content, encoding="utf-8")
    run_git(root, "add", "--", Path(relpath).as_posix())
    return path


def track_git_symlink(root, relpath, target):
    path = root / relpath
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(target.encode("utf-8"))
    posix = Path(relpath).as_posix()
    blob = run_git(root, "hash-object", "-w", "--", posix).stdout.strip()
    run_git(root, "update-index", "--add", "--cacheinfo", f"120000,{blob},{posix}")
    return path


def run_cli(*args, cwd=None):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


def finding_codes(findings):
    return {item.code for item in findings}


def test_finding_is_an_immutable_dataclass():
    finding = Finding(code="broken-markdown-link", path="README.md", detail="missing")
    assert (finding.code, finding.path, finding.detail) == (
        "broken-markdown-link",
        "README.md",
        "missing",
    )
    with pytest.raises(FrozenInstanceError):
        finding.code = "other"


def test_repository_checker_reports_broken_links_generated_logs_and_dangling_links(tmp_path):
    root = init_repo(tmp_path)
    track(root, "README.md", "See the [missing page](missing.md).\n")
    track(root, "audit/build-correctness.log", "generated log\n")
    track_git_symlink(root, "dangling", "/this/symlink/target/does/not/exist")
    findings = check_repository(root)
    assert finding_codes(findings) >= {"broken-markdown-link", "tracked-generated-artifact", "dangling-symlink"}
    assert any(item.path == "README.md" and item.code == "broken-markdown-link" for item in findings)
    assert any(item.path == "audit/build-correctness.log" for item in findings)
    assert any(item.path == "dangling" and item.code == "dangling-symlink" for item in findings)


def test_markdown_link_resolution_and_ignores(tmp_path):
    root = init_repo(tmp_path)
    track(root, "present.md", "# present\n")
    track(root, "subdir/kept.txt", "kept\n")
    track(root, "folder name/a.md", "# spaced\n")
    track(root, "docs/inner.md", "[up](../present.md) [bad](../nope.md)\n")
    track(
        root,
        "README.md",
        "\n".join(
            [
                "See [ok](present.md), [dir](subdir/), [spaces](folder%20name/a.md),",
                "[quoted](<folder name/a.md>), [frag](present.md#ignored-anchor),",
                "[query](present.md?raw=1).",
                "External: [http](http://example.com/missing.md), [https](https://example.com/x),",
                "[mailto](mailto:dev@example.com), [data](data:text/plain,hello),",
                "[proto](//cdn.example.com/x.png), [anchor](#local-only).",
                "Image: ![missing](no-such.png)",
                "```",
                "[fenced](also-missing.md)",
                "```",
                "Broken: [gone](absent.md)",
                "",
                "[ref-ok]: present.md",
                "[ref-bad]: missing-ref.md",
                "",
                "Use [okref][ref-ok] and [badref][ref-bad].",
                "",
            ]
        )
        + "\n",
    )
    track(root, "untracked-hint.md", "ignore me")
    (root / "untracked.md").write_text("[broken](nope.md)\n", encoding="utf-8")
    findings = [item for item in check_repository(root) if item.code == "broken-markdown-link"]
    missing = {(item.path, item.detail) for item in findings}
    assert any(path == "README.md" and "absent.md" in detail for path, detail in missing)
    assert any(path == "README.md" and "missing-ref.md" in detail for path, detail in missing)
    assert any(path == "docs/inner.md" and "nope.md" in detail for path, detail in missing)
    details = " ".join(detail for _, detail in missing)
    for ignored in (
        "no-such.png",
        "also-missing.md",
        "http://example.com/missing.md",
        "https://example.com/x",
        "mailto:dev@example.com",
        "data:text/plain,hello",
        "//cdn.example.com/x.png",
        "#local-only",
        "present.md",
        "subdir/",
        "folder%20name/a.md",
        "folder name/a.md",
        "untracked.md",
        "nope.md",
    ):
        if ignored == "nope.md":
            continue
        assert ignored not in details
    assert not any(item.path == "untracked.md" for item in findings)
    assert not any("nope.md" in item.detail and item.path == "untracked.md" for item in findings)


def test_image_only_reference_definitions_are_ignored(tmp_path):
    root = init_repo(tmp_path)
    track(
        root,
        "README.md",
        "\n".join(
            [
                "![plot][Figure]",
                "![Collapsed][]",
                "![shared][Both Uses]",
                "See the [shared image][both  uses] and the [guide][Missing Guide].",
                "",
                "[figure]: figure-only.png",
                "[collapsed]: collapsed-only.png",
                "[both uses]: shared-target.png",
                "[missing guide]: gone.md",
                "",
            ]
        )
        + "\n",
    )
    findings = [item for item in check_repository(root) if item.code == "broken-markdown-link"]
    details = " ".join(item.detail for item in findings)
    assert "gone.md" in details
    assert "shared-target.png" in details
    assert "figure-only.png" not in details
    assert "collapsed-only.png" not in details


def test_generated_artifact_allowlist(tmp_path):
    root = init_repo(tmp_path)
    allowlisted = [
        "audit/build-correctness.log",
        "audit/correctness-tests.log",
        "audit/correctness-tests.xml",
        "audit/docking-dry-run.log",
        "audit/wheel-verification.log",
    ]
    for relpath in allowlisted:
        track(root, relpath, "generated\n")
    track(root, "audit/docking-probe-results.json", "{}\n")
    track(root, "audit/notes.log", "research note\n")
    track(root, "tmp/wheel-verification.log", "wrong place\n")
    (root / "audit" / "untracked-correctness-tests.log").write_text("not tracked\n", encoding="utf-8")
    findings = [item for item in check_repository(root) if item.code == "tracked-generated-artifact"]
    assert sorted(item.path for item in findings) == allowlisted
    assert finding_codes(findings) == {"tracked-generated-artifact"}


def test_findings_are_sorted_deterministically(tmp_path):
    root = init_repo(tmp_path)
    track(root, "z.md", "[a](missing-z.md)\n")
    track(root, "a.md", "[b](missing-a.md) [c](missing-c.md)\n")
    track(root, "audit/wheel-verification.log", "log\n")
    track(root, "audit/build-correctness.log", "log\n")
    findings = check_repository(root)
    assert findings == sorted(findings, key=lambda item: (item.path, item.code, item.detail))
    assert [item.path for item in findings] == sorted(item.path for item in findings)


def test_cli_exits_zero_only_when_clean(tmp_path):
    clean = init_repo(tmp_path / "clean")
    track(clean, "README.md", "See [ok](present.md).\n")
    track(clean, "present.md", "ok\n")
    clean_result = run_cli(str(clean))
    assert clean_result.returncode == 0, clean_result.stderr
    assert clean_result.stdout == ""

    dirty = init_repo(tmp_path / "dirty")
    track(dirty, "README.md", "[missing](gone.md)\n")
    track(dirty, "audit/docking-dry-run.log", "log\n")
    dirty_result = run_cli(str(dirty))
    assert dirty_result.returncode != 0
    lines = [line for line in dirty_result.stdout.splitlines() if line.strip()]
    assert lines == sorted(lines)
    assert any("broken-markdown-link" in line and "README.md" in line for line in lines)
    assert any("tracked-generated-artifact" in line and "audit/docking-dry-run.log" in line for line in lines)


def test_cli_defaults_to_repository_containing_the_script():
    repo = Path(__file__).resolve().parents[1]
    explicit = run_cli(str(repo))
    default = run_cli()
    assert explicit.returncode == default.returncode
    assert explicit.stdout == default.stdout
    assert explicit.stderr == default.stderr


def test_non_git_root_fails_clearly(tmp_path):
    missing = tmp_path / "not-a-repo"
    missing.mkdir()
    with pytest.raises(RuntimeError, match="[Gg]it"):
        check_repository(missing)
    result = run_cli(str(missing))
    assert result.returncode != 0
    assert "git" in result.stderr.lower()


def test_dangling_symlink_when_supported(tmp_path):
    root = init_repo(tmp_path)
    track(root, "target.txt", "here\n")
    track_git_symlink(root, "valid-link", "target.txt")
    track_git_symlink(root, "missing-link", "no-such-target")
    findings = [item for item in check_repository(root) if item.code == "dangling-symlink"]
    assert {item.path for item in findings} == {"missing-link"}

    os_root = init_repo(tmp_path / "os-symlinks")
    track(os_root, "kept.txt", "kept\n")
    link = os_root / "os-dangling"
    try:
        os.symlink("definitely-missing-os-target", link)
    except OSError as exc:
        pytest.skip(f"creating symlinks is unavailable: {exc}")
    if not link.is_symlink():
        pytest.skip("creating symlinks is unavailable")
    run_git(os_root, "add", "--", "os-dangling")
    os_findings = check_repository(os_root)
    assert any(item.code == "dangling-symlink" and item.path == "os-dangling" for item in os_findings)


def test_cli_parser_module_accepts_representative_public_commands():
    from workflow.cli_parser import build_parser

    parser = build_parser()
    assert parser.parse_args(["workflow", "status", "--project-dir", "p"]).workflow_command == "status"
    assert parser.parse_args(["dock", "run", "--project-dir", "p"]).dock_command == "run"
    assert parser.parse_args(["analyze", "comparative", "--project-dir", "p"]).analyze_command == "comparative"
