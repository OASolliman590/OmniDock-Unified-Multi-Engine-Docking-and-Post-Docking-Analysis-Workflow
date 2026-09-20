"""Check tracked Markdown links, generated artifacts, and dangling symlinks."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse
import os
import re
import subprocess
import sys
import urllib.parse


GENERATED_ARTIFACTS = frozenset(
    {
        "audit/build-correctness.log",
        "audit/correctness-tests.log",
        "audit/correctness-tests.xml",
        "audit/docking-dry-run.log",
        "audit/wheel-verification.log",
    }
)
MARKDOWN_SUFFIXES = {".md", ".markdown"}
_OPEN_FENCE = re.compile(r"^( {0,3})(`{3,}|~{3,})")
_REF_DEF = re.compile(r"^ {0,3}\[([^\]]+)\]:\s*(<[^>\n]*>|\S+)", re.MULTILINE)


@dataclass(frozen=True)
class Finding:
    code: str
    path: str
    detail: str


def check_repository(root: Path) -> list[Finding]:
    root = Path(root).expanduser().resolve()
    entries = _tracked_entries(root)
    findings = [
        *_broken_markdown_links(root, entries),
        *_tracked_generated_artifacts(entries),
        *_dangling_symlinks(root, entries),
    ]
    return sorted(set(findings), key=lambda item: (item.path, item.code, item.detail))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Check repository contracts for Markdown links, generated artifacts, and dangling symlinks.",
    )
    parser.add_argument(
        "root",
        nargs="?",
        default=None,
        help="Repository root (default: repository containing this script)",
    )
    args = parser.parse_args(argv)
    root = Path(args.root) if args.root is not None else Path(__file__).resolve().parent.parent
    try:
        findings = check_repository(root)
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    for finding in findings:
        print(f"{finding.path}: {finding.code}: {finding.detail}")
    return 1 if findings else 0


def _tracked_entries(root: Path) -> list[tuple[str, str]]:
    probe = _run_git(root, "rev-parse", "--is-inside-work-tree", text=True)
    if probe.returncode != 0 or probe.stdout.strip() != "true":
        message = (probe.stderr or probe.stdout or "").strip() or f"not a Git work tree: {root}"
        raise RuntimeError(message)
    listed = _run_git(root, "ls-files", "-s", "-z", "--full-name")
    if listed.returncode != 0:
        message = (listed.stderr or listed.stdout).decode("utf-8", "replace").strip()
        raise RuntimeError(message or f"git ls-files failed in {root}")
    entries: list[tuple[str, str]] = []
    for raw in listed.stdout.split(b"\0"):
        if not raw:
            continue
        try:
            meta, path_bytes = raw.split(b"\t", 1)
        except ValueError as exc:
            raise RuntimeError("unexpected git ls-files -s output") from exc
        mode = meta.split(b" ", 1)[0].decode("ascii")
        entries.append((mode, os.fsdecode(path_bytes).replace("\\", "/")))
    return entries


def _run_git(root: Path, *args: str, text: bool = False) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(
            ["git", "-C", str(root), *args],
            capture_output=True,
            text=text,
        )
    except FileNotFoundError as exc:
        raise RuntimeError("git executable not found") from exc


def _broken_markdown_links(root: Path, entries: list[tuple[str, str]]) -> list[Finding]:
    findings: list[Finding] = []
    for _, path in entries:
        if Path(path).suffix.lower() not in MARKDOWN_SUFFIXES:
            continue
        source = root / path
        try:
            text = source.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            text = source.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        seen: set[str] = set()
        for dest in _markdown_destinations(text):
            if dest in seen:
                continue
            if not _missing_local_target(source, dest):
                continue
            seen.add(dest)
            findings.append(
                Finding(
                    code="broken-markdown-link",
                    path=path,
                    detail=f"missing local target {dest!r}",
                )
            )
    return findings


def _tracked_generated_artifacts(entries: list[tuple[str, str]]) -> list[Finding]:
    findings: list[Finding] = []
    for _, path in entries:
        if path in GENERATED_ARTIFACTS:
            findings.append(
                Finding(
                    code="tracked-generated-artifact",
                    path=path,
                    detail="generated verification artifact must not be tracked",
                )
            )
    return findings


def _dangling_symlinks(root: Path, entries: list[tuple[str, str]]) -> list[Finding]:
    findings: list[Finding] = []
    for mode, path in entries:
        full = root / path
        try:
            fs_symlink = full.is_symlink()
        except OSError:
            fs_symlink = False
        if mode != "120000" and not fs_symlink:
            continue
        target_text = _symlink_target_text(full)
        if target_text is None:
            findings.append(
                Finding(
                    code="dangling-symlink",
                    path=path,
                    detail="symbolic link target does not exist",
                )
            )
            continue
        target = Path(target_text)
        if not target.is_absolute():
            target = full.parent / target
        try:
            exists = target.exists()
        except OSError:
            exists = False
        if not exists:
            findings.append(
                Finding(
                    code="dangling-symlink",
                    path=path,
                    detail=f"symbolic link target does not exist: {target_text}",
                )
            )
    return findings


def _symlink_target_text(path: Path) -> str | None:
    try:
        if path.is_symlink():
            return os.fsdecode(os.readlink(path))
    except OSError:
        return None
    if not path.is_file():
        return None
    try:
        raw = path.read_bytes()
    except OSError:
        return None
    if raw.endswith(b"\r\n"):
        raw = raw[:-2]
    elif raw.endswith((b"\n", b"\r")):
        raw = raw[:-1]
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return os.fsdecode(raw)


def _markdown_destinations(text: str) -> list[str]:
    text = _without_fenced_code(text.replace("\r\n", "\n").replace("\r", "\n"))
    link_labels: set[str] = set()
    image_labels: set[str] = set()
    destinations = list(_inline_destinations(text, link_labels, image_labels))
    for label, dest in _reference_definitions(text):
        key = _normalize_label(label)
        if dest and not (key in image_labels and key not in link_labels):
            destinations.append(dest)
    return destinations


def _without_fenced_code(text: str) -> str:
    kept: list[str] = []
    fence_char = ""
    fence_len = 0
    for line in text.splitlines(keepends=True):
        body = line.splitlines()[0] if line else ""
        if not fence_char:
            match = _OPEN_FENCE.match(body)
            if match:
                marker = match.group(2)
                fence_char = marker[0]
                fence_len = len(marker)
                kept.append("\n" if line.endswith("\n") else "")
                continue
            kept.append(line)
            continue
        if re.match(rf"^ {{0,3}}{re.escape(fence_char)}{{{fence_len},}}\s*$", body):
            fence_char = ""
            fence_len = 0
        kept.append("\n" if line.endswith("\n") else "")
    return "".join(kept)


def _inline_destinations(text: str, link_labels: set[str], image_labels: set[str]):
    i = 0
    n = len(text)
    while i < n:
        char = text[i]
        if char == "`":
            i = _skip_code_span(text, i)
            continue
        is_image = char == "!" and i + 1 < n and text[i + 1] == "["
        if is_image or char == "[":
            open_idx = i + 1 if is_image else i
            close = _matching_bracket(text, open_idx)
            if close is None:
                i += 1
                continue
            inner = text[open_idx + 1 : close]
            inline_dest, ref_label, end = _parse_link_suffix(text, close + 1)
            if inline_dest and not is_image:
                yield inline_dest
            elif ref_label is not None:
                key = _normalize_label(ref_label or inner)
                if key:
                    (image_labels if is_image else link_labels).add(key)
            i = end
            continue
        i += 1


def _reference_definitions(text: str):
    for match in _REF_DEF.finditer(text):
        yield match.group(1), match.group(2)


def _parse_link_suffix(text: str, index: int) -> tuple[str | None, str | None, int]:
    if index < len(text) and text[index] == "(":
        dest, end = _paren_destination(text, index)
        return dest, None, end
    if index < len(text) and text[index] == "[":
        close = _matching_bracket(text, index)
        if close is None:
            return None, None, index
        return None, text[index + 1 : close], close + 1
    return None, None, index


def _normalize_label(label: str) -> str:
    return " ".join(label.split()).casefold()


def _paren_destination(text: str, open_paren: int) -> tuple[str | None, int]:
    i = open_paren + 1
    n = len(text)
    while i < n and text[i] in " \t\n":
        i += 1
    if i >= n:
        return None, open_paren + 1
    if text[i] == "<":
        j = i + 1
        while j < n and text[j] != ">":
            j += 2 if text[j] == "\\" and j + 1 < n else 1
        if j >= n:
            return None, i
        dest = text[i + 1 : j]
        j += 1
    else:
        j = i
        depth = 0
        while j < n:
            char = text[j]
            if char == "\\" and j + 1 < n:
                j += 2
                continue
            if char in " \t\n":
                break
            if char == "(":
                depth += 1
            elif char == ")":
                if depth == 0:
                    break
                depth -= 1
            j += 1
        dest = text[i:j]
    while j < n and text[j] != ")":
        j += 1
    if j < n:
        j += 1
    dest = dest.strip()
    return (dest or None), j


def _matching_bracket(text: str, open_idx: int) -> int | None:
    depth = 0
    i = open_idx
    n = len(text)
    while i < n:
        char = text[i]
        if char == "`":
            i = _skip_code_span(text, i)
            continue
        if char == "\\" and i + 1 < n:
            i += 2
            continue
        if char == "[":
            depth += 1
        elif char == "]":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return None


def _skip_code_span(text: str, index: int) -> int:
    ticks = 0
    while index + ticks < len(text) and text[index + ticks] == "`":
        ticks += 1
    close = text.find("`" * ticks, index + ticks)
    if close < 0:
        return index + ticks
    return close + ticks


def _missing_local_target(source: Path, dest: str) -> bool:
    local = _local_path_from_destination(dest)
    if local is None:
        return False
    candidate = Path(local)
    if not candidate.is_absolute():
        candidate = source.parent / candidate
    try:
        return not candidate.exists()
    except OSError:
        return True


def _local_path_from_destination(dest: str) -> str | None:
    dest = dest.strip()
    if dest.startswith("<") and dest.endswith(">") and len(dest) >= 2:
        dest = dest[1:-1].strip()
    if not dest or dest.startswith("//"):
        return None
    parsed = urllib.parse.urlsplit(dest)
    if parsed.scheme or parsed.netloc:
        return None
    path = urllib.parse.unquote(parsed.path)
    if not path:
        return None
    return path


if __name__ == "__main__":
    raise SystemExit(main())
