#!/usr/bin/env python3
"""Verify the delivered files without installing scientific dependencies."""
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
IGNORED_DIRS = {".git", ".venv", "generated", "__pycache__"}


def verify():
    manifest = json.loads((ROOT / "MANIFEST.json").read_text(encoding="utf-8"))
    expected = manifest["files"]
    actual = {
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file()
        and not any(part in IGNORED_DIRS for part in path.relative_to(ROOT).parts)
        and path.name != "MANIFEST.json"
    }
    errors = [f"Missing: {name}" for name in sorted(set(expected) - actual)]
    errors += [f"Unexpected: {name}" for name in sorted(actual - set(expected))]
    for name in sorted(actual & set(expected)):
        path = ROOT / name
        if path.is_symlink():
            errors.append(f"Symlink: {name}")
        elif hashlib.sha256(path.read_bytes()).hexdigest() != expected[name]:
            errors.append(f"SHA-256 mismatch: {name}")
    if errors:
        raise RuntimeError("\n".join(errors))
    return len(expected)


if __name__ == "__main__":
    try:
        count = verify()
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        print(f"Integrity check failed: {error}", file=sys.stderr)
        raise SystemExit(1)
    print(f"Integrity check passed: {count} files.")
