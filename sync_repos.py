#!/usr/bin/env python3
"""
sync_repos.py

Flattens all matching files from a "normal" (nested) repo checkout into a
"bots" (flat, no-subfolders) repo checkout.

- Walks SOURCE_DIR recursively, collects files with the target extensions.
- Copies each one into TARGET_DIR's root, using just the basename (flattened).
- Skips copying a file if the content is byte-identical to what's already
  there, so you don't get noisy "no-op" commits.
- Deletes files in TARGET_DIR that have a target extension but no longer
  correspond to any file in SOURCE_DIR (keeps the bot repo a true mirror).
- Refuses to silently overwrite when two different source files would
  flatten to the same basename (name collision) -- exits non-zero so the
  GitHub Action run fails loudly instead of quietly losing a file.

Uses only the Python standard library (os, shutil, filecmp, pathlib,
hashlib) -- no pip installs required, per project constraints.

Usage:
    python3 sync_repos.py <source_dir> <target_dir>
"""

import filecmp
import shutil
import sys
from pathlib import Path

# Extensions to mirror. Add/remove as needed.
EXTENSIONS = {".py", ".ipynb", ".pgn", ".tcl"}

# Directory names to never descend into when walking the source repo.
IGNORE_DIRS = {".git", ".github", "__pycache__", ".ipynb_checkpoints", ".venv", "venv"}

# Files that live in the target repo but are NOT managed by this sync
# (so they are never auto-deleted even if their extension matches).
PROTECTED_TARGET_NAMES = set()  # e.g. {"keep_this.py"} if you ever need an exception


def find_source_files(source_dir: Path) -> dict[str, Path]:
    """
    Walk source_dir recursively and return {basename: full_path} for every
    file whose extension is in EXTENSIONS. Raises on basename collisions.
    """
    found: dict[str, list[Path]] = {}
    for path in source_dir.rglob("*"):
        if path.is_dir():
            continue
        if any(part in IGNORE_DIRS for part in path.relative_to(source_dir).parts):
            continue
        if path.suffix.lower() in EXTENSIONS:
            found.setdefault(path.name, []).append(path)

    collisions = {name: paths for name, paths in found.items() if len(paths) > 1}
    if collisions:
        print("ERROR: filename collisions found -- flattening would overwrite files:", file=sys.stderr)
        for name, paths in collisions.items():
            print(f"  {name}:", file=sys.stderr)
            for p in paths:
                print(f"    - {p.relative_to(source_dir)}", file=sys.stderr)
        print("\nRename one of the conflicting files in the source repo, or", file=sys.stderr)
        print("extend sync_repos.py to namespace them, then re-run.", file=sys.stderr)
        sys.exit(1)

    return {name: paths[0] for name, paths in found.items()}


def sync(source_dir: Path, target_dir: Path) -> bool:
    """Returns True if any files were added, updated, or removed."""
    source_files = find_source_files(source_dir)
    changed = False

    added, updated, removed = [], [], []

    # Copy / update
    for name, src_path in source_files.items():
        dst_path = target_dir / name
        if dst_path.exists() and filecmp.cmp(src_path, dst_path, shallow=False):
            continue  # identical content, nothing to do
        is_new = not dst_path.exists()
        shutil.copy2(src_path, dst_path)
        (added if is_new else updated).append(name)
        changed = True

    # Remove target files that no longer have a source counterpart
    expected_names = set(source_files.keys())
    for path in target_dir.iterdir():
        if path.is_dir():
            continue
        if path.name in PROTECTED_TARGET_NAMES:
            continue
        if path.suffix.lower() in EXTENSIONS and path.name not in expected_names:
            path.unlink()
            removed.append(path.name)
            changed = True

    if added:
        print(f"Added ({len(added)}):")
        for n in sorted(added):
            print(f"  + {n}")
    if updated:
        print(f"Updated ({len(updated)}):")
        for n in sorted(updated):
            print(f"  ~ {n}")
    if removed:
        print(f"Removed ({len(removed)}):")
        for n in sorted(removed):
            print(f"  - {n}")
    if not changed:
        print("No changes -- source and target are already in sync.")

    return changed


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 sync_repos.py <source_dir> <target_dir>", file=sys.stderr)
        sys.exit(2)

    source_dir = Path(sys.argv[1]).resolve()
    target_dir = Path(sys.argv[2]).resolve()

    if not source_dir.is_dir():
        print(f"ERROR: source dir does not exist: {source_dir}", file=sys.stderr)
        sys.exit(2)
    if not target_dir.is_dir():
        print(f"ERROR: target dir does not exist: {target_dir}", file=sys.stderr)
        sys.exit(2)

    changed = sync(source_dir, target_dir)

    # Emit a flag file the calling shell/workflow can check cheaply,
    # avoiding the need to parse stdout.
    (target_dir / ".sync_changed").write_text("1" if changed else "0")


if __name__ == "__main__":
    main()
