#!/usr/bin/env python3
"""Rebuild the release SHA-256 manifest after all other files are final."""
import argparse
import csv
import hashlib
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-root", type=Path, required=True)
    args = parser.parse_args()
    release = args.release_root.resolve()
    manifest = release / "validation" / "SHA256SUMS.csv"
    rows = []
    for path in sorted(release.rglob("*")):
        if not path.is_file() or path == manifest:
            continue
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                digest.update(block)
        rows.append((path.relative_to(release).as_posix(), path.stat().st_size,
                     digest.hexdigest()))
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["path", "size_bytes", "sha256"])
        writer.writerows(rows)
    print(f"Wrote {len(rows)} entries to {manifest}")


if __name__ == "__main__":
    main()
