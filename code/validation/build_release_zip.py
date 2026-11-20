#!/usr/bin/env python3
"""Build and independently read back the final ScienceDB release ZIP."""
import argparse
import csv
import hashlib
import json
from pathlib import Path, PurePosixPath
import zipfile


def hash_stream(handle):
    digest = hashlib.sha256()
    for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
        digest.update(block)
    return digest.hexdigest()


def hash_file(path):
    with path.open("rb") as handle:
        return hash_stream(handle)


def load_manifest(release):
    with (release / "validation" / "SHA256SUMS.csv").open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def verify_manifest(release, rows):
    listed = {PurePosixPath(row["path"]).as_posix() for row in rows}
    actual = {path.relative_to(release).as_posix() for path in release.rglob("*")
              if path.is_file() and path != release / "validation" / "SHA256SUMS.csv"}
    errors = []
    if listed != actual:
        errors.append({"manifest_set_difference": {"unlisted": sorted(actual - listed),
                                                    "missing_on_disk": sorted(listed - actual)}})
    for row in rows:
        path = release / PurePosixPath(row["path"])
        if not path.is_file():
            continue
        if path.stat().st_size != int(row["size_bytes"]) or hash_file(path) != row["sha256"]:
            errors.append({"source_mismatch": row["path"]})
    return actual, errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    args = parser.parse_args()
    release = args.release_root.resolve()
    output = args.output.resolve()
    audit_output = args.audit_output.resolve()
    rows = load_manifest(release)
    source_files, errors = verify_manifest(release, rows)

    quick = json.loads((release / "validation" / "release_verification.json").read_text(encoding="utf-8"))
    full = json.loads((release / "validation" / "release_verification_full.json").read_text(encoding="utf-8"))
    if quick.get("status") != "PASS" or quick.get("hash_verification", {}).get("n_entries") != len(rows):
        errors.append("Final quick verification is not a complete PASS")
    if full.get("status") != "PASS" or full.get("headline_metrics", {}).get("mode") != "recomputed":
        errors.append("Saved full verification is not a recomputed PASS")
    if errors:
        raise RuntimeError(json.dumps(errors, ensure_ascii=False))

    output.parent.mkdir(parents=True, exist_ok=True)
    all_files = sorted(path for path in release.rglob("*") if path.is_file())
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED,
                         compresslevel=1, allowZip64=True) as archive:
        for path in all_files:
            archive.write(path, path.relative_to(release).as_posix())

    zip_errors = []
    with zipfile.ZipFile(output, "r") as archive:
        bad_member = archive.testzip()
        if bad_member:
            zip_errors.append({"crc_failure": bad_member})
        names = archive.namelist()
        expected_names = {path.relative_to(release).as_posix() for path in all_files}
        if set(names) != expected_names or len(names) != len(expected_names):
            zip_errors.append("ZIP member set differs from release")
        for name in names:
            pure = PurePosixPath(name)
            if pure.is_absolute() or ".." in pure.parts:
                zip_errors.append({"unsafe_member": name})
        for row in rows:
            name = PurePosixPath(row["path"]).as_posix()
            with archive.open(name) as handle:
                if hash_stream(handle) != row["sha256"]:
                    zip_errors.append({"zip_hash_mismatch": name})
        with archive.open("validation/SHA256SUMS.csv") as zipped, (release / "validation" / "SHA256SUMS.csv").open("rb") as source:
            if hash_stream(zipped) != hash_stream(source):
                zip_errors.append("ZIP manifest differs from source manifest")

    report = {
        "status": "PASS" if not zip_errors else "FAIL",
        "release_files": len(all_files),
        "release_bytes": sum(path.stat().st_size for path in all_files),
        "manifest_entries": len(rows),
        "manifest_source_recheck": "PASS",
        "quick_verification": quick["status"],
        "full_verification": full["status"],
        "full_headline_mode": full["headline_metrics"]["mode"],
        "zip_path": str(output),
        "zip_bytes": output.stat().st_size,
        "zip_sha256": hash_file(output),
        "zip_members": len(all_files),
        "zip_crc_test": "PASS" if not bad_member else "FAIL",
        "zip_manifest_hash_recheck": "PASS" if not zip_errors else "FAIL",
        "errors": zip_errors,
    }
    audit_output.parent.mkdir(parents=True, exist_ok=True)
    audit_output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if zip_errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
