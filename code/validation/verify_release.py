# -*- coding: utf-8 -*-
"""Independent release read-back, coverage, hash and optional metric verification.

Quick mode reads only the release and its validation tables. --full additionally
recomputes direct and station-grouped CV headline metrics with the formal V2.1
kernel, using prepared intermediates in --work-dir, and asserts BASE equivalence.
"""
from pathlib import Path
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import sys

import netCDF4
import numpy as np


def open_nc(path):
    with open(path, "rb") as fh:
        return netCDF4.Dataset(os.path.basename(path), memory=fh.read())


def default_release():
    if os.environ.get("SM_RELEASE_DIR"):
        return Path(os.environ["SM_RELEASE_DIR"])
    for p in (Path.cwd(), *Path(__file__).resolve().parents):
        if (p / "data" / "fusion").is_dir():
            return p
        candidate = p / "datarelease-scienceDB-v2.1"
        if (candidate / "data" / "fusion").is_dir():
            return candidate
        candidate = p / "datarelease-scienceDB"
        if (candidate / "data" / "fusion").is_dir():
            return candidate
    raise FileNotFoundError("Use --release-root or SM_RELEASE_DIR")


def pct(n, total):
    return float(100.0 * n / total) if total else 0.0


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def refresh_report_hashes(manifest, release, targets):
    """Keep self-generated verification reports consistent with the manifest."""
    targets = set(targets)
    with manifest.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        rows = list(reader)
    found = set()
    for row in rows:
        relative = Path(row["path"]).as_posix()
        if relative not in targets:
            continue
        path = release / Path(relative)
        row["size_bytes"] = str(path.stat().st_size)
        row["sha256"] = file_sha256(path)
        found.add(relative)
    if found != targets:
        raise RuntimeError(f"Manifest is missing verification report rows: {sorted(targets - found)}")
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-root", type=Path, default=default_release())
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--full", action="store_true", help="recompute direct and strict-CV headline metrics")
    parser.add_argument("--full-output", type=Path)
    parser.add_argument("--skip-checksum", action="store_true",
                        help="defer manifest checking during a staged manifest rebuild")
    args = parser.parse_args()
    release = args.release_root.resolve()
    validation = release / "validation"
    validation.mkdir(exist_ok=True)
    fusion_path = release / "data" / "fusion" / "fusion_sm_36km_3596_20150331_20261001.nc"
    checks = []

    def check(name, condition, detail=""):
        checks.append({"check": name, "pass": bool(condition), "detail": str(detail)})

    with open_nc(fusion_path) as ds:
        sm = np.ma.filled(ds.variables["sm_oi"][:], np.nan)
        ae = np.ma.filled(ds.variables["Aerr"][:], np.nan)
        flag = np.asarray(ds.variables["do_oi"][:])
        time = np.asarray(ds.variables["time"][:])
        dates = np.asarray(ds.variables["time_ymd"][:]).astype(int)
        lat = np.asarray(ds.variables["lat"][:])
        lon = np.asarray(ds.variables["lon"][:])
    check("grid_shape", sm.shape == ae.shape == flag.shape, sm.shape)
    check("time_monotonic", np.all(np.diff(time) == 1))
    check("finite_range", np.all((sm[np.isfinite(sm)] >= 0) & (sm[np.isfinite(sm)] <= 1)))
    check("aerr_nonnegative", np.all(ae[np.isfinite(ae)] >= 0))
    check("aerr_matches_sm_mask", np.array_equal(np.isfinite(ae), np.isfinite(sm)))
    check("background_error_where_no_oi", np.allclose(ae[np.isfinite(ae) & (flag == 0)], 0.046, atol=1e-7))

    release_summary = json.loads((validation / "release_fusion_summary.json").read_text(encoding="utf-8"))
    noncommon = set(int(x) for x in release_summary["noncommon_fallback_dates"])
    common = ~np.isin(dates, list(noncommon))
    check("common_analysis_day_count", int(common.sum()) == release_summary["common_analysis_days"], int(common.sum()))
    coverage_rows = []
    coverage = {}
    analysis_scope = (f"analysis_{int(common.sum())}_days_complete" if not noncommon else
                      f"analysis_{int(common.sum())}_days_excluding_{len(noncommon)}_fallbacks")
    scopes = [(f"full_{len(dates)}_days", slice(None)), (analysis_scope, common)]
    for scope, sel in scopes:
        sf, af, ff = sm[sel], ae[sel], flag[sel]
        oi = ff > 0
        finite_background = (~oi) & np.isfinite(sf)
        missing = ~np.isfinite(sf)
        total = int(sf.size)
        row = {
            "scope": scope, "n_total": total, "n_oi": int(oi.sum()),
            "n_finite_background": int(finite_background.sum()), "n_missing": int(missing.sum()),
            "pct_oi": pct(int(oi.sum()), total),
            "pct_finite_background": pct(int(finite_background.sum()), total),
            "pct_missing": pct(int(missing.sum()), total),
            "n_finite_sm_oi": int(np.isfinite(sf).sum()), "n_finite_Aerr": int(np.isfinite(af).sum()),
        }
        coverage_rows.append(row)
        coverage[scope] = row
        check(scope + "_partition", row["n_oi"] + row["n_finite_background"] + row["n_missing"] == total)
    with (validation / "coverage_summary.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(coverage_rows[0]))
        writer.writeheader()
        writer.writerows(coverage_rows)

    # Verify the complete station calendar and its CF time independently.
    station_dir = release / "data" / "station"
    with open_nc(station_dir / "staD10_qc_3596_20150101_20261001.nc") as ds:
        station_time = np.asarray(ds.variables["time"][:], dtype=float)
        station_dates = np.asarray(ds.variables["time_ymd"][:]).astype(int)
    expected_dates = np.array([
        int((dt.date(2015, 1, 1) + dt.timedelta(days=i)).strftime("%Y%m%d"))
        for i in range((dt.date(2026, 10, 1) - dt.date(2015, 1, 1)).days + 1)
    ], dtype=int)
    check("station_complete_4292_day_axis", np.array_equal(station_dates, expected_dates), len(station_dates))
    check("station_time_matches_time_ymd", np.array_equal(station_time, np.arange(len(expected_dates), dtype=float)))
    raw = np.load(station_dir / "raw3596.npy", mmap_mode="r")
    cleaned = np.load(station_dir / "cleaned3596.npy", mmap_mode="r")
    check("station_array_shapes", raw.shape == cleaned.shape == (3596, len(expected_dates)), raw.shape)

    special = {}
    for text in ("2023-01-01", "2024-01-01"):
        date = int(text.replace("-", ""))
        hit = np.flatnonzero(dates == date)
        station_hit = np.flatnonzero(station_dates == date)
        check(text + "_present", len(hit) == 1 and len(station_hit) == 1)
        if len(hit) != 1 or len(station_hit) != 1:
            continue
        i, si = int(hit[0]), int(station_hit[0])
        raw_count = int(np.isfinite(raw[:, si]).sum())
        qc_count = int(np.isfinite(cleaned[:, si]).sum())
        special[text] = {
            "raw_reporting_stations": raw_count, "qc_valid_stations": qc_count,
            "do_oi_sum": int(flag[i].sum()), "n_finite_sm_oi": int(np.isfinite(sm[i]).sum()),
            "explanation": "The restored complete station calendar includes this date; valid station observations participate in the dynamic-R OI analysis.",
        }
        check(text + "_restored_station_data", raw_count > 0 and qc_count > 0, special[text])
        check(text + "_actual_oi", int(flag[i].sum()) > 0, int(flag[i].sum()))

    hash_report = {"n_entries": 0, "n_missing": 0, "n_size_mismatch": 0,
                   "n_mismatch": 0, "n_unlisted": 0}
    manifest = validation / "SHA256SUMS.csv"
    if manifest.is_file() and not args.skip_checksum:
        listed = set()
        with manifest.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                expected, relative = row["sha256"], row["path"]
                listed.add(Path(relative).as_posix())
                hash_report["n_entries"] += 1
                path = release / Path(relative)
                if not path.is_file():
                    hash_report["n_missing"] += 1
                    continue
                if path.stat().st_size != int(row["size_bytes"]):
                    hash_report["n_size_mismatch"] += 1
                if file_sha256(path).lower() != expected.lower():
                    hash_report["n_mismatch"] += 1
        actual = {p.relative_to(release).as_posix() for p in release.rglob("*")
                  if p.is_file() and p != manifest}
        unlisted = sorted(actual - listed)
        hash_report["n_unlisted"] = len(unlisted)
        check("checksum_manifest_complete", not unlisted, unlisted[:20])
        check("checksum_manifest", hash_report["n_missing"] == 0 and
              hash_report["n_size_mismatch"] == 0 and hash_report["n_mismatch"] == 0, hash_report)
    elif not args.skip_checksum:
        check("checksum_manifest_present", False)

    headline = {"mode": "supplied_tables_checked"}
    for name in ("cv_station_grouped_summary.json", "direct_gridcell_metrics.csv"):
        check(name + "_present", (validation / name).is_file())
    if args.full:
        if not args.work_dir:
            parser.error("--full requires --work-dir")
        os.environ["SM_DATA_DIR"] = str(release / "data")
        os.environ["SM_WORK_DIR"] = str(args.work_dir.resolve())
        code_root = Path(__file__).resolve().parents[1]
        sys.path.insert(0, str(code_root / "python" / "full"))
        sys.path.insert(0, str(code_root / "validation"))
        from s15_sensitivity import Sensitivity
        output = args.full_output or (args.work_dir / "verification_headline_v7")
        engine = Sensitivity(output)
        engine.run("BASE")
        engine.verify_base()
        headline = {"mode": "recomputed", "output": str(output), "base_matches_released_metrics": True}

    passed = all(c["pass"] for c in checks)
    result = {
        "status": "PASS" if passed else "FAIL", "release": str(release),
        "checks": checks, "coverage": coverage, "special_dates": special,
        "hash_verification": hash_report, "headline_metrics": headline,
    }
    json_name = "release_verification_full.json" if args.full else "release_verification.json"
    text_name = "verify_report_full.txt" if args.full else "verify_report.txt"
    with (validation / json_name).open("w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)
    lines = ["Release verification: " + result["status"], ""]
    lines += [("%s: %s %s" % ("PASS" if c["pass"] else "FAIL", c["check"], c["detail"])).rstrip()
              for c in checks]
    lines += ["", "Coverage (full axis / common analysis axis):"]
    for row in coverage_rows:
        lines.append("%s: OI=%d (%.6f%%), finite background=%d (%.6f%%), missing=%d (%.6f%%)" %
                     (row["scope"], row["n_oi"], row["pct_oi"], row["n_finite_background"],
                      row["pct_finite_background"], row["n_missing"], row["pct_missing"]))
    (validation / text_name).write_text("\n".join(lines) + "\n", encoding="utf-8")
    if manifest.is_file() and not args.skip_checksum:
        refresh_report_hashes(manifest, release, {f"validation/{json_name}", f"validation/{text_name}"})
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
