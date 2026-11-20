# -*- coding: utf-8 -*-
"""Synchronize, independently verify and checksum the ScienceDB V2.1 package."""
from paths_full import (DATA, REFERENCE, SMAP, WORK, STATION3596_NC,
                        STA_FALLBACK_ERR, SMAP_TCA_ERR, open_nc, station_axis_ymd)

import csv
import hashlib
import json
import os
import shutil
from pathlib import Path
import numpy as np

RELEASE = Path(os.environ.get("SM_RELEASE_DIR", os.path.dirname(DATA)))
VALIDATION = RELEASE / "validation"
VALIDATION.mkdir(parents=True, exist_ok=True)
FUSION_NAME = "fusion_sm_36km_3596_20150331_20261001.nc"
FUSION_PATH = RELEASE / "data" / "fusion" / FUSION_NAME

copy_files = [
    "direct_gridcell_metrics.csv", "station_scope_metrics.csv", "station_scope_counts.json",
    "station_fusion_metrics.csv", "cv_station_grouped_summary.json", "cv_fold_metrics.csv",
    "cv_fold_assignments.csv", "reference_station_metrics_by_scope.csv",
    "reference_station_metrics_joint.csv", "reference_roles.json", "tca_sensitivity.json",
    "tca_sensitivity_by_station.csv", "aerr_summary.json", "fusion_diagnostics.json",
    "release_fusion_summary.json", "station_release_summary.json", "station_network_summary.json",
    "station_axis_audit.json", "station_daily_counts.csv", "raw_axis_comparison.json",
    "station_product_verification.json", "v2_key_metrics.json", "stratification_counts.json",
    "station_annual_metrics_inside.csv", "station_seasonal_metrics_inside.csv",
    "station_monthly_metrics_inside.csv", "station_density_metrics_inside.csv",
    "station_ndvi_metrics_inside.csv", "station_province_metrics_inside.csv",
    "gridcell_annual_metrics.csv", "gridcell_monthly_metrics.csv",
    "gridcell_seasonal_metrics.csv", "gridcell_daily_metrics.csv",
    "cv_annual_metrics.csv", "cv_seasonal_metrics.csv",
    "reference_station_metrics_monthly.csv", "reference_station_metrics_seasonal.csv",
]
for name in copy_files:
    src = Path(WORK) / name
    if not src.is_file():
        raise FileNotFoundError(src)
    shutil.copy2(src, VALIDATION / name)
for name in ["stats_smap_fusion.csv", "stats_smap_fusion_pergrid.csv"]:
    shutil.copy2(Path(WORK) / name, RELEASE / "data" / "reference" / name)
for name in ["within_cell_pairs.csv", "within_cell_distance.csv", "within_cell_summary.json",
             "hl_background_correlation.csv", "hl_background_correlation_summary.json",
             "aerr_deciles.csv", "aerr_cell_metrics.csv", "aerr_reliability_summary.json"]:
    src = Path(WORK) / "diagnostics_v7" / name
    if not src.is_file():
        raise FileNotFoundError(src)
    shutil.copy2(src, VALIDATION / name)
for name in ["sensitivity_common_metrics.csv", "sensitivity_summary.json"]:
    src = Path(WORK) / "sensitivity_v7" / name
    if not src.is_file():
        raise FileNotFoundError(src)
    shutil.copy2(src, VALIDATION / name)

# Read the actual release files, not only the in-memory working arrays.
ds = open_nc(str(FUSION_PATH))
lat = np.asarray(ds.variables["lat"][:])
lon = np.asarray(ds.variables["lon"][:])
grid_index = np.asarray(ds.variables["grid_index"][:])
ymd = np.asarray(ds.variables["time_ymd"][:]).astype(int)
sm = np.ma.filled(ds.variables["sm_oi"][:], np.nan).astype(float)
ae = np.ma.filled(ds.variables["Aerr"][:], np.nan).astype(float)
flag = np.asarray(ds.variables["do_oi"][:])
ds.close()

with (VALIDATION / "release_fusion_summary.json").open(encoding="utf-8") as fh:
    expected_fusion = json.load(fh)
with (VALIDATION / "station_release_summary.json").open(encoding="utf-8") as fh:
    expected_station = json.load(fh)

checks = []
def check(name, condition, detail):
    checks.append((name, bool(condition), detail))
    if not condition:
        raise AssertionError(f"{name}: {detail}")

check("dimensions", sm.shape == (expected_fusion["days"], expected_fusion["grids"]) and ae.shape == sm.shape and flag.shape == sm.shape,
      f"sm={sm.shape}, Aerr={ae.shape}, do_oi={flag.shape}")
finite_sm = np.isfinite(sm)
finite_ae = np.isfinite(ae)
check("finite_sm", finite_sm.sum() == expected_fusion["finite_sm"], f"N={finite_sm.sum()}")
check("Aerr_availability", np.array_equal(finite_sm, finite_ae), "Aerr finite exactly where sm_oi is finite")
check("do_oi_values", set(np.unique(flag)).issubset({0, 1}), f"values={np.unique(flag).tolist()}")
check("actual_oi_grid_days", int(flag.sum()) == expected_fusion["actual_oi_grid_days"], f"N={int(flag.sum())}")
check("Aerr_background_bound", not np.any(ae[finite_ae] > SMAP_TCA_ERR + 1e-6),
      f"max={float(ae[finite_ae].max()):.9f}")
check("Aerr_positive", np.all(ae[finite_ae] >= 0), f"min={float(ae[finite_ae].min()):.9f}")
check("physical_range", np.all((sm[finite_sm] >= 0) & (sm[finite_sm] <= 0.7)),
      f"range={float(sm[finite_sm].min()):.6f}..{float(sm[finite_sm].max()):.6f}")
check("no_oi_zero_error", np.allclose(ae[(flag == 0) & finite_ae], SMAP_TCA_ERR, atol=1e-7),
      "all background-only finite cells retain SMAP error")
check("oi_error_reduction", np.all(ae[(flag == 1) & finite_ae] <= SMAP_TCA_ERR + 1e-6),
      "all OI cells have posterior error no greater than background")

# Compare every common production day against the actual working arrays.
g = np.load(Path(WORK) / "grid_prep_full.npz", allow_pickle=True)
common = g["common_ymd"].astype(int)
work_fusion = np.load(Path(WORK) / "fusion_full.npy")
work_aerr = np.load(Path(WORK) / "fusion_full_err.npy")
work_flag = np.load(Path(WORK) / "fusion_do_oi.npy")
common_rows = np.array([int(np.flatnonzero(ymd == d)[0]) for d in common])
check("common_day_count", len(common) == expected_fusion["common_analysis_days"], f"N={len(common)}")
check("release_matches_work_sm", np.allclose(sm[common_rows], np.clip(work_fusion, 0, 0.7), atol=1e-7, equal_nan=True),
      "all common-day sm_oi values match clipped working analysis")
check("release_matches_work_Aerr", np.allclose(ae[common_rows], work_aerr, atol=1e-7, equal_nan=True),
      "all common-day posterior errors match working analysis")
check("release_matches_work_flags", np.array_equal(flag[common_rows], work_flag),
      "all common-day OI flags match working analysis")
check("oi_condition_grids", int((flag.sum(axis=0) > 0).sum()) == int((work_flag.sum(axis=0) > 0).sum()),
      f"N={int((flag.sum(axis=0) > 0).sum())}")

# Verify the complete station calendar, CF time and all three released arrays.
station_days = station_axis_ymd().astype(np.int32)
sds = open_nc(str(STATION3596_NC))
station_time = np.asarray(sds.variables["time"][:], dtype=float)
station_ymd = np.asarray(sds.variables["time_ymd"][:]).astype(int)
station_sm = np.ma.filled(sds.variables["sm"][:], np.nan).astype(np.float32)
station_flag = np.asarray(sds.variables["qc_flag"][:])
sds.close()
raw_release = np.load(RELEASE / "data" / "station" / "raw3596.npy")
clean_release = np.load(RELEASE / "data" / "station" / "cleaned3596.npy")
flag_release = np.load(RELEASE / "data" / "station" / "qc_flag3596.npy")
check("station_complete_calendar", np.array_equal(station_ymd, station_days),
      f"{station_ymd[0]}..{station_ymd[-1]}, N={len(station_ymd)}")
check("station_CF_time", np.array_equal(station_time, np.arange(len(station_days), dtype=float)),
      "time equals elapsed days since 2015-01-01 and is consistent with time_ymd")
check("station_raw_shape", raw_release.shape == (3596, len(station_days)), str(raw_release.shape))
check("station_clean_shape", clean_release.shape == raw_release.shape, str(clean_release.shape))
check("station_flag_shape", flag_release.shape == raw_release.shape, str(flag_release.shape))
check("station_nc_matches_cleaned", np.array_equal(station_sm, clean_release.T, equal_nan=True),
      "NetCDF sm is the transpose of cleaned3596.npy including NaN masks")
check("station_nc_matches_flags", np.array_equal(station_flag, flag_release.T),
      "NetCDF qc_flag is the transpose of qc_flag3596.npy")
check("station_finite_cleaned", int(np.isfinite(clean_release).sum()) == expected_station["finite_cleaned"],
      f"N={int(np.isfinite(clean_release).sum())}")

# The two restored dates must contain real station data and actual OI updates.
for date in [20230101, 20240101]:
    di = int(np.flatnonzero(ymd == date)[0])
    si = int(np.flatnonzero(station_ymd == date)[0])
    raw_count = int(np.isfinite(raw_release[:, si]).sum())
    qc_count = int(np.isfinite(clean_release[:, si]).sum())
    check(f"restored_{date}_raw", raw_count > 0, f"raw stations={raw_count}")
    check(f"restored_{date}_qc", qc_count > 0, f"QC stations={qc_count}")
    check(f"restored_{date}_in_common", date in set(common.tolist()), "date enters station-SMAP analysis")
    check(f"restored_{date}_oi", int(flag[di].sum()) > 0, f"OI grid-cells={int(flag[di].sum())}")

# Per-grid coverage is derived from the verified release array.
coverage_path = RELEASE / "data" / "coverage" / "coverage_fusion.csv"
with coverage_path.open("w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["grid", "lat", "lon", "grid_index", "nvalid_days", "valid_pct"])
    for i in range(sm.shape[1]):
        n = int(finite_sm[:, i].sum())
        w.writerow([i + 1, float(lat[i]), float(lon[i]), int(grid_index[i]), n, round(100*n/sm.shape[0], 2)])

# Regenerate the release coverage preview from the same array.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8.2, 5.8), constrained_layout=True)
pct = 100 * finite_sm.sum(axis=0) / sm.shape[0]
sc = ax.scatter(lon, lat, c=pct, s=5, cmap="viridis", vmin=0, vmax=100, linewidths=0)
ax.set(xlabel="Longitude (°E)", ylabel="Latitude (°N)", xlim=(73, 136), ylim=(17, 54),
       title="V2.1 fusion valid-data coverage (2015-03-31 to 2026-10-01)")
ax.set_aspect("equal", adjustable="box")
cb = fig.colorbar(sc, ax=ax, pad=0.02)
cb.set_label("Valid days (%)")
fig.savefig(RELEASE / "data" / "coverage" / "coverage_fusion.png", dpi=300)
plt.close(fig)

summary = {
    "days": int(sm.shape[0]), "grids": int(sm.shape[1]), "finite_sm": int(finite_sm.sum()),
    "finite_Aerr": int(finite_ae.sum()), "actual_oi_grid_days": int(flag.sum()),
    "oi_condition_grids": int((flag.sum(axis=0) > 0).sum()),
    "Aerr_median": float(np.median(ae[finite_ae])), "Aerr_mean": float(np.mean(ae[finite_ae])),
    "Aerr_P10": float(np.percentile(ae[finite_ae], 10)),
    "Aerr_P90": float(np.percentile(ae[finite_ae], 90)),
    "Aerr_max": float(ae[finite_ae].max()),
    "station_days": int(len(station_days)), "station_finite_raw": int(np.isfinite(raw_release).sum()),
    "station_finite_cleaned": int(np.isfinite(clean_release).sum()),
}
with (VALIDATION / "release_readback_summary.json").open("w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=2)
with (VALIDATION / "verify_report.txt").open("w", encoding="utf-8") as fh:
    fh.write("ScienceDB V2.1 independent release read-back verification\n")
    fh.write("Generated from the actual NetCDF and NumPy products after export.\n\n")
    for name, passed, detail in checks:
        fh.write(f"PASS | {name} | {detail}\n")
    fh.write("\nSummary\n" + json.dumps(summary, indent=2) + "\n")

# The manifest is regenerated again after final diagnostics and metadata updates.
manifest = RELEASE / "validation" / "SHA256SUMS.csv"
rows = []
for path in sorted(RELEASE.rglob("*")):
    if not path.is_file() or path == manifest:
        continue
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            h.update(block)
    rows.append((path.relative_to(RELEASE).as_posix(), path.stat().st_size, h.hexdigest()))
with manifest.open("w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["path", "size_bytes", "sha256"])
    w.writerows(rows)
print(json.dumps({"release": str(RELEASE), "checks": len(checks), "manifest_files": len(rows), **summary}, indent=2))
