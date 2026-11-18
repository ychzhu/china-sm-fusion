# -*- coding: utf-8 -*-
"""Synchronize, independently verify and checksum the ScienceDB V2 package."""
from paths_full import DATA, REFERENCE, SMAP, WORK, open_nc

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
    "release_fusion_summary.json", "v2_key_metrics.json", "stratification_counts.json",
    "station_annual_metrics_inside.csv", "station_seasonal_metrics_inside.csv",
    "station_monthly_metrics_inside.csv", "station_density_metrics_inside.csv",
    "station_ndvi_metrics_inside.csv", "station_province_metrics_inside.csv",
    "gridcell_annual_metrics.csv", "gridcell_monthly_metrics.csv",
    "gridcell_seasonal_metrics.csv", "gridcell_daily_metrics.csv",
    "cv_annual_metrics.csv", "cv_seasonal_metrics.csv",
]
for name in copy_files:
    src = Path(WORK) / name
    if not src.is_file():
        raise FileNotFoundError(src)
    shutil.copy2(src, VALIDATION / name)
for name in ["stats_smap_fusion.csv", "stats_smap_fusion_pergrid.csv"]:
    shutil.copy2(Path(WORK) / name, RELEASE / "data" / "reference" / name)

# Read the actual release file, not the in-memory working arrays.
ds = open_nc(str(FUSION_PATH))
lat = np.asarray(ds.variables["lat"][:])
lon = np.asarray(ds.variables["lon"][:])
grid_index = np.asarray(ds.variables["grid_index"][:])
ymd = np.asarray(ds.variables["time_ymd"][:]).astype(int)
sm = np.ma.filled(ds.variables["sm_oi"][:], np.nan).astype(float)
ae = np.ma.filled(ds.variables["Aerr"][:], np.nan).astype(float)
flag = np.asarray(ds.variables["do_oi"][:])
fill_sm = float(ds.variables["sm_oi"]._FillValue)
ds.close()

checks = []
def check(name, condition, detail):
    checks.append((name, bool(condition), detail))
    if not condition:
        raise AssertionError(f"{name}: {detail}")

check("dimensions", sm.shape == (4116, 7313) and ae.shape == sm.shape and flag.shape == sm.shape,
      f"sm={sm.shape}, Aerr={ae.shape}, do_oi={flag.shape}")
finite_sm = np.isfinite(sm)
finite_ae = np.isfinite(ae)
check("finite_sm", finite_sm.sum() == 21084332, f"N={finite_sm.sum()}")
check("Aerr_availability", np.array_equal(finite_sm, finite_ae), "Aerr finite exactly where sm_oi is finite")
check("do_oi_values", set(np.unique(flag)).issubset({0, 1}), f"values={np.unique(flag).tolist()}")
check("actual_oi_grid_days", int(flag.sum()) == 14017073, f"N={int(flag.sum())}")
check("oi_condition_grids", int((flag.sum(axis=0) > 0).sum()) == 4689,
      f"N={int((flag.sum(axis=0) > 0).sum())}")
check("Aerr_background_bound", not np.any(ae[finite_ae] > 0.046 + 1e-6),
      f"max={float(ae[finite_ae].max()):.9f}")
check("Aerr_positive", np.all(ae[finite_ae] >= 0), f"min={float(ae[finite_ae].min()):.9f}")
check("physical_range", np.all((sm[finite_sm] >= 0) & (sm[finite_sm] <= 0.7)),
      f"range={float(sm[finite_sm].min()):.6f}..{float(sm[finite_sm].max()):.6f}")
check("no_oi_zero_error", np.allclose(ae[(flag == 0) & finite_ae], 0.046, atol=1e-7),
      "all background-only finite cells retain SMAP error")
check("oi_error_reduction", np.all(ae[(flag == 1) & finite_ae] <= 0.046 + 1e-6),
      "all OI cells have posterior error no greater than background")

# Verify both non-common dates against the actual raw SMAP two-day background.
smap_path = Path(SMAP) / "smapOI_2day_20150331_20261017.nc"
raw_ds = open_nc(str(smap_path))
raw_ymd = np.asarray(raw_ds.variables["time_ymd"][:]).astype(int)
raw_sm = np.ma.filled(raw_ds.variables["sm"][:], np.nan).astype(float)
raw_ds.close()
raw_pos = {d: i for i, d in enumerate(raw_ymd)}
for date in [20230101, 20240101]:
    di = int(np.where(ymd == date)[0][0])
    raw = np.clip(raw_sm[raw_pos[date]], 0, 0.7)
    valid = np.isfinite(raw)
    check(f"fallback_{date}_flag", not flag[di].any(), f"sum={int(flag[di].sum())}")
    check(f"fallback_{date}_sm", np.array_equal(sm[di][valid], raw[valid]),
          f"finite raw={int(valid.sum())}")
    check(f"fallback_{date}_aerr", np.allclose(ae[di][valid], 0.046, atol=1e-7),
          "raw-background finite cells use 0.046")
    check(f"fallback_{date}_missing", not finite_sm[di][~valid].any() and not finite_ae[di][~valid].any(),
          "missing raw background remains missing")

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
       title="V2 fusion valid-data coverage (2015-03-31 to 2026-10-01)")
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
}
with (VALIDATION / "release_readback_summary.json").open("w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=2)
with (VALIDATION / "verify_report.txt").open("w", encoding="utf-8") as fh:
    fh.write("ScienceDB V2 independent release read-back verification\n")
    fh.write("Generated from the actual NetCDF after export.\n\n")
    for name, passed, detail in checks:
        fh.write(f"PASS | {name} | {detail}\n")
    fh.write("\nSummary\n" + json.dumps(summary, indent=2) + "\n")

# The manifest is generated last so it covers every final document and table.
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
