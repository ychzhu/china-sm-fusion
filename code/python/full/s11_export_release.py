# -*- coding: utf-8 -*-
"""Export V2.1 fusion NetCDF with posterior Aerr and actual daily OI flags."""
from paths_full import *
from oi_utils import build_geometry

import datetime
import json
import os
from pathlib import Path
import shutil
import tempfile
import numpy as np
import netCDF4

A = 10.0
THRESH_NUM, MIN_HIST = 3, 30
source_path = Path(FUSION3596_NC)
backup_dir = Path(WORK).parent / "v2_baseline"
backup_dir.mkdir(parents=True, exist_ok=True)
backup_path = backup_dir / source_path.name
if not backup_path.exists():
    shutil.copyfile(source_path, backup_path)

source = open_nc(str(source_path))
time_ymd = source.variables["time_ymd"][:].astype(np.int32)
base_variables = {}
for name in ["lat", "lon", "grid_index", "time", "time_ymd"]:
    var = source.variables[name]
    base_variables[name] = {
        "data": var[:], "dims": var.dimensions, "dtype": var.dtype,
        "attrs": {a: var.getncattr(a) for a in var.ncattrs() if a != "_FillValue"},
        "fill": getattr(var, "_FillValue", None),
    }
global_attrs = {a: source.getncattr(a) for a in source.ncattrs()}
source.close()

g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
common = g["common_ymd"].astype(str)
common_int = common.astype(np.int32)
fusion = np.load(os.path.join(WORK, "fusion_full.npy"))
aerr = np.load(os.path.join(WORK, "fusion_full_err.npy"))
assert fusion.shape == aerr.shape == (len(common), len(base_variables["lat"]["data"]))

n_time, n_grid = len(time_ymd), fusion.shape[1]
sm_full = np.full((n_time, n_grid), np.nan, dtype=np.float32)
aerr_full = np.full((n_time, n_grid), np.nan, dtype=np.float32)
do_full = np.zeros((n_time, n_grid), dtype=np.int8)
common_rows = np.array([int(np.where(time_ymd == d)[0][0]) for d in common_int])
sm_full[common_rows] = np.clip(fusion, 0.0, 0.7).astype(np.float32)
aerr_full[common_rows] = aerr.astype(np.float32)

# Dates outside the common station/SMAP analysis axis, if any, retain raw SMAP.
# The rebuilt V2.1 calendar has no non-common dates, but the generic policy is
# retained for future temporal extensions rather than hard-coding two dates.
noncommon = np.setdiff1d(np.arange(n_time), common_rows)
if len(noncommon):
    smap = open_nc(os.path.join(SMAP, "smapOI_2day_20150331_20261017.nc"))
    smap_dates = smap.variables["time_ymd"][:].astype(np.int32)
    smap_values = np.ma.filled(smap.variables["sm"][:], np.nan)
    smap.close()
    pos = {int(d): i for i, d in enumerate(smap_dates)}
    for row in noncommon:
        d = int(time_ymd[row])
        if d in pos:
            sm_full[row] = np.clip(smap_values[pos[d]], 0.0, 0.7).astype(np.float32)
            aerr_full[row] = np.where(np.isfinite(sm_full[row]), 0.046, np.nan).astype(np.float32)

flag_path = os.path.join(WORK, "fusion_do_oi.npy")
if os.path.exists(flag_path):
    do_full[common_rows] = np.load(flag_path).astype(np.int8)
else:
    # Reconstruct actual information-use flags without recomputing OI weights.
    nb = np.load(os.path.join(WORK, "neighbors_full.npz"), allow_pickle=True)
    smap_cdf = np.load(os.path.join(WORK, "smap_cdf_full.npy"))
    bearing = g["bearing"]
    glat, glon = g["glat"], g["glon"]
    grid_obs = g["grid_obs"]
    grid_err_daily = g["grid_err_daily"]
    active, neighbors, distances = [], [], []
    for gi in range(n_grid):
        if nb["num_around"][gi] < THRESH_NUM or np.isfinite(smap_cdf[:, gi]).sum() <= MIN_HIST:
            continue
        bp = nb["cap_idx"][gi]
        bp = bp[bp >= 0]
        active.append(gi)
        neighbors.append(bp)
        distances.append(nb["cap_dist"][gi][:len(bp)])
    active = np.asarray(active, dtype=int)
    bp_pad, _, _ = build_geometry(glat, glon, bearing, active, neighbors, distances, A)
    bp_safe = np.where(bp_pad >= 0, bp_pad, 0)
    valid_pad = bp_pad >= 0
    for di in range(len(common)):
        target = smap_cdf[di, active]
        usable = (valid_pad & np.isfinite(grid_obs[bp_safe, di]) &
                  np.isfinite(smap_cdf[di, bearing[bp_safe]]) &
                  np.isfinite(grid_err_daily[bp_safe, di]) &
                  np.isfinite(target)[:, None])
        do_full[common_rows[di], active] = np.any(usable, axis=1).astype(np.int8)

fd, temp_path = tempfile.mkstemp(suffix=".nc", prefix="china_sm_v2_")
os.close(fd)
os.unlink(temp_path)
try:
    out = netCDF4.Dataset(temp_path, "w", format="NETCDF4")
    out.createDimension("grid", n_grid)
    out.createDimension("time", n_time)
    for name, spec in base_variables.items():
        create = {"zlib": True, "complevel": 4, "shuffle": True}
        if spec["fill"] is not None:
            create["fill_value"] = spec["fill"]
        var = out.createVariable(name, spec["dtype"], spec["dims"], **create)
        for attr, value in spec["attrs"].items():
            var.setncattr(attr, value)
        var[:] = spec["data"]
    fill = np.float32(-9999.0)
    sm = out.createVariable("sm_oi", "f4", ("time", "grid"), zlib=True,
                            complevel=4, shuffle=True, fill_value=fill,
                            chunksizes=(30, n_grid))
    sm.long_name = "daily surface soil moisture fused by dynamic-error optimal interpolation"
    sm.units = "m3 m-3"
    sm.valid_range = np.array([0.0, 0.7], dtype=np.float32)
    sm.comment = "CDF-adjusted SMAP plus OI where daily station information is available; non-updated cells retain the available background"
    ae = out.createVariable("Aerr", "f4", ("time", "grid"), zlib=True,
                            complevel=4, shuffle=True, fill_value=fill,
                            chunksizes=(30, n_grid))
    ae.long_name = "standard optimal-interpolation posterior analysis error standard deviation"
    ae.units = "m3 m-3"
    ae.comment = "sigma_b*sqrt(1-c^T w); R uses actual valid stations per cell-day; 0.046 background SD where no OI information is used"
    do = out.createVariable("do_oi", "i1", ("time", "grid"), zlib=True,
                            complevel=4, shuffle=True, chunksizes=(30, n_grid))
    do.long_name = "actual daily optimal-interpolation application flag"
    do.flag_values = np.array([0, 1], dtype=np.int8)
    do.flag_meanings = "background_only daily_station_information_applied"
    for attr, value in global_attrs.items():
        out.setncattr(attr, value)
    out.setncattr("title", f"China 36 km daily surface soil moisture fusion dataset, Version {RELEASE_VERSION}")
    out.setncattr("version", RELEASE_VERSION)
    out.setncattr("oi_model", "dynamic daily observation R and standard posterior Aerr; a=10; R=1.08 degrees; MAX_NEAR=9")
    out.setncattr("analysis_axis", f"{len(common)} common station-SMAP days; full product axis {n_time} days")
    out.setncattr("station_background_error_fallback", STA_FALLBACK_ERR)
    out.setncattr("station_fallback_basis", "three-decimal median of valid station TCA error estimates")
    if len(noncommon):
        out.setncattr("noncommon_date_policy", "raw SMAP 2-day background, Aerr=0.046 where finite, do_oi=0")
    else:
        out.setncattr("noncommon_date_policy", "none; all product dates belong to the common station-SMAP analysis axis")
    out.setncattr("date_created", datetime.datetime.now(datetime.timezone.utc).isoformat())
    sm[:] = np.ma.masked_invalid(sm_full)
    ae[:] = np.ma.masked_invalid(aerr_full)
    do[:] = do_full
    out.close()
    shutil.copyfile(temp_path, source_path)
finally:
    if os.path.exists(temp_path):
        os.unlink(temp_path)

valid_sm = np.isfinite(sm_full)
valid_ae = np.isfinite(aerr_full)
summary = {
    "file": str(source_path), "version": RELEASE_VERSION, "days": int(n_time),
    "common_analysis_days": int(len(common)), "noncommon_fallback_dates": time_ymd[noncommon].astype(int).tolist(),
    "grids": int(n_grid), "finite_sm": int(valid_sm.sum()),
    "actual_oi_grid_days": int(do_full.sum()),
    "Aerr_median": float(np.median(aerr_full[valid_ae])),
    "Aerr_mean": float(np.mean(aerr_full[valid_ae])),
    "Aerr_P10": float(np.percentile(aerr_full[valid_ae], 10)),
    "Aerr_P90": float(np.percentile(aerr_full[valid_ae], 90)),
    "Aerr_above_background": int(np.sum(aerr_full[valid_ae] > 0.046 + 1e-6)),
}
with open(os.path.join(WORK, "release_fusion_summary.json"), "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)
print(json.dumps(summary, ensure_ascii=False, indent=2))
