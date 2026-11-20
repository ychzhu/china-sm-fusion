# -*- coding: utf-8 -*-
"""Build the 3596-station raw matrix on the complete 4292-day calendar."""
from paths_full import *

import glob
import json
import os
import numpy as np
import pandas as pd
import netCDF4


# Complete declared calendar; never infer it from a potentially incomplete NC.
axis_ymd = station_axis_ymd().tolist()
n_days = len(axis_ymd)
day_pos = {d: k for k, d in enumerate(axis_ymd)}
print("declared complete day axis:", n_days, axis_ymd[0], "->", axis_ymd[-1])

# gather stations across the 12 yearly CSVs
source_dir = os.environ.get("SM_STATION_SOURCE_DIR", os.path.join(STATION, "staD10_yearly"))
files = sorted(glob.glob(os.path.join(source_dir, "staD10_2*.csv")))
if not files:
    raise FileNotFoundError(
        f"No annual station CSVs in {source_dir}; set SM_STATION_SOURCE_DIR to the authoritative source tables")
meta = {}
raw = None
for fi, f in enumerate(files):
    d = pd.read_csv(f, na_values=["--", ""])
    dcols = [c for c in d.columns if c.startswith("D")]
    ids = d["Station_ID_d"].astype(int).to_numpy()
    for st_id, la, lo, pr in zip(ids, d.Lat.to_numpy(float), d.Lon.to_numpy(float),
                                 d.Province.astype(str).to_numpy()):
        if st_id not in meta:
            meta[st_id] = (float(la), float(lo), pr)
    if raw is None:
        all_ids = sorted(meta.keys() | set(ids.tolist()))
        # rebuild full id list only at end; allocate after first pass below
    print("read %s rows=%d days=%d" % (os.path.basename(f), len(d), len(dcols)))

# canonical station order: sorted by cma id
all_ids = sorted(meta.keys())
n_sta = len(all_ids)
sta_pos = {s: k for k, s in enumerate(all_ids)}
lat = np.array([meta[s][0] for s in all_ids])
lon = np.array([meta[s][1] for s in all_ids])
prov = np.array([meta[s][2] for s in all_ids])
print("unique stations:", n_sta)

raw = np.full((n_sta, n_days), np.nan, dtype=np.float32)
present_year = np.zeros((n_sta, len(files)), dtype=bool)
source_audit = []
for yi, f in enumerate(files):
    d = pd.read_csv(f, na_values=["--", ""])
    dcols = [c for c in d.columns if c.startswith("D")]
    source_dates = [c[1:] for c in dcols]
    if len(source_dates) != len(set(source_dates)):
        raise ValueError(f"Duplicate daily columns in {f}")
    years = sorted({x[:4] for x in source_dates})
    expected = {x for x in axis_ymd if x[:4] in years}
    unknown = sorted(set(source_dates) - set(day_pos))
    missing = sorted(expected - set(source_dates))
    if unknown or missing:
        raise ValueError(
            f"Annual source calendar mismatch in {os.path.basename(f)}: "
            f"unknown={unknown[:10]}, missing={missing[:10]}")
    ids = d["Station_ID_d"].astype(int).to_numpy()
    X = d[dcols].to_numpy(float)
    col_idx = np.array([day_pos[x] for x in source_dates], dtype=int)
    rows = np.array([sta_pos[s] for s in ids])
    present_year[rows, yi] = True
    # Fill raw only where a source observation is finite.  Every source date
    # has already been validated, so no date can be silently discarded.
    for ri, rr in enumerate(rows):
        vals = X[ri]
        fin = np.isfinite(vals)
        raw[rr, col_idx[fin]] = vals[fin]
    source_audit.append({
        "file": os.path.basename(f), "years": years,
        "n_source_days": len(source_dates), "n_expected_days": len(expected),
        "n_station_rows": int(len(ids)),
        "n_finite_source_values": int(np.isfinite(X).sum()),
        "calendar_complete": True,
    })

raw_daily = np.isfinite(raw).sum(axis=0)
pd.DataFrame({"date": axis_ymd, "raw_reporting_stations": raw_daily}).to_csv(
    os.path.join(WORK, "station_raw_daily_counts.csv"), index=False)
audit = {
    "axis_start": axis_ymd[0], "axis_end": axis_ymd[-1], "n_days": n_days,
    "calendar_rule": "inclusive daily calendar derived from declared start/end",
    "source_files": source_audit,
    "restored_dates": {
        x: int(raw_daily[day_pos[x]]) for x in ("20230101", "20240101")
    },
    "raw_finite_values": int(np.isfinite(raw).sum()),
}
with open(os.path.join(WORK, "station_axis_audit.json"), "w", encoding="utf-8") as fh:
    json.dump(audit, fh, ensure_ascii=False, indent=2)
print("restored raw reporting counts:", audit["restored_dates"], flush=True)

print("raw finite values:", np.isfinite(raw).sum())
print("raw range: %.4f .. %.4f" % (np.nanmin(raw), np.nanmax(raw)))
# sanity: number of values outside [0,0.6]
out_range = np.sum((raw < 0) | (raw > 0.6))
print("raw values outside [0,0.6]:", out_range)

np.save(os.path.join(WORK, "raw3596.npy"), raw)
np.save(os.path.join(WORK, "ids3596.npy"), np.array(all_ids))
np.save(os.path.join(WORK, "lat3596.npy"), lat)
np.save(os.path.join(WORK, "lon3596.npy"), lon)
np.save(os.path.join(WORK, "present_year.npy"), present_year)
pd.DataFrame({"cma_id": all_ids, "lat": lat, "lon": lon, "province": prov}
             ).to_csv(os.path.join(WORK, "stations3596_meta.csv"), index=False)
print("saved work_full raw artifacts")
