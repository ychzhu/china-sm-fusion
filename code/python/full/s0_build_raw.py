# -*- coding: utf-8 -*-
"""Build the full 3596-station raw matrix (authoritative 4290-day axis)."""
from paths_full import *

import glob
import os
import numpy as np
import pandas as pd
import netCDF4


# authoritative day axis from the released 3596 station product
axis_ymd = station_axis_ymd().tolist()
n_days = len(axis_ymd)
day_pos = {d: k for k, d in enumerate(axis_ymd)}
print("authoritative day axis:", n_days, axis_ymd[0], "->", axis_ymd[-1])

# gather stations across the 12 yearly CSVs
files = sorted(glob.glob(os.path.join(STATION, "staD10_yearly", "staD10_2*.csv")))
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
for yi, f in enumerate(files):
    d = pd.read_csv(f, na_values=["--", ""])
    dcols = [c for c in d.columns if c.startswith("D")]
    ids = d["Station_ID_d"].astype(int).to_numpy()
    X = d[dcols].to_numpy(float)
    col_idx = np.array([day_pos[c[1:]] if c[1:] in day_pos else -1 for c in dcols])
    keep_cols = col_idx >= 0
    X = X[:, keep_cols]
    col_idx = col_idx[keep_cols]
    rows = np.array([sta_pos[s] for s in ids])
    present_year[rows, yi] = True
    # fill raw (only finite values)
    for ri, rr in enumerate(rows):
        vals = X[ri]
        fin = np.isfinite(vals)
        raw[rr, col_idx[fin]] = vals[fin]

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
