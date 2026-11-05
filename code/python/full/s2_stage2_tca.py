# -*- coding: utf-8 -*-
"""Stage 2: full-domain Triple Collocation on objective-cleaned 3596 stations."""
from paths_full import *

import os
import csv
import numpy as np
import netCDF4
from scipy.spatial import cKDTree

MIN_COMMON = 100


def open_nc(path):
    with open(path, "rb") as fh:
        return netCDF4.Dataset(os.path.basename(path), memory=fh.read())


# objective-cleaned station matrix + geometry
clean = np.load(os.path.join(WORK, "obj_clean3596.npy"))
slat = np.load(os.path.join(WORK, "lat3596.npy"))
slon = np.load(os.path.join(WORK, "lon3596.npy"))
sids = np.load(os.path.join(WORK, "ids3596.npy"))
s_ymd = station_axis_ymd()

# grid geometry (3596 fusion product, same 7313-cell grid)
glat, glon = grid_geometry()

# SMAP 2-day
sm = open_nc(os.path.join(SMAP, "smapOI_2day_20150331_20261017.nc"))
m_ymd = sm.variables["time_ymd"][:].astype(str)
m_sm = np.ma.filled(sm.variables["sm"][:], np.nan)
sm.close()

# ERA5
er = open_nc(os.path.join(REFERENCE, "era5_gcs_swvl1_36km_20150331_20261001.nc"))
e_ymd = er.variables["time_ymd"][:].astype(str)
e_sm = np.ma.filled(er.variables["sm"][:], np.nan)
er.close()

# nearest grid for each station
tree = cKDTree(np.column_stack([glat, glon]))
_, s2g = tree.query(np.column_stack([slat, slon]), k=1)

si = {d: k for k, d in enumerate(s_ymd)}
mi = {d: k for k, d in enumerate(m_ymd)}
ei = {d: k for k, d in enumerate(e_ymd)}
common = sorted(set(s_ymd.tolist()) & set(mi) & set(ei))
s_idx = np.array([si[d] for d in common])
m_idx = np.array([mi[d] for d in common])
e_idx = np.array([ei[d] for d in common])
print("common days:", len(common), common[0], "->", common[-1])

rows = []
fail_n = 0
for k in range(len(sids)):
    g = int(s2g[k])
    x1 = clean[k, s_idx]
    x2 = m_sm[m_idx, g]
    x3 = e_sm[e_idx, g]
    ok = np.isfinite(x1) & np.isfinite(x2) & np.isfinite(x3)
    if ok.sum() < MIN_COMMON:
        fail_n += 1
        continue
    c = np.cov(np.column_stack([x1[ok], x2[ok], x3[ok]]), rowvar=False)
    q1 = c[0, 0] - c[0, 1] * c[0, 2] / c[1, 2]
    q2 = c[1, 1] - c[0, 1] * c[1, 2] / c[0, 2]
    q3 = c[2, 2] - c[0, 2] * c[1, 2] / c[0, 1]
    if q1 <= 0 or q2 <= 0 or q3 <= 0:
        fail_n += 1
        continue
    rows.append([k, g, np.sqrt(q1), np.sqrt(q2), np.sqrt(q3), int(ok.sum())])

with open(os.path.join(WORK, "tca3596.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["station_index", "grid_index", "station_sigma",
                "smap_sigma", "era5_sigma", "n_common"])
    w.writerows(rows)

arr = np.array(rows)
print("stations with TCA: %d (failed/insufficient: %d)" % (len(rows), fail_n))
for j, name in [(2, "station"), (3, "smap"), (4, "era5")]:
    v = arr[:, j]
    print("  %-7s median=%.4f mean=%.4f p10=%.4f p90=%.4f min=%.4f max=%.4f" %
          (name, np.median(v), v.mean(), np.percentile(v, 10),
           np.percentile(v, 90), v.min(), v.max()))
print("common days median=%.0f min=%.0f max=%.0f" %
      (np.median(arr[:, 5]), arr[:, 5].min(), arr[:, 5].max()))
print("saved tca3596.csv")
