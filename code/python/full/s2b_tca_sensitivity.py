# -*- coding: utf-8 -*-
"""TCA sensitivity to SMAP two-day backfill and production-grid support."""
from paths_full import *
from validation_utils import map_stations

import csv
import json
import os
import numpy as np
import pandas as pd

MIN_COMMON = 100
clean = np.load(os.path.join(WORK, "obj_clean3596.npy"))
slat = np.load(os.path.join(WORK, "lat3596.npy"))
slon = np.load(os.path.join(WORK, "lon3596.npy"))
sids = np.load(os.path.join(WORK, "ids3596.npy"))
s_ymd = station_axis_ymd()
glat, glon = grid_geometry()
g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
s2g, assigned, _ = map_stations(g["members"], glat, glon, slat, slon)

sm = open_nc(os.path.join(SMAP, "smapOI_daily_20150331_20261017.nc"))
m_ymd = sm.variables["time_ymd"][:].astype(str)
m_sm = np.ma.filled(sm.variables["sm"][:], np.nan)
sm.close()
er = open_nc(os.path.join(REFERENCE, "era5_gcs_swvl1_36km_20150331_20261001.nc"))
e_ymd = er.variables["time_ymd"][:].astype(str)
e_sm = np.ma.filled(er.variables["sm"][:], np.nan)
er.close()
si = {d: k for k, d in enumerate(s_ymd)}
mi = {d: k for k, d in enumerate(m_ymd)}
ei = {d: k for k, d in enumerate(e_ymd)}
common = sorted(set(s_ymd.tolist()) & set(mi) & set(ei))
s_idx = np.array([si[d] for d in common])
m_idx = np.array([mi[d] for d in common])
e_idx = np.array([ei[d] for d in common])
rows = []
for k in range(len(sids)):
    gi = int(s2g[k])
    x1 = clean[k, s_idx]
    x2 = m_sm[m_idx, gi]
    x3 = e_sm[e_idx, gi]
    ok = np.isfinite(x1) & np.isfinite(x2) & np.isfinite(x3)
    if ok.sum() < MIN_COMMON:
        continue
    cov = np.cov(np.column_stack([x1[ok], x2[ok], x3[ok]]), rowvar=False)
    q1 = cov[0, 0] - cov[0, 1] * cov[0, 2] / cov[1, 2]
    q2 = cov[1, 1] - cov[0, 1] * cov[1, 2] / cov[0, 2]
    q3 = cov[2, 2] - cov[0, 2] * cov[1, 2] / cov[0, 1]
    if q1 <= 0 or q2 <= 0 or q3 <= 0:
        continue
    rows.append([k, gi, np.sqrt(q1), np.sqrt(q2), np.sqrt(q3), int(ok.sum()), int(assigned[k])])
with open(os.path.join(WORK, "tca3596_daily_only.csv"), "w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh)
    writer.writerow(["station_index", "grid_index", "station_sigma", "smap_sigma",
                     "era5_sigma", "n_common", "inside_production_grid"])
    writer.writerows(rows)

daily = pd.DataFrame(rows, columns=["station_index", "grid_index", "station_sigma_daily",
                                    "smap_sigma_daily", "era5_sigma_daily", "n_common_daily",
                                    "inside_production_grid"])
base = pd.read_csv(os.path.join(WORK, "tca3596.csv"))
base["inside_production_grid"] = assigned[base.station_index.to_numpy()].astype(int)
merged = base.merge(daily, on="station_index", how="inner")
merged.to_csv(os.path.join(WORK, "tca_sensitivity_by_station.csv"), index=False)
summary = {
    "baseline_all": {
        "N": int(len(base)),
        "station_sigma_median": float(base.station_sigma.median()),
        "smap_sigma_median": float(base.smap_sigma.median()),
        "era5_sigma_median": float(base.era5_sigma.median()),
    },
    "baseline_inside_only": {
        "N": int(base.inside_production_grid.sum()),
        "station_sigma_median": float(base.loc[base.inside_production_grid == 1, "station_sigma"].median()),
        "smap_sigma_median": float(base.loc[base.inside_production_grid == 1, "smap_sigma"].median()),
        "era5_sigma_median": float(base.loc[base.inside_production_grid == 1, "era5_sigma"].median()),
    },
    "daily_only_all": {
        "N": int(len(daily)),
        "station_sigma_median": float(daily.station_sigma_daily.median()),
        "smap_sigma_median": float(daily.smap_sigma_daily.median()),
        "era5_sigma_median": float(daily.era5_sigma_daily.median()),
    },
    "paired_daily_minus_2day": {
        "N": int(len(merged)),
        "station_sigma_median_delta": float((merged.station_sigma_daily - merged.station_sigma).median()),
        "smap_sigma_median_delta": float((merged.smap_sigma_daily - merged.smap_sigma).median()),
        "era5_sigma_median_delta": float((merged.era5_sigma_daily - merged.era5_sigma).median()),
    },
}
with open(os.path.join(WORK, "tca_sensitivity.json"), "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)
print(json.dumps(summary, indent=2))
