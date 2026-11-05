# -*- coding: utf-8 -*-
"""Full-domain station-level comparison against independent references."""
from paths_full import *

import os
import numpy as np
import netCDF4
from scipy.spatial import cKDTree


g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
glat, glon, common_ymd = g["glat"], g["glon"], g["common_ymd"]
smap_bg = g["smap_bg"]
smap_cdf = np.load(os.path.join(WORK, "smap_cdf_full.npy"))
fusion = np.load(os.path.join(WORK, "fusion_full.npy"))

slat = np.load(os.path.join(WORK, "lat3596.npy"))
slon = np.load(os.path.join(WORK, "lon3596.npy"))
cleaned = np.load(os.path.join(WORK, "obj_clean3596.npy"))
sta_ymd = station_axis_ymd()
sta_pos = {d: i for i, d in enumerate(sta_ymd)}
sta_obs = cleaned[:, np.array([sta_pos[d] for d in common_ymd])]

_, ng = cKDTree(np.column_stack([glat, glon])).query(
    np.column_stack([slat, slon]), k=1)


def open_nc(path):
    with open(path, "rb") as fh:
        return netCDF4.Dataset(os.path.basename(path), memory=fh.read())


def load_ref(rel):
    rds = open_nc(rel)
    rymd = rds.variables["time_ymd"][:].astype(str)
    rsm = np.ma.filled(rds.variables["sm"][:], np.nan).astype(float)
    rds.close()
    rmap = {d: i for i, d in enumerate(rymd)}
    return rsm[np.array([rmap[d] for d in common_ymd])]


era5 = load_ref(os.path.join(REFERENCE, "era5_gcs_swvl1_36km_20150331_20261001.nc"))
gldas = load_ref(os.path.join(REFERENCE, "gldas_noah_soilm_36km_20150331_20261001.nc"))

systems = {
    "GLDAS": gldas[:, ng], "ERA5": era5[:, ng],
    "SMAP_raw": smap_bg[:, ng], "SMAP_CDF": smap_cdf[:, ng],
    "Fusion": fusion[:, ng]}

print("system vs station observations (pooled station-day pairs)")
for name, arr in systems.items():
    a = arr.T
    valid = np.isfinite(a) & np.isfinite(sta_obs)
    x, y = a[valid], sta_obs[valid]
    d = x - y
    print("%-9s N=%9d MAE=%.4f RMSE=%.4f BIAS=%+.4f CORR=%.4f" %
          (name, len(x), np.mean(np.abs(d)), np.sqrt(np.mean(d ** 2)),
           np.mean(d), np.corrcoef(x, y)[0, 1]))
