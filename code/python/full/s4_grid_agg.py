# -*- coding: utf-8 -*-
"""Station-to-grid aggregation with day-specific observation error (rho=0)."""
from paths_full import *

import os
import numpy as np
import pandas as pd
import netCDF4

GRID_RES = 0.18
FALLBACK = 0.047
SMAP_ERR_NEW = 0.046


def open_nc(path):
    with open(path, "rb") as fh:
        return netCDF4.Dataset(os.path.basename(path), memory=fh.read())


glat, glon = grid_geometry()
n_grid = len(glat)
cleaned = np.load(os.path.join(WORK, "obj_clean3596.npy"))
slat = np.load(os.path.join(WORK, "lat3596.npy"))
slon = np.load(os.path.join(WORK, "lon3596.npy"))
sids = np.load(os.path.join(WORK, "ids3596.npy"))
sta_ymd = station_axis_ymd()

sm = open_nc(os.path.join(SMAP, "smapOI_2day_20150331_20261017.nc"))
smap_ymd = sm.variables["time_ymd"][:].astype(str)
smap_bg = np.ma.filled(sm.variables["sm"][:], np.nan)
sm.close()

members = [np.where((np.abs(slat - glat[gi]) <= GRID_RES) &
                    (np.abs(slon - glon[gi]) <= GRID_RES))[0]
           for gi in range(n_grid)]
n_with = np.array([len(m) for m in members])
bearing = np.where(n_with > 0)[0]
print("bearing grids=%d; unique assigned stations=%d" %
      (len(bearing), len(np.unique(np.concatenate([members[i] for i in bearing])))))

tca = pd.read_csv(os.path.join(WORK, "tca3596.csv"))
sta_sigma = np.full(len(sids), np.nan)
sta_sigma[tca.station_index.to_numpy()] = tca.station_sigma.to_numpy()

# Retain the old static estimate only for audit/backward comparison.
grid_err = np.full(n_grid, np.nan)
for gi in bearing:
    s = sta_sigma[members[gi]].copy()
    s[~np.isfinite(s) | (s <= 0)] = FALLBACK
    grid_err[gi] = np.sqrt(np.sum(s ** 2)) / len(members[gi])

common = np.array(sorted(set(sta_ymd.tolist()) & set(smap_ymd.tolist())))
sp = {d: i for i, d in enumerate(sta_ymd)}
mp = {d: i for i, d in enumerate(smap_ymd)}
sta_aligned = cleaned[:, np.array([sp[d] for d in common])]
smap_bg = smap_bg[np.array([mp[d] for d in common])]

grid_obs = np.full((len(bearing), len(common)), np.nan, dtype=np.float64)
grid_err_daily = np.full((len(bearing), len(common)), np.nan, dtype=np.float32)
grid_nvalid = np.zeros((len(bearing), len(common)), dtype=np.int16)
for k, gi in enumerate(bearing):
    values = sta_aligned[members[gi]]
    sigma = sta_sigma[members[gi]].copy()
    sigma[~np.isfinite(sigma) | (sigma <= 0)] = FALLBACK
    valid = np.isfinite(values)
    n = valid.sum(axis=0)
    grid_nvalid[k] = n
    has = n > 0
    grid_obs[k, has] = np.sum(np.where(valid, values, 0.0), axis=0)[has] / n[has]
    # sigma_g(t) = sqrt(sum_{i in g(t)} sigma_i^2) / n_g(t).
    sum_var = np.sum(np.where(valid, sigma[:, None] ** 2, 0.0), axis=0)
    grid_err_daily[k, has] = np.sqrt(sum_var[has]) / n[has]

print("common analysis days=%d grid_obs finite=%d" %
      (len(common), np.isfinite(grid_obs).sum()))
print("daily active stations: median=%.2f mean=%.2f max=%d" %
      (np.median(grid_nvalid[grid_nvalid > 0]),
       np.mean(grid_nvalid[grid_nvalid > 0]), grid_nvalid.max()))
print("daily cell error median=%.6f mean=%.6f" %
      (np.nanmedian(grid_err_daily), np.nanmean(grid_err_daily)))

np.savez(os.path.join(WORK, "grid_prep_full.npz"),
         members=np.array(members, dtype=object), bearing=bearing,
         n_with=n_with, grid_err=grid_err, grid_obs=grid_obs,
         grid_err_daily=grid_err_daily, grid_nvalid=grid_nvalid,
         glat=glat, glon=glon, common_ymd=common,
         smap_bg=smap_bg, smap_err=np.full(len(common), SMAP_ERR_NEW))
print("saved grid_prep_full.npz with dynamic daily R")
