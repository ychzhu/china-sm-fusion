# -*- coding: utf-8 -*-
"""Full-domain (3596) station -> 36-km grid aggregation + cell error (rho=0)."""
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


# grid geometry (3596 fusion product)
glat, glon = grid_geometry()
n_grid = len(glat)

# final cleaned station matrix + geometry
cleaned = np.load(os.path.join(WORK, "obj_clean3596.npy"))
slat = np.load(os.path.join(WORK, "lat3596.npy"))
slon = np.load(os.path.join(WORK, "lon3596.npy"))
sids = np.load(os.path.join(WORK, "ids3596.npy"))
sta_ymd = station_axis_ymd()

# SMAP 2-day background
sm = open_nc(os.path.join(SMAP, "smapOI_2day_20150331_20261017.nc"))
smap_ymd = sm.variables["time_ymd"][:].astype(str)
smap_bg = np.ma.filled(sm.variables["sm"][:], np.nan)
sm.close()

# box membership
members = [np.where((np.abs(slat - glat[gi]) <= GRID_RES) &
                    (np.abs(slon - glon[gi]) <= GRID_RES))[0] for gi in range(n_grid)]
n_with = np.array([len(m) for m in members])
bearing = np.where(n_with > 0)[0]
print("bearing grids=%d (old 1714), assigned stations=%d" %
      (len(bearing), (n_with > 0).sum()))

# new TCA per-station sigma
tca = pd.read_csv(os.path.join(WORK, "tca3596.csv"))
sta_sigma = np.full(len(sids), np.nan)
sta_sigma[tca.station_index.to_numpy()] = tca.station_sigma.to_numpy()

# static cell error rho=0
grid_err = np.full(n_grid, np.nan)
for gi in bearing:
    s = sta_sigma[members[gi]].copy()
    s[~np.isfinite(s) | (s <= 0)] = FALLBACK
    grid_err[gi] = np.sqrt(np.sum(s ** 2)) / len(members[gi])
print("cell error median=%.4f mean=%.4f (old median 0.0386)" %
      (np.median(grid_err[bearing]), np.mean(grid_err[bearing])))

# common axis
common = np.array(sorted(set(sta_ymd.tolist()) & set(smap_ymd.tolist())))
sp = {d: i for i, d in enumerate(sta_ymd)}
mp = {d: i for i, d in enumerate(smap_ymd)}
sel_s = np.array([sp[d] for d in common])
sel_m = np.array([mp[d] for d in common])
sta_aligned = cleaned[:, sel_s]
smap_bg = smap_bg[sel_m]

grid_obs = np.full((len(bearing), len(common)), np.nan)
for k, gi in enumerate(bearing):
    with np.errstate(all="ignore"):
        grid_obs[k] = np.nanmean(sta_aligned[members[gi]], axis=0)
print("common analysis days=%d grid_obs finite=%d" %
      (len(common), np.isfinite(grid_obs).sum()))

# stations per bearing cell distribution
print("stations per bearing cell: median=%.0f mean=%.1f max=%d" %
      (np.median(n_with[bearing]), n_with[bearing].mean(), n_with[bearing].max()))

np.savez(os.path.join(WORK, "grid_prep_full.npz"),
         members=np.array(members, dtype=object), bearing=bearing,
         n_with=n_with, grid_err=grid_err, grid_obs=grid_obs,
         glat=glat, glon=glon, common_ymd=common,
         smap_bg=smap_bg, smap_err=np.full(len(common), SMAP_ERR_NEW))
print("saved grid_prep_full.npz")
