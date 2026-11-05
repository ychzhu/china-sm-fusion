# -*- coding: utf-8 -*-
"""Full-domain OI fusion (SMAP background error constant 0.046)."""
from paths_full import *

import os
import numpy as np

A = 10.0
THRESH_NUM, MIN_HIST = 3, 30

g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
nb = np.load(os.path.join(WORK, "neighbors_full.npz"), allow_pickle=True)
bearing = g["bearing"]
glat, glon = g["glat"], g["glon"]
grid_err = g["grid_err"]
grid_obs = g["grid_obs"]
SMAP_ERR = float(g["smap_err"][0])
smap_cdf = np.load(os.path.join(WORK, "smap_cdf_full.npy"))
num_around, cap_idx, cap_dist = nb["num_around"], nb["cap_idx"], nb["cap_dist"]
n_days, n_grid = smap_cdf.shape

err_O_bearing = grid_err[bearing]
blat, blon = glat[bearing], glon[bearing]


def haversine(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dl = np.radians(lon2 - lon1)
    aa = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(aa)))


doOI = np.zeros(n_grid, bool)
geo = [None] * n_grid
for gi in range(n_grid):
    days_valid = np.isfinite(smap_cdf[:, gi]).sum()
    if num_around[gi] < THRESH_NUM or days_valid <= MIN_HIST:
        continue
    bp = cap_idx[gi]
    bp = bp[bp >= 0]
    nA = len(bp)
    la, lo = blat[bp], blon[bp]
    dmat = np.zeros((nA, nA))
    for r in range(nA):
        dmat[r] = haversine(la[r], lo[r], la, lo)
    B = 1.0 / np.power(A, dmat)
    corIk = 1.0 / np.power(A, cap_dist[gi][:nA])
    errCon = err_O_bearing[bp] / SMAP_ERR
    bad = ~np.isfinite(errCon) | (errCon < 0)
    errCon[bad] = 1.0
    static_w = np.linalg.pinv(B + np.diag(errCon ** 2)) @ corIk
    doOI[gi] = True
    geo[gi] = {"bp": bp, "static_w": static_w}
print("doOI grids=%d" % doOI.sum())

fusion = smap_cdf.copy()
fusion_err = np.full((n_days, n_grid), SMAP_ERR)
for di in range(n_days):
    for gi in range(n_grid):
        if not doOI[gi]:
            continue
        bp = geo[gi]["bp"]
        w = geo[gi]["static_w"]
        resid = grid_obs[bp, di] - smap_cdf[di, bearing[bp]]
        resid = np.where(np.isfinite(resid), resid, 0.0)
        fusion[di, gi] = smap_cdf[di, gi] + np.sum(resid * w)
        e_diff = err_O_bearing[bp] - SMAP_ERR
        fusion_err[di, gi] = SMAP_ERR + np.sum(e_diff * w)

print("fusion finite=%d (smap_cdf finite=%d)" %
      (np.isfinite(fusion).sum(), np.isfinite(smap_cdf).sum()))
np.save(os.path.join(WORK, "fusion_full.npy"), fusion)
np.save(os.path.join(WORK, "fusion_full_err.npy"), fusion_err)
print("saved fusion_full.npy and fusion_full_err.npy")
