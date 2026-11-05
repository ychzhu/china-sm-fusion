# -*- coding: utf-8 -*-
"""Full-domain CDF bias adjustment of the SMAP background."""
from paths_full import *

import os
import numpy as np

THRESH_NUM, MIN_HIST = 3, 30

g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
n = np.load(os.path.join(WORK, "neighbors_full.npz"), allow_pickle=True)
smap_bg = g["smap_bg"]
grid_obs = g["grid_obs"]
num_around = n["num_around"]
cap_idx = n["cap_idx"]
n_days, n_grid = smap_bg.shape


def interp1_extrap(x, y, xi):
    out = np.interp(xi, x, y)
    left = xi < x[0]
    if left.any():
        slope = (y[1] - y[0]) / (x[1] - x[0])
        out[left] = y[0] + slope * (xi[left] - x[0])
    right = xi > x[-1]
    if right.any():
        slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
        out[right] = y[-1] + slope * (xi[right] - x[-1])
    return out


def cdf_adjusting2(smap, sta):
    s = np.asarray(smap, float).copy()
    t = np.asarray(sta, float).copy()
    s[(s > 0.6) | (s < 0.04)] = np.nan
    t[(t > 0.6) | (t < 0.04)] = np.nan
    sv, tv = s[~np.isnan(s)], t[~np.isnan(t)]
    if len(sv) < 2 or len(tv) < 2:
        return s
    ps = np.arange(1, len(sv) + 1) / (len(sv) + 1)
    pt = np.arange(1, len(tv) + 1) / (len(tv) + 1)
    sint = interp1_extrap(ps, np.sort(sv), pt)
    diff = np.sort(tv) - sint
    coeff = np.polyfit(sint, diff, 2)
    out = s + np.polyval(coeff, s)
    out[(out > 1) | (out < 0)] = np.nan
    return out


smap_cdf = smap_bg.copy()
applied = 0
for gi in range(n_grid):
    valid = np.isfinite(smap_bg[:, gi])
    if num_around[gi] >= THRESH_NUM and valid.sum() > MIN_HIST:
        ci = cap_idx[gi]
        ci = ci[ci >= 0]
        sta_ref = grid_obs[ci, :].flatten()
        smap_cdf[:, gi] = cdf_adjusting2(smap_bg[:, gi], sta_ref)
        applied += 1
print("CDF applied to %d / %d grids" % (applied, n_grid))
print("smap background finite=%d CDF finite=%d" %
      (np.isfinite(smap_bg).sum(), np.isfinite(smap_cdf).sum()))
np.save(os.path.join(WORK, "smap_cdf_full.npy"), smap_cdf)
print("saved smap_cdf_full.npy")
