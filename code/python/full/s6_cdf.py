# -*- coding: utf-8 -*-
"""Full-domain CDF bias adjustment of the SMAP background."""
from paths_full import *
from cdf_utils import cdf_adjusting2

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
