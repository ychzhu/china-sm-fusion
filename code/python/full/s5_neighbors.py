# -*- coding: utf-8 -*-
"""Full-domain capped neighbour selection (R1=1.08, MAX_NEAR=9)."""
from paths_full import *

import os
import numpy as np

R1 = 1.08
MAX_NEAR = 9

d = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
bearing = d["bearing"]
glat, glon = d["glat"], d["glon"]
blat, blon = glat[bearing], glon[bearing]
n_grid = len(glat)


def haversine_deg(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dl = np.radians(lon2 - lon1)
    aa = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(aa)))


num_around = np.zeros(n_grid, int)
cap_idx = np.full((n_grid, MAX_NEAR), -1, int)
cap_dist = np.full((n_grid, MAX_NEAR), np.nan)
for gi in range(n_grid):
    cand = np.where((np.abs(blat - glat[gi]) <= R1) &
                    (np.abs(blon - glon[gi]) <= R1))[0]
    num_around[gi] = len(cand)
    if len(cand) == 0:
        continue
    dist = haversine_deg(glat[gi], glon[gi], blat[cand], blon[cand])
    # deterministic tie-break: distance first, then bearing position
    keep = np.lexsort((cand, dist))[:MAX_NEAR]
    cap_idx[gi, :len(keep)] = cand[keep]
    cap_dist[gi, :len(keep)] = dist[keep]

print("grids with >=3 neighbors:", (num_around >= 3).sum())
print("num_around min=%d median=%d max=%d" %
      (num_around.min(), np.median(num_around), num_around.max()))
np.savez(os.path.join(WORK, "neighbors_full.npz"),
         num_around=num_around, cap_idx=cap_idx, cap_dist=cap_dist)
print("saved neighbors_full.npz")
