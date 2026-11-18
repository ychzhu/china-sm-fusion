# -*- coding: utf-8 -*-
"""Dynamic-R OI fusion with standard posterior analysis uncertainty."""
from paths_full import *
from oi_utils import build_geometry, dynamic_oi_day

import json
import os
import numpy as np

A = 10.0
THRESH_NUM, MIN_HIST = 3, 30

g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
nb = np.load(os.path.join(WORK, "neighbors_full.npz"), allow_pickle=True)
bearing = g["bearing"]
glat, glon = g["glat"], g["glon"]
grid_obs = g["grid_obs"]
grid_err_daily = g["grid_err_daily"].astype(np.float64)
SMAP_ERR = float(g["smap_err"][0])
smap_cdf = np.load(os.path.join(WORK, "smap_cdf_full.npy"))
num_around, cap_idx, cap_dist = nb["num_around"], nb["cap_idx"], nb["cap_dist"]
n_days, n_grid = smap_cdf.shape

active_grids = []
neighbor_lists = []
distance_lists = []
for gi in range(n_grid):
    if num_around[gi] < THRESH_NUM or np.isfinite(smap_cdf[:, gi]).sum() <= MIN_HIST:
        continue
    bp = cap_idx[gi]
    bp = bp[bp >= 0]
    active_grids.append(gi)
    neighbor_lists.append(bp)
    distance_lists.append(cap_dist[gi][:len(bp)])
active_grids = np.asarray(active_grids, dtype=np.int32)
bp_pad, b_mat, c_vec = build_geometry(
    glat, glon, bearing, active_grids, neighbor_lists, distance_lists, A)
bp_safe = np.where(bp_pad >= 0, bp_pad, 0)
print("dynamic OI grids=%d" % len(active_grids))

fusion = smap_cdf.copy()
fusion_err = np.where(np.isfinite(smap_cdf), SMAP_ERR, np.nan)
fusion_applied = np.zeros((n_days, n_grid), dtype=bool)
clip_count = 0
for di in range(n_days):
    target_bg = smap_cdf[di, active_grids]
    obs = grid_obs[bp_safe, di]
    neighbor_bg = smap_cdf[di, bearing[bp_safe]]
    neighbor_err = grid_err_daily[bp_safe, di]
    inc, post, has_info, reduction = dynamic_oi_day(
        target_bg, obs, neighbor_bg, neighbor_err, bp_pad,
        b_mat, c_vec, SMAP_ERR)
    fusion[di, active_grids] = target_bg + inc
    fusion_err[di, active_grids] = post
    fusion_applied[di, active_grids] = has_info
    clip_count += int(np.sum((reduction < -1e-10) | (reduction > 1 + 1e-10)))
    if di % 250 == 0 or di == n_days - 1:
        print("day %d/%d" % (di + 1, n_days), flush=True)

valid = np.isfinite(fusion_err)
active_err = fusion_err[:, active_grids]
active_valid = np.isfinite(active_err)
diag = {
    "model": "dynamic daily observation R; standard OI posterior sigma_b*sqrt(1-c^T w)",
    "a": A,
    "sigma_b": SMAP_ERR,
    "n_oi_grids": int(len(active_grids)),
    "actual_oi_grid_days": int(fusion_applied.sum()),
    "n_days": int(n_days),
    "aerr_median_all": float(np.median(fusion_err[valid])),
    "aerr_mean_all": float(np.mean(fusion_err[valid])),
    "aerr_p10_all": float(np.percentile(fusion_err[valid], 10)),
    "aerr_p90_all": float(np.percentile(fusion_err[valid], 90)),
    "aerr_median_oi": float(np.median(active_err[active_valid])),
    "aerr_mean_oi": float(np.mean(active_err[active_valid])),
    "aerr_above_background": int(np.sum(fusion_err[valid] > SMAP_ERR + 1e-12)),
    "variance_reduction_outside_tolerance": clip_count,
}
with open(os.path.join(WORK, "fusion_diagnostics.json"), "w", encoding="utf-8") as fh:
    json.dump(diag, fh, ensure_ascii=False, indent=2)
print(json.dumps(diag, ensure_ascii=False, indent=2))
np.save(os.path.join(WORK, "fusion_full.npy"), fusion)
np.save(os.path.join(WORK, "fusion_full_err.npy"), fusion_err)
np.save(os.path.join(WORK, "fusion_do_oi.npy"), fusion_applied)
print("saved dynamic fusion, posterior Aerr and actual daily OI flag")
