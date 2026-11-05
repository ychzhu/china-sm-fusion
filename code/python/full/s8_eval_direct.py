# -*- coding: utf-8 -*-
"""Full-domain direct (in-sample) validation."""
from paths_full import *

import os
import numpy as np


g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
bearing = g["bearing"]
grid_obs = g["grid_obs"]
smap_cdf = np.load(os.path.join(WORK, "smap_cdf_full.npy"))
fusion = np.load(os.path.join(WORK, "fusion_full.npy"))
n_days = grid_obs.shape[1]

bg_b = smap_cdf[:, bearing]
fu_b = fusion[:, bearing]


def pooled_metrics(pred, obs):
    diff = pred - obs
    pv = pred - np.mean(pred)
    ov = obs - np.mean(obs)
    corr = np.sum(pv * ov) / np.sqrt(np.sum(pv ** 2) * np.sum(ov ** 2))
    return np.mean(np.abs(diff)), np.sqrt(np.mean(diff ** 2)), np.mean(diff), corr


out = {}
for name, pred in (("background", bg_b.T), ("fusion", fu_b.T)):
    valid = np.isfinite(grid_obs) & np.isfinite(pred)
    out[name] = (valid.sum(), pooled_metrics(pred[valid], grid_obs[valid]))
    dmae, dcorr = [], []
    for di in range(n_days):
        v = np.isfinite(grid_obs[:, di]) & np.isfinite(pred[:, di])
        if v.sum() < 3:
            continue
        dmae.append(np.mean(np.abs(pred[v, di] - grid_obs[v, di])))
        if np.std(pred[v, di]) > 0 and np.std(grid_obs[v, di]) > 0:
            dcorr.append(np.corrcoef(pred[v, di], grid_obs[v, di])[0, 1])
    out[name + "_daily"] = (np.mean(dmae), np.mean(dcorr))

for name in ("background", "fusion"):
    n, m = out[name]
    print("=== %s === pairs=%d" % (name, n))
    print("  pooled MAE=%.4f RMSE=%.4f BIAS=%+.4f CORR=%.4f" % m)
    print("  daily  MAE=%.4f CORR=%.4f" % out[name + "_daily"])
b, f = out["background"][1], out["fusion"][1]
print("in-sample improvement: MAE %.1f%%  RMSE %.1f%%" %
      (100 * (1 - f[0] / b[0]), 100 * (1 - f[1] / b[1])))
