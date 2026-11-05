# -*- coding: utf-8 -*-
"""Full-domain leak-free 10-fold CV (SMAP error constant 0.046)."""
from paths_full import *

import os
import numpy as np

A, THRESH, SEED, NFOLDS = 10.0, 3, 2026, 10

g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
nb = np.load(os.path.join(WORK, "neighbors_full.npz"), allow_pickle=True)
bearing, glat, glon = g["bearing"], g["glat"], g["glon"]
grid_err, grid_obs = g["grid_err"], g["grid_obs"]
SMAP_ERR = float(g["smap_err"][0])
smap_cdf = np.load(os.path.join(WORK, "smap_cdf_full.npy"))
cap_idx, cap_dist = nb["cap_idx"], nb["cap_dist"]
nEff = len(bearing)
n_days = grid_obs.shape[1]
blat, blon = glat[bearing], glon[bearing]
err_O = grid_err[bearing]

rng = np.random.RandomState(SEED)
perm = rng.permutation(nEff)
edges = np.round(np.linspace(0, nEff, NFOLDS + 1)).astype(int)
fold_of = np.full(nEff, -1, int)
foldV = []
for k in range(NFOLDS):
    m = perm[edges[k]:edges[k + 1]]
    foldV.append(m)
    fold_of[m] = k
print("fold sizes:", [len(f) for f in foldV])

full_nbr, full_B, full_corIk = [], [], []
for v in range(nEff):
    gv = bearing[v]
    ci = cap_idx[gv]
    ci = ci[ci >= 0]
    cd = cap_dist[gv][:len(ci)]
    la, lo = blat[ci], blon[ci]
    nA = len(ci)
    dmat = np.zeros((nA, nA))
    for r in range(nA):
        p1, p2 = np.radians(la[r]), np.radians(la)
        dphi = np.radians(la - la[r])
        dl = np.radians(lo - lo[r])
        aa = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
        dmat[r] = np.degrees(2 * np.arcsin(np.sqrt(aa)))
    full_nbr.append(ci)
    full_B.append(1.0 / np.power(A, dmat))
    full_corIk.append(1.0 / np.power(A, cd))

allB, allF, allO = [], [], []
for k in range(NFOLDS):
    V = foldV[k]
    Bv = np.full((len(V), n_days), np.nan)
    Fv = np.full((len(V), n_days), np.nan)
    Ov = np.full((len(V), n_days), np.nan)
    for ii, v in enumerate(V):
        gv = bearing[v]
        nbr = full_nbr[v]
        train_sel = fold_of[nbr] != k
        nbrT = nbr[train_sel]
        Bv[ii] = smap_cdf[:, gv]
        Ov[ii] = grid_obs[v]
        if len(nbrT) < THRESH:
            Fv[ii] = Bv[ii]
            continue
        Bsub = full_B[v][np.ix_(train_sel, train_sel)]
        cIk = full_corIk[v][train_sel]
        errCon = err_O[nbrT] / SMAP_ERR
        bad = ~np.isfinite(errCon) | (errCon < 0)
        errCon[bad] = 1.0
        w = np.linalg.pinv(Bsub + np.diag(errCon ** 2)) @ cIk
        for di in range(n_days):
            resid = grid_obs[nbrT, di] - smap_cdf[di, bearing[nbrT]]
            resid = np.where(np.isfinite(resid), resid, 0.0)
            Fv[ii, di] = smap_cdf[di, gv] + np.sum(resid * w)
    ok = np.isfinite(Ov) & np.isfinite(Bv) & np.isfinite(Fv)
    allB.append(Bv[ok]); allF.append(Fv[ok]); allO.append(Ov[ok])
    print("fold %d: pairs=%d MAE bg %.4f -> OI %.4f" %
          (k + 1, ok.sum(), np.mean(np.abs(Bv[ok] - Ov[ok])),
           np.mean(np.abs(Fv[ok] - Ov[ok]))))

allB = np.concatenate(allB); allF = np.concatenate(allF); allO = np.concatenate(allO)
print("===== 10-fold CV overall =====")
print("total pairs=%d" % len(allO))
print("MAE : bg %.4f -> OI %.4f (improve %.1f%%)" %
      (np.mean(np.abs(allB - allO)), np.mean(np.abs(allF - allO)),
       100 * (1 - np.mean(np.abs(allF - allO)) / np.mean(np.abs(allB - allO)))))
print("RMSE: bg %.4f -> OI %.4f (improve %.1f%%)" %
      (np.sqrt(np.mean((allB - allO) ** 2)), np.sqrt(np.mean((allF - allO) ** 2)),
       100 * (1 - np.sqrt(np.mean((allF - allO) ** 2)) / np.sqrt(np.mean((allB - allO) ** 2)))))
print("CORR: bg %.4f -> OI %.4f" %
      (np.corrcoef(allB, allO)[0, 1], np.corrcoef(allF, allO)[0, 1]))
print("BIAS: bg %+.4f -> OI %+.4f" %
      (np.mean(allB - allO), np.mean(allF - allO)))
print("pairs improved: %.1f%%" %
      (100 * np.mean(np.abs(allF - allO) < np.abs(allB - allO))))
np.savez(os.path.join(WORK, "cv_results_full.npz"), allB=allB, allF=allF, allO=allO)
