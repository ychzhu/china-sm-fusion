# -*- coding: utf-8 -*-
"""Shared validation metrics and station-to-production-grid mapping."""
import numpy as np
from scipy.spatial import cKDTree


def metrics(pred, obs):
    pred = np.asarray(pred, dtype=np.float64)
    obs = np.asarray(obs, dtype=np.float64)
    d = pred - obs
    return {
        "N": int(len(obs)),
        "MAE": float(np.mean(np.abs(d))),
        "RMSE": float(np.sqrt(np.mean(d ** 2))),
        "BIAS": float(np.mean(d)),
        "CORR": float(np.corrcoef(pred, obs)[0, 1]),
    }


def map_stations(members, glat, glon, slat, slon):
    """Map member stations to their nearest member cell; others to nearest cell."""
    nearest = cKDTree(np.column_stack([glat, glon])).query(
        np.column_stack([slat, slon]), k=1)[1]
    owners = [[] for _ in range(len(slat))]
    for gi, sta in enumerate(members):
        for s in sta:
            owners[int(s)].append(gi)
    assigned = np.array([len(x) > 0 for x in owners])
    grid = nearest.copy()
    for s, candidates in enumerate(owners):
        if not candidates:
            continue
        candidates = np.asarray(candidates, dtype=int)
        d2 = ((glat[candidates] - slat[s]) ** 2 +
              (glon[candidates] - slon[s]) ** 2)
        grid[s] = candidates[np.argmin(d2)]
    distance = np.sqrt((glat[grid] - slat) ** 2 + (glon[grid] - slon) ** 2)
    return grid, assigned, distance
