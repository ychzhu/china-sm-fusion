# -*- coding: utf-8 -*-
"""Shared dynamic-R optimal-interpolation geometry and numerical routines."""
import numpy as np


def haversine(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dl = np.radians(lon2 - lon1)
    aa = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(aa)))


def build_geometry(glat, glon, bearing, target_grids, neighbor_lists,
                   distance_lists, a=10.0):
    """Build padded batched geometry for arbitrary target/neighbour subsets."""
    n_target = len(target_grids)
    max_near = max((len(x) for x in neighbor_lists), default=0)
    bp_pad = np.full((n_target, max_near), -1, dtype=np.int32)
    b_mat = np.zeros((n_target, max_near, max_near), dtype=np.float64)
    c_vec = np.zeros((n_target, max_near), dtype=np.float64)
    for i, (bp, dist) in enumerate(zip(neighbor_lists, distance_lists)):
        bp = np.asarray(bp, dtype=np.int32)
        dist = np.asarray(dist, dtype=np.float64)
        n = len(bp)
        bp_pad[i, :n] = bp
        la, lo = glat[bearing[bp]], glon[bearing[bp]]
        dmat = np.zeros((n, n), dtype=np.float64)
        for r in range(n):
            dmat[r] = haversine(la[r], lo[r], la, lo)
        b_mat[i, :n, :n] = np.power(a, -dmat)
        c_vec[i, :n] = np.power(a, -dist)
    return bp_pad, b_mat, c_vec


def dynamic_oi_day(target_bg, obs, neighbor_bg, neighbor_err, bp_pad,
                   b_mat, c_vec, sigma_b):
    """Return OI increment and standard posterior SD for one day.

    R is rebuilt from the actual stations contributing to each cell-day. Missing
    neighbours receive zero c and residual plus a very large R, so they cannot
    affect either the increment or posterior variance.
    """
    target_bg = np.asarray(target_bg, dtype=np.float64)
    valid_pad = bp_pad >= 0
    usable = (valid_pad & np.isfinite(obs) & np.isfinite(neighbor_bg) &
              np.isfinite(neighbor_err) & (neighbor_err >= 0) &
              np.isfinite(target_bg)[:, None])
    ratio = np.where(usable, neighbor_err / sigma_b, 1.0e6)
    system = b_mat.copy()
    diag = np.arange(system.shape[1])
    system[:, diag, diag] += ratio ** 2
    c_eff = np.where(usable, c_vec, 0.0)
    try:
        weights = np.linalg.solve(system, c_eff[..., None])[..., 0]
    except np.linalg.LinAlgError:
        weights = np.einsum("nij,nj->ni", np.linalg.pinv(system), c_eff)
    residual = np.where(usable, obs - neighbor_bg, 0.0)
    increment = np.sum(weights * residual, axis=1)
    reduction = np.sum(c_eff * weights, axis=1)
    posterior = sigma_b * np.sqrt(np.clip(1.0 - reduction, 0.0, 1.0))
    has_information = np.any(usable, axis=1)
    increment = np.where(has_information, increment, 0.0)
    increment = np.where(np.isfinite(target_bg), increment, np.nan)
    posterior = np.where(np.isfinite(target_bg), posterior, np.nan)
    return increment, posterior, has_information, reduction
