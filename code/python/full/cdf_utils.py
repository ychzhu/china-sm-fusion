# -*- coding: utf-8 -*-
"""Shared CDF-matching utilities used by production and fold-specific CV."""
import numpy as np


def interp1_extrap(x, y, xi):
    """Linear interpolation with linear extrapolation at both ends."""
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
    """Quadratic CDF adjustment; station reference is supplied by the caller."""
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
