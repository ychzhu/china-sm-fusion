# -*- coding: utf-8 -*-
"""Stage 1: objective hard-fault screening on the full 3596 domain (TCA-independent).

Per-year neighbor filtering: the active-window and spatial checks use only
stations that actually appear in that year, so later-installed stations do not
dilute the neighbor reporting fraction.
"""
from paths_full import *

import os
import numpy as np
import pandas as pd
import netCDF4
from scipy.spatial import cKDTree


NEIGH_RADIUS = 0.75
MIN_NEIGH = 3
SEASON_MIN_DAYS = 10
RANGE_LO, RANGE_HI = 0.0, 0.6
RANGE_BAD_FRAC = 0.05
FLAT_STD = 0.005
NEIGH_ACTIVE_STD = 0.01
JUMP_ABS = 0.15
JUMP_NEIGH = 0.05
JUMP_COUNT = 3
COMP_FLOOR = 0.50
L2_MEAN_LO, L2_MEAN_HI = -0.10, 0.15
R_THRESH = 0.5
KEEP, MASK, DROP = 0, 1, 2

# Synthetic / test identifiers that are not assigned CMA station numbers.
# Excluded by an objective identifier-validity (provenance) screen, independent
# of any error estimate. These are the 1-digit codes and the all-repeated-digit
# 5-digit codes present in the raw archive.
NON_STATION_IDS = {1, 2, 3, 11111, 22222, 33333, 55555, 66666, 99999}

SEASONS = {"DJF": (12, 1, 2), "MAM": (3, 4, 5),
           "JJA": (6, 7, 8), "SON": (9, 10, 11)}
GROWING = {"MAM", "JJA", "SON"}

ymd = station_axis_ymd()
n_day = len(ymd)
year = np.array([int(d[:4]) for d in ymd])
month = np.array([int(d[4:6]) for d in ymd])
years = np.arange(year.min(), year.max() + 1)

raw = np.load(os.path.join(WORK, "raw3596.npy"))
lat = np.load(os.path.join(WORK, "lat3596.npy"))
lon = np.load(os.path.join(WORK, "lon3596.npy"))
ids = np.load(os.path.join(WORK, "ids3596.npy"))
present = np.load(os.path.join(WORK, "present_year.npy"))
n_sta = len(ids)

global_neigh = cKDTree(np.column_stack([lat, lon])).query_ball_point(
    np.column_stack([lat, lon]), NEIGH_RADIUS)

obj_flag = np.zeros((n_sta, n_day), dtype=np.int8)
rows = []
for yi, yr in enumerate(years):
    ymask = year == yr
    yidx = np.where(ymask)[0]
    for i in range(n_sta):
        if not present[i, yi]:
            continue
        non_station = int(ids[i]) in NON_STATION_IDS
        # year-restricted neighbors (present this year, excluding self)
        nbr = [j for j in global_neigh[i] if j != i and present[j, yi]]
        x = raw[i, ymask].astype(float)
        reported = np.isfinite(x)
        n_reported = int(reported.sum())
        has_nbr = len(nbr) > MIN_NEIGH
        day_flag = np.full(len(yidx), KEEP, dtype=np.int8)
        reasons = []

        if has_nbr:
            nbr_rep = np.isfinite(raw[np.array(nbr)][:, ymask])
            expected = reported | (nbr_rep.mean(axis=0) >= 0.5)
            n_expected = int(expected.sum())
            outage = int((expected & ~reported).sum())
            comp = n_reported / n_expected if n_expected else 0.0
        else:
            n_expected = len(yidx)
            outage = int((~reported).sum())
            comp = n_reported / n_expected if n_expected else 0.0
        l0_fail = n_expected >= 10 and comp < COMP_FLOOR

        l1_fail = False
        growing_stuck = False
        nmean_full = None
        if has_nbr:
            with np.errstate(all="ignore"):
                nmean_full = np.nanmean(raw[np.array(nbr)][:, ymask], axis=0)
        for sn, mlist in SEASONS.items():
            msel = np.isin(month[ymask], mlist)
            if msel.sum() == 0:
                continue
            xs = x[msel]
            valid = np.isfinite(xs)
            if valid.sum() < SEASON_MIN_DAYS:
                continue
            xv = xs[valid]
            sidx_local = np.where(msel)[0][valid]

            bad = (xv < RANGE_LO) | (xv > RANGE_HI)
            if bad.any():
                day_flag[sidx_local[bad]] = MASK
                if bad.sum() / valid.sum() > RANGE_BAD_FRAC:
                    l1_fail = True
                    reasons.append("range_bad")

            if xv.std() < FLAT_STD:
                if sn == "DJF":
                    day_flag[sidx_local] = MASK
                    reasons.append("winter_freeze")
                elif sn in GROWING:
                    if has_nbr:
                        if np.nanstd(nmean_full[msel][valid]) > NEIGH_ACTIVE_STD:
                            growing_stuck = True
                            reasons.append("growing_stuck_" + sn)
                    else:
                        growing_stuck = True
                        reasons.append("growing_stuck_" + sn)

            if has_nbr:
                jump = (np.abs(np.diff(x, prepend=np.nan)) > JUMP_ABS) & \
                       (np.abs(np.diff(nmean_full, prepend=np.nan)) < JUMP_NEIGH)
                if jump.sum() >= JUMP_COUNT:
                    l1_fail = True
                    reasons.append("jump_bad")
                else:
                    day_flag[np.where(jump)[0]] = MASK

        l2_fail = False
        if has_nbr:
            with np.errstate(all="ignore"):
                mean3 = np.nanmean(x) - np.nanmean(nmean_full)
            if not (L2_MEAN_LO <= mean3 <= L2_MEAN_HI):
                bothm = reported & np.isfinite(nmean_full)
                if bothm.sum() >= 2 and np.std(x[bothm]) > 0 and np.std(nmean_full[bothm]) > 0:
                    if np.corrcoef(x[bothm], nmean_full[bothm])[0, 1] < R_THRESH:
                        l2_fail = True
                        reasons.append("spatial_bad")
                elif np.abs(mean3) > 0.30:
                    l2_fail = True
                    reasons.append("spatial_bad")

        if non_station:
            reasons.append("non_station_id")
        drop = non_station or l0_fail or l1_fail or l2_fail or growing_stuck
        if l0_fail:
            reasons.append("completeness_bad")
        if drop:
            day_flag[:] = DROP
        obj_flag[i, ymask] = day_flag
        rows.append({"cma_id": int(ids[i]), "year": int(yr),
                     "n_reported": n_reported, "expected_days": n_expected,
                     "individual_outage": outage,
                     "active_completeness": round(comp, 4),
                     "masked_days": int((day_flag == MASK).sum()),
                     "decision": "drop" if drop else "keep",
                     "reasons": ";".join(dict.fromkeys(reasons))})

res = pd.DataFrame(rows)
res.to_csv(os.path.join(WORK, "obj_station_year.csv"), index=False)
np.save(os.path.join(WORK, "obj_flag3596.npy"), obj_flag)
clean = np.where(obj_flag == KEEP, raw, np.nan).astype(np.float32)
np.save(os.path.join(WORK, "obj_clean3596.npy"), clean)

summ = []
for yr in years:
    sub = res[res["year"] == yr]
    rd = (sub["decision"] == "drop").sum()
    summ.append({"year": int(yr), "n_present": len(sub),
                 "n_keep": int((sub["decision"] == "keep").sum()),
                 "n_drop": int(rd),
                 "drop_pct": round(100 * rd / len(sub), 1),
                 "masked_days": int(sub["masked_days"].sum())})
summ = pd.DataFrame(summ)
summ.to_csv(os.path.join(WORK, "obj_yearly_summary.csv"), index=False)
print(summ.to_string(index=False))
print("TOTAL station-years=%d keep=%d drop=%d (%.1f%%)" %
      (len(res), (res.decision == 'keep').sum(), (res.decision == 'drop').sum(),
       100 * (res.decision == 'drop').mean()))
print("masked station-days=%d" % res.masked_days.sum())
print("clean-for-TCA finite=%d" % np.isfinite(clean).sum())
