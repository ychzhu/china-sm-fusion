# -*- coding: utf-8 -*-
"""Strict station-grouped 10-fold CV with fold-specific CDF and dynamic R."""
from paths_full import *
from cdf_utils import cdf_adjusting2
from oi_utils import build_geometry, dynamic_oi_day

import csv
import json
import os
import numpy as np

A, THRESH, MIN_HIST, SEED, NFOLDS = 10.0, 3, 30, 2026, 10

g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
nb = np.load(os.path.join(WORK, "neighbors_full.npz"), allow_pickle=True)
bearing, glat, glon = g["bearing"], g["glat"], g["glon"]
grid_obs = g["grid_obs"]
grid_err_daily = g["grid_err_daily"].astype(np.float64)
members = g["members"]
SMAP_ERR = float(g["smap_err"][0])
smap_bg = g["smap_bg"]
cap_idx, cap_dist = nb["cap_idx"], nb["cap_dist"]
nEff = len(bearing)
n_days, n_grid = smap_bg.shape

# Union all validation cells that share at least one archived station.
parent = np.arange(nEff)


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(x, y):
    rx, ry = find(x), find(y)
    if rx != ry:
        parent[ry] = rx


station_owner = {}
for v, gi in enumerate(bearing):
    for station in members[gi]:
        station = int(station)
        if station in station_owner:
            union(v, station_owner[station])
        else:
            station_owner[station] = v
components = {}
for v in range(nEff):
    components.setdefault(find(v), []).append(v)

# Balance held-out cell-day counts while keeping every component intact.
rng = np.random.RandomState(SEED)
keys = list(components)
rng.shuffle(keys)
keys.sort(key=lambda root: -int(np.isfinite(grid_obs[components[root]]).sum()))
fold_load = np.zeros(NFOLDS, dtype=np.int64)
fold_size = np.zeros(NFOLDS, dtype=int)
fold_of = np.full(nEff, -1, dtype=int)
component_of = np.full(nEff, -1, dtype=int)
for component_id, root in enumerate(keys):
    comp = components[root]
    load = int(np.isfinite(grid_obs[comp]).sum())
    k = int(np.lexsort((fold_size, fold_load))[0])
    fold_of[comp] = k
    component_of[comp] = component_id
    fold_load[k] += load
    fold_size[k] += len(comp)
assert np.all(fold_of >= 0)
print("fold grid sizes:", fold_size.tolist())
print("fold cell-day loads:", fold_load.tolist())

with open(os.path.join(WORK, "cv_fold_assignments.csv"), "w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh)
    writer.writerow(["grid_index", "bearing_position", "fold", "component", "station_indices"])
    for v, gi in enumerate(bearing):
        writer.writerow([int(gi), v, int(fold_of[v] + 1), int(component_of[v]),
                         ";".join(str(int(x)) for x in members[gi])])

# Prove that no station occurs on both sides of any fold.
for k in range(NFOLDS):
    val_stations = set(int(s) for v in np.where(fold_of == k)[0]
                       for s in members[bearing[v]])
    train_stations = set(int(s) for v in np.where(fold_of != k)[0]
                         for s in members[bearing[v]])
    assert val_stations.isdisjoint(train_stations), "station leakage in fold %d" % (k + 1)


def metrics(obs, *systems):
    out = {}
    for name, values in systems:
        d = values - obs
        out[name] = {
            "N": int(len(obs)),
            "MAE": float(np.mean(np.abs(d))),
            "RMSE": float(np.sqrt(np.mean(d ** 2))),
            "BIAS": float(np.mean(d)),
            "CORR": float(np.corrcoef(values, obs)[0, 1]),
        }
    return out


all_b, all_f, all_o = [], [], []
all_fold, all_grid, all_time = [], [], []
fold_metrics = {}
cdf_audit = []
for k in range(NFOLDS):
    V = np.where(fold_of == k)[0]
    fold_cdf = smap_bg.copy()
    cdf_applied = 0
    cdf_training_neighbors = []
    for gi in range(n_grid):
        valid = np.isfinite(smap_bg[:, gi])
        if valid.sum() <= MIN_HIST:
            continue
        ci = cap_idx[gi]
        ci = ci[ci >= 0]
        ci = ci[fold_of[ci] != k]
        if len(ci) < THRESH:
            continue
        fold_cdf[:, gi] = cdf_adjusting2(smap_bg[:, gi], grid_obs[ci, :].flatten())
        cdf_applied += 1
        cdf_training_neighbors.append(len(ci))
    cdf_audit.append({
        "fold": k + 1,
        "cdf_applied_grids": int(cdf_applied),
        "min_training_neighbors": int(min(cdf_training_neighbors, default=0)),
        "max_training_neighbors": int(max(cdf_training_neighbors, default=0)),
    })
    print("fold %d fold-specific CDF grids=%d" % (k + 1, cdf_applied), flush=True)

    grid_v = bearing[V]
    Bv = fold_cdf[:, grid_v].T.copy()
    Ov = grid_obs[V].copy()
    Fv = Bv.copy()
    active = []
    neighbor_lists = []
    distance_lists = []
    for row, v in enumerate(V):
        gi = bearing[v]
        ci = cap_idx[gi]
        ci = ci[ci >= 0]
        cd = cap_dist[gi][:len(ci)]
        train = fold_of[ci] != k
        ci, cd = ci[train], cd[train]
        if len(ci) >= THRESH and np.isfinite(fold_cdf[:, gi]).sum() > MIN_HIST:
            active.append(row)
            neighbor_lists.append(ci)
            distance_lists.append(cd)
    active = np.asarray(active, dtype=int)
    if len(active):
        bp_pad, b_mat, c_vec = build_geometry(
            glat, glon, bearing, grid_v[active], neighbor_lists, distance_lists, A)
        bp_safe = np.where(bp_pad >= 0, bp_pad, 0)
        for di in range(n_days):
            target_bg = fold_cdf[di, grid_v[active]]
            obs = grid_obs[bp_safe, di]
            neighbor_bg = fold_cdf[di, bearing[bp_safe]]
            neighbor_err = grid_err_daily[bp_safe, di]
            inc, _, _, _ = dynamic_oi_day(
                target_bg, obs, neighbor_bg, neighbor_err, bp_pad,
                b_mat, c_vec, SMAP_ERR)
            Fv[active, di] = target_bg + inc
        print("fold %d dynamic OI targets=%d" % (k + 1, len(active)), flush=True)

    ok = np.isfinite(Ov) & np.isfinite(Bv) & np.isfinite(Fv)
    fold_metrics[str(k + 1)] = metrics(Ov[ok], ("background", Bv[ok]), ("fusion", Fv[ok]))
    all_b.append(Bv[ok])
    all_f.append(Fv[ok])
    all_o.append(Ov[ok])
    all_fold.append(np.full(ok.sum(), k + 1, dtype=np.int16))
    all_grid.append(np.repeat(grid_v[:, None], n_days, axis=1)[ok])
    all_time.append(np.tile(np.arange(n_days, dtype=np.int32), (len(V), 1))[ok])
    print("fold %d pairs=%d MAE %.6f -> %.6f" %
          (k + 1, ok.sum(), fold_metrics[str(k + 1)]["background"]["MAE"],
           fold_metrics[str(k + 1)]["fusion"]["MAE"]), flush=True)

all_b = np.concatenate(all_b)
all_f = np.concatenate(all_f)
all_o = np.concatenate(all_o)
all_fold = np.concatenate(all_fold)
all_grid = np.concatenate(all_grid)
all_time = np.concatenate(all_time)
overall = metrics(all_o, ("background", all_b), ("fusion", all_f))
summary = {
    "validation": "strict station-grouped 10-fold CV",
    "fold_specific_cdf": True,
    "dynamic_daily_observation_error": True,
    "a": A,
    "radius_degrees": 1.08,
    "max_near": 9,
    "seed": SEED,
    "n_validation_grids": int(nEff),
    "n_components": int(len(components)),
    "station_overlap_between_train_and_validation": 0,
    "fold_grid_sizes": fold_size.tolist(),
    "fold_cell_day_loads": fold_load.tolist(),
    "cdf_audit": cdf_audit,
    "overall": overall,
    "by_fold": fold_metrics,
}
with open(os.path.join(WORK, "cv_station_grouped_summary.json"), "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)
with open(os.path.join(WORK, "cv_fold_metrics.csv"), "w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh)
    writer.writerow(["fold", "system", "N", "MAE", "RMSE", "BIAS", "CORR"])
    for fold, values in fold_metrics.items():
        for system, m in values.items():
            writer.writerow([fold, system, m["N"], m["MAE"], m["RMSE"], m["BIAS"], m["CORR"]])
np.savez(os.path.join(WORK, "cv_results_full.npz"),
         allB=all_b, allF=all_f, allO=all_o, sample_fold=all_fold,
         sample_grid=all_grid, sample_time=all_time, fold_of=fold_of,
         component_of=component_of)
print("===== strict station-grouped 10-fold CV =====")
print(json.dumps(overall, indent=2))
print("saved strict CV arrays, assignments, fold metrics and summary")
