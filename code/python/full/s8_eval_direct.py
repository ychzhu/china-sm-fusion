# -*- coding: utf-8 -*-
"""Direct validation with explicit in-sample and station-grid support scopes."""
from paths_full import *
from validation_utils import metrics, map_stations

import csv
import json
import os
import numpy as np


def write_rows(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
bearing, members = g["bearing"], g["members"]
glat, glon, common_ymd = g["glat"], g["glon"], g["common_ymd"]
grid_obs = g["grid_obs"]
smap_cdf = np.load(os.path.join(WORK, "smap_cdf_full.npy"))
fusion = np.load(os.path.join(WORK, "fusion_full.npy"))
fusion_err = np.load(os.path.join(WORK, "fusion_full_err.npy"))

rows = []
for name, pred in (("SMAP_CDF", smap_cdf[:, bearing].T), ("Fusion", fusion[:, bearing].T)):
    valid = np.isfinite(grid_obs) & np.isfinite(pred)
    row = {"scope": "station-bearing cell-days (in-sample)", "system": name}
    row.update(metrics(pred[valid], grid_obs[valid]))
    rows.append(row)
write_rows(os.path.join(WORK, "direct_gridcell_metrics.csv"),
           ["scope", "system", "N", "MAE", "RMSE", "BIAS", "CORR"], rows)

slat = np.load(os.path.join(WORK, "lat3596.npy"))
slon = np.load(os.path.join(WORK, "lon3596.npy"))
cleaned = np.load(os.path.join(WORK, "obj_clean3596.npy"))
sta_ymd = station_axis_ymd()
sta_pos = {d: i for i, d in enumerate(sta_ymd)}
sta_obs = cleaned[:, np.array([sta_pos[d] for d in common_ymd])]
s2g, assigned, distance = map_stations(members, glat, glon, slat, slon)
np.savez_compressed(os.path.join(WORK, "station_grid_mapping.npz"),
                    grid_index=s2g, assigned=assigned, distance_degrees=distance)
systems = {"SMAP_CDF": smap_cdf[:, s2g], "Fusion": fusion[:, s2g]}
rows = []
per_station = []
for scope, select in (("inside", assigned), ("outside", ~assigned), ("all", np.ones(len(slat), bool))):
    for name, arr in systems.items():
        values = arr.T[select]
        obs = sta_obs[select]
        valid = np.isfinite(values) & np.isfinite(obs)
        row = {"scope": scope, "system": name}
        row.update(metrics(values[valid], obs[valid]))
        rows.append(row)
    values = systems["Fusion"].T
    for s in np.where(select)[0]:
        valid = np.isfinite(values[s]) & np.isfinite(sta_obs[s])
        if valid.sum() < 30:
            continue
        row = {"station_index": int(s), "grid_index": int(s2g[s]),
               "scope": "inside" if assigned[s] else "outside"}
        row.update(metrics(values[s, valid], sta_obs[s, valid]))
        per_station.append(row)
write_rows(os.path.join(WORK, "station_scope_metrics.csv"),
           ["scope", "system", "N", "MAE", "RMSE", "BIAS", "CORR"], rows)
write_rows(os.path.join(WORK, "station_fusion_metrics.csv"),
           ["station_index", "grid_index", "scope", "N", "MAE", "RMSE", "BIAS", "CORR"], per_station)
counts = {
    "archived_station_identifiers": int(len(slat)),
    "inside_station_identifiers": int(assigned.sum()),
    "outside_station_identifiers": int((~assigned).sum()),
    "inside_valid_station_days": int(np.isfinite(sta_obs[assigned]).sum()),
    "outside_valid_station_days": int(np.isfinite(sta_obs[~assigned]).sum()),
    "outside_distance_median_degrees": float(np.median(distance[~assigned])),
    "outside_distance_max_degrees": float(np.max(distance[~assigned])),
}
with open(os.path.join(WORK, "station_scope_counts.json"), "w", encoding="utf-8") as fh:
    json.dump(counts, fh, ensure_ascii=False, indent=2)
valid_err = fusion_err[np.isfinite(fusion_err)]
aerr = {
    "N": int(valid_err.size), "median": float(np.median(valid_err)),
    "mean": float(np.mean(valid_err)),
    "P10": float(np.percentile(valid_err, 10)),
    "P90": float(np.percentile(valid_err, 90)),
    "above_background": int(np.sum(valid_err > 0.046 + 1e-12)),
}
with open(os.path.join(WORK, "aerr_summary.json"), "w", encoding="utf-8") as fh:
    json.dump(aerr, fh, indent=2)
print(json.dumps({"gridcell": rows[:2], "station_counts": counts, "Aerr": aerr}, indent=2))
