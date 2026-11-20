# -*- coding: utf-8 -*-
"""Station and gridded comparisons with ERA5/GLDAS roles explicitly separated."""
from paths_full import *
from validation_utils import metrics, map_stations

import csv
import json
import os
import numpy as np
import netCDF4


def open_nc(path):
    with open(path, "rb") as fh:
        return netCDF4.Dataset(os.path.basename(path), memory=fh.read())


def load_ref(path, common_ymd):
    ds = open_nc(path)
    ymd = ds.variables["time_ymd"][:].astype(str)
    values = np.ma.filled(ds.variables["sm"][:], np.nan).astype(float)
    ds.close()
    pos = {d: i for i, d in enumerate(ymd)}
    return values[np.array([pos[d] for d in common_ymd])]


def write_rows(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
glat, glon, common_ymd, members = g["glat"], g["glon"], g["common_ymd"], g["members"]
smap_bg = g["smap_bg"]
smap_cdf = np.load(os.path.join(WORK, "smap_cdf_full.npy"))
fusion = np.load(os.path.join(WORK, "fusion_full.npy"))
era5 = load_ref(os.path.join(REFERENCE, "era5_gcs_swvl1_36km_20150331_20261001.nc"), common_ymd)
gldas = load_ref(os.path.join(REFERENCE, "gldas_noah_soilm_36km_20150331_20261001.nc"), common_ymd)

slat = np.load(os.path.join(WORK, "lat3596.npy"))
slon = np.load(os.path.join(WORK, "lon3596.npy"))
cleaned = np.load(os.path.join(WORK, "obj_clean3596.npy"))
sta_ymd = station_axis_ymd()
pos = {d: i for i, d in enumerate(sta_ymd)}
sta_obs = cleaned[:, np.array([pos[d] for d in common_ymd])]
s2g, assigned, _ = map_stations(members, glat, glon, slat, slon)
systems = {
    "GLDAS": gldas[:, s2g].T, "ERA5": era5[:, s2g].T,
    "SMAP_raw": smap_bg[:, s2g].T, "SMAP_CDF": smap_cdf[:, s2g].T,
    "Fusion": fusion[:, s2g].T,
}
rows = []
for scope, select in (("inside", assigned), ("outside", ~assigned), ("all", np.ones(len(slat), bool))):
    obs = sta_obs[select]
    for name, arr in systems.items():
        values = arr[select]
        valid = np.isfinite(values) & np.isfinite(obs)
        row = {"scope": scope, "system": name, "mask": "native pair availability"}
        row.update(metrics(values[valid], obs[valid]))
        rows.append(row)
write_rows(os.path.join(WORK, "reference_station_metrics_by_scope.csv"),
           ["scope", "system", "mask", "N", "MAE", "RMSE", "BIAS", "CORR"], rows)

rows = []
for scope, select in (("inside", assigned), ("outside", ~assigned), ("all", np.ones(len(slat), bool))):
    obs = sta_obs[select]
    joint = np.isfinite(obs)
    for arr in systems.values():
        joint &= np.isfinite(arr[select])
    for name, arr in systems.items():
        values = arr[select]
        row = {"scope": scope, "system": name, "mask": "joint five-system availability"}
        row.update(metrics(values[joint], obs[joint]))
        rows.append(row)
write_rows(os.path.join(WORK, "reference_station_metrics_joint.csv"),
           ["scope", "system", "mask", "N", "MAE", "RMSE", "BIAS", "CORR"], rows)

# Native in-grid monthly and seasonal reference metrics for the main figures.
# Each system retains its own finite-pair availability; no joint mask is imposed.
obs_inside = sta_obs[assigned]
month = np.array([int(d[4:6]) for d in common_ymd])
season = np.array(["Winter (DJF)" if m == 12 or m <= 2 else
                   "Spring (MAM)" if m <= 5 else
                   "Summer (JJA)" if m <= 8 else "Autumn (SON)" for m in month])
for filename, values, labels in (
        ("reference_station_metrics_monthly.csv", month, list(range(1, 13))),
        ("reference_station_metrics_seasonal.csv", season,
         ["Winter (DJF)", "Spring (MAM)", "Summer (JJA)", "Autumn (SON)"])):
    temporal_rows = []
    for label in labels:
        time_mask = values == label
        for name, arr in systems.items():
            pred = arr[assigned][:, time_mask]
            obs = obs_inside[:, time_mask]
            valid = np.isfinite(pred) & np.isfinite(obs)
            row = {"Period": label, "scope": "inside", "system": name,
                   "mask": "native pair availability"}
            row.update(metrics(pred[valid], obs[valid]))
            temporal_rows.append(row)
    write_rows(os.path.join(WORK, filename),
               ["Period", "scope", "system", "mask", "N", "MAE", "RMSE", "BIAS", "CORR"],
               temporal_rows)

# Gridded reference comparisons: BIAS = reference - product.
pairs = [("ERA5_vs_SMAP_raw", era5, smap_bg), ("ERA5_vs_Fusion", era5, fusion),
         ("GLDAS_vs_SMAP_raw", gldas, smap_bg), ("GLDAS_vs_Fusion", gldas, fusion),
         ("SMAP_raw_vs_Fusion", smap_bg, fusion)]
rows, pergrid = [], []
for name, ref, product in pairs:
    valid = np.isfinite(ref) & np.isfinite(product)
    diff = ref[valid] - product[valid]
    row = {"Pair": name, "N": int(valid.sum()),
           "MAE": float(np.mean(np.abs(diff))),
           "RMSE": float(np.sqrt(np.mean(diff ** 2))),
           "BIAS": float(np.mean(diff)),
           "CORR": float(np.corrcoef(ref[valid], product[valid])[0, 1])}
    rows.append(row)
    cell_metrics = []
    for gi in range(ref.shape[1]):
        ok = valid[:, gi]
        if ok.sum() < 60:
            continue
        d = ref[ok, gi] - product[ok, gi]
        cell_metrics.append((np.mean(np.abs(d)), np.sqrt(np.mean(d ** 2)),
                             np.mean(d), np.corrcoef(ref[ok, gi], product[ok, gi])[0, 1]))
    cell_metrics = np.asarray(cell_metrics)
    pergrid.append({"Pair": name, "ValidGrids": int(len(cell_metrics)),
                    "MAE_median": float(np.median(cell_metrics[:, 0])),
                    "RMSE_median": float(np.median(cell_metrics[:, 1])),
                    "BIAS_median": float(np.median(cell_metrics[:, 2])),
                    "CORR_median": float(np.median(cell_metrics[:, 3]))})
write_rows(os.path.join(WORK, "stats_smap_fusion.csv"),
           ["Pair", "N", "MAE", "RMSE", "BIAS", "CORR"], rows)
write_rows(os.path.join(WORK, "stats_smap_fusion_pergrid.csv"),
           ["Pair", "ValidGrids", "MAE_median", "RMSE_median", "BIAS_median", "CORR_median"], pergrid)
roles = {
    "ERA5": "third-party reanalysis; used in TCA calibration, therefore not a fully independent validation reference",
    "GLDAS": "third-party reanalysis; not used in TCA calibration or fusion; relative independence is clearer",
    "validation_masks": "native pair availability and joint five-system availability are reported separately",
}
with open(os.path.join(WORK, "reference_roles.json"), "w", encoding="utf-8") as fh:
    json.dump(roles, fh, ensure_ascii=False, indent=2)
print("saved station native/joint metrics and gridded reference tables")
