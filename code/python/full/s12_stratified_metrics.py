# -*- coding: utf-8 -*-
"""Unified V2 annual, seasonal, density, vegetation and regional summaries."""
from paths_full import *
from validation_utils import metrics, map_stations

import csv
import json
import os
from datetime import datetime
import numpy as np
import pandas as pd
import netCDF4


def write_csv(name, rows):
    if not rows:
        return
    with open(os.path.join(WORK, name), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def row_base(**kwargs):
    return dict(kwargs)


def add_metrics(row, pred, obs, prefix=""):
    valid = np.isfinite(pred) & np.isfinite(obs)
    if valid.sum() < 2:
        row.update({prefix + "N": int(valid.sum()), prefix + "MAE": np.nan,
                    prefix + "RMSE": np.nan, prefix + "BIAS": np.nan, prefix + "CORR": np.nan})
        return
    row.update({prefix + k: v for k, v in metrics(pred[valid], obs[valid]).items()})


def load_reference(path, common_ymd):
    ds = open_nc(path)
    ymd = ds.variables["time_ymd"][:].astype(str)
    sm = np.ma.filled(ds.variables["sm"][:], np.nan).astype(float)
    ds.close()
    pos = {d: i for i, d in enumerate(ymd)}
    return sm[np.array([pos[d] for d in common_ymd])]


g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
glat, glon, bearing = g["glat"], g["glon"], g["bearing"]
common_ymd = g["common_ymd"].astype(str)
grid_obs = g["grid_obs"].T
smap_cdf = np.load(os.path.join(WORK, "smap_cdf_full.npy"))
fusion = np.load(os.path.join(WORK, "fusion_full.npy"))
slat = np.load(os.path.join(WORK, "lat3596.npy"))
slon = np.load(os.path.join(WORK, "lon3596.npy"))
cleaned = np.load(os.path.join(WORK, "obj_clean3596.npy"))
sta_ymd = station_axis_ymd()
pos = {d: i for i, d in enumerate(sta_ymd)}
sta_obs = cleaned[:, np.array([pos[d] for d in common_ymd])]
s2g, assigned, _ = map_stations(g["members"], glat, glon, slat, slon)
sta_fusion = fusion[:, s2g].T
sta_cdf = smap_cdf[:, s2g].T
era5 = load_reference(os.path.join(REFERENCE, "era5_gcs_swvl1_36km_20150331_20261001.nc"), common_ymd)
sta_era5 = era5[:, s2g].T

dates = np.array([datetime.strptime(str(d), "%Y%m%d") for d in common_ymd])
years = np.array([d.year for d in dates])
months = np.array([d.month for d in dates])
seasons = np.array(["Winter (DJF)" if m in (12, 1, 2) else
                    "Spring (MAM)" if m in (3, 4, 5) else
                    "Summer (JJA)" if m in (6, 7, 8) else "Autumn (SON)" for m in months])

# Station metadata in the released NetCDF and NDVI CSV.
ds = open_nc(STATION3596_NC)
station_id = np.array([str(x) for x in netCDF4.chartostring(ds.variables["station_id"][:])])
province = np.array([str(x).strip() for x in netCDF4.chartostring(ds.variables["province"][:])])
ds.close()
ndvi_table = pd.read_csv(os.path.join(DATA, "ndvi", "station_ndvi.csv"), dtype={"Station_ID": str})
ndvi_table["Station_ID"] = ndvi_table["Station_ID"].str.strip()
ndvi_map = ndvi_table.set_index("Station_ID")["NDVI_growing"].to_dict()
ndvi = np.array([ndvi_map.get(s, np.nan) for s in station_id], dtype=float)

# Number of station-bearing cells inside the operational radius at each cell.
neighbor_count = np.zeros(len(glat), dtype=int)
blat, blon = glat[bearing], glon[bearing]
for gi in bearing:
    cand_local = np.where((np.abs(blat - glat[gi]) <= 1.08) &
                          (np.abs(blon - glon[gi]) <= 1.08))[0]
    p1 = np.radians(glat[gi]); p2 = np.radians(blat[cand_local])
    dp = np.radians(blat[cand_local] - glat[gi]); dl = np.radians(blon[cand_local] - glon[gi])
    h = np.sin(dp/2)**2 + np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2
    dist = np.degrees(2*np.arcsin(np.sqrt(np.clip(h, 0, 1))))
    neighbor_count[gi] = int((dist <= 1.08 + 1e-12).sum())
density = np.array([neighbor_count[c] if c >= 0 else 0 for c in s2g])
density_labels = ["Sparse (0–3)", "Moderate (4–6)", "High (7–9)", "Very high (≥10)"]
density_masks = [density <= 3, (density >= 4) & (density <= 6),
                 (density >= 7) & (density <= 9), density >= 10]
ndvi_labels = ["V1 sparse", "V2 grass/crop", "V3 mixed", "V4 forest"]
ndvi_masks = [ndvi < 0.2, (ndvi >= 0.2) & (ndvi < 0.4),
              (ndvi >= 0.4) & (ndvi < 0.6), ndvi >= 0.6]

# Annual station-level results; ERA5 is paired to the fusion-valid samples.
qc = pd.read_csv(os.path.join(STATION, "qc_station_year_3596.csv"))
annual = []
for year in range(2015, 2027):
    tm = years == year
    r = {"Year": year,
         "Raw stations": int((qc.year == year).sum()),
         "QC-passing": int(((qc.year == year) & (qc.decision == "keep")).sum())}
    pred, obs = sta_fusion[assigned][:, tm], sta_obs[assigned][:, tm]
    add_metrics(r, pred, obs, "Fusion_")
    valid = np.isfinite(pred) & np.isfinite(obs)
    ep, eo = sta_era5[assigned][:, tm], sta_obs[assigned][:, tm]
    add_metrics(r, np.where(valid, ep, np.nan), eo, "ERA5_")
    annual.append(r)
write_csv("station_annual_metrics_inside.csv", annual)

season_rows = []
for label in ["Winter (DJF)", "Spring (MAM)", "Summer (JJA)", "Autumn (SON)"]:
    tm = seasons == label
    r = {"Season": label}
    add_metrics(r, sta_fusion[assigned][:, tm], sta_obs[assigned][:, tm])
    season_rows.append(r)
write_csv("station_seasonal_metrics_inside.csv", season_rows)

month_rows = []
for month in range(1, 13):
    tm = months == month
    r = {"Month": month}
    add_metrics(r, sta_fusion[assigned][:, tm], sta_obs[assigned][:, tm], "Fusion_")
    add_metrics(r, sta_cdf[assigned][:, tm], sta_obs[assigned][:, tm], "SMAP_CDF_")
    month_rows.append(r)
write_csv("station_monthly_metrics_inside.csv", month_rows)

for filename, labels, masks, kind in [
        ("station_density_metrics_inside.csv", density_labels, density_masks, "Density class"),
        ("station_ndvi_metrics_inside.csv", ndvi_labels, ndvi_masks, "Vegetation class")]:
    rows = []
    for label, sm in zip(labels, masks):
        select = sm & assigned
        r = {kind: label, "Stations": int(sm.sum()), "Stations with valid pairs": int(select.sum())}
        add_metrics(r, sta_fusion[select], sta_obs[select])
        rows.append(r)
    write_csv(filename, rows)

province_rows = []
for name in sorted(set(province)):
    select = assigned & (province == name)
    r = {"Province": name, "Stations": int((province == name).sum()),
         "Stations with valid pairs": int(select.sum())}
    add_metrics(r, sta_fusion[select], sta_obs[select])
    province_rows.append(r)
write_csv("station_province_metrics_inside.csv", province_rows)

# Grid-cell in-sample summaries for the corresponding main/SI figures.
grid_annual, grid_seasonal, grid_monthly, grid_daily = [], [], [], []
bg_b, fu_b = smap_cdf[:, bearing], fusion[:, bearing]
for year in range(2015, 2027):
    tm = years == year
    r = {"Year": year}
    add_metrics(r, bg_b[tm], grid_obs[tm], "SMAP_CDF_")
    add_metrics(r, fu_b[tm], grid_obs[tm], "Fusion_")
    grid_annual.append(r)
for month in range(1, 13):
    tm = months == month
    r = {"Month": month}
    add_metrics(r, bg_b[tm], grid_obs[tm], "SMAP_CDF_")
    add_metrics(r, fu_b[tm], grid_obs[tm], "Fusion_")
    grid_monthly.append(r)
for label in ["Winter (DJF)", "Spring (MAM)", "Summer (JJA)", "Autumn (SON)"]:
    tm = seasons == label
    r = {"Season": label}
    add_metrics(r, bg_b[tm], grid_obs[tm], "SMAP_CDF_")
    add_metrics(r, fu_b[tm], grid_obs[tm], "Fusion_")
    grid_seasonal.append(r)
for di, day in enumerate(common_ymd):
    r = {"Date": day}
    add_metrics(r, bg_b[di], grid_obs[di], "SMAP_CDF_")
    add_metrics(r, fu_b[di], grid_obs[di], "Fusion_")
    grid_daily.append(r)
write_csv("gridcell_annual_metrics.csv", grid_annual)
write_csv("gridcell_monthly_metrics.csv", grid_monthly)
write_csv("gridcell_seasonal_metrics.csv", grid_seasonal)
write_csv("gridcell_daily_metrics.csv", grid_daily)

# Strict-CV temporal summaries.
cv = np.load(os.path.join(WORK, "cv_results_full.npz"))
all_b, all_f, all_o = cv["allB"], cv["allF"], cv["allO"]
time_index = cv["sample_time"]
cv_annual, cv_seasonal = [], []
for year in range(2015, 2027):
    select = years[time_index] == year
    r = {"Year": year}
    add_metrics(r, all_b[select], all_o[select], "SMAP_CDF_")
    add_metrics(r, all_f[select], all_o[select], "Fusion_")
    cv_annual.append(r)
for label in ["Winter (DJF)", "Spring (MAM)", "Summer (JJA)", "Autumn (SON)"]:
    select = seasons[time_index] == label
    r = {"Season": label}
    add_metrics(r, all_b[select], all_o[select], "SMAP_CDF_")
    add_metrics(r, all_f[select], all_o[select], "Fusion_")
    cv_seasonal.append(r)
write_csv("cv_annual_metrics.csv", cv_annual)
write_csv("cv_seasonal_metrics.csv", cv_seasonal)

counts = {
    "archived_stations": len(station_id),
    "inside_station_identifiers": int(assigned.sum()),
    "outside_station_identifiers": int((~assigned).sum()),
    "density_station_sum": int(sum(int(r["Stations"]) for r in csv.DictReader(open(os.path.join(WORK, "station_density_metrics_inside.csv"), encoding="utf-8")))),
    "ndvi_station_sum": int(sum(int(r["Stations"]) for r in csv.DictReader(open(os.path.join(WORK, "station_ndvi_metrics_inside.csv"), encoding="utf-8")))),
    "missing_ndvi": int((~np.isfinite(ndvi)).sum()),
    "annual_pair_sum_fusion": int(sum(r["Fusion_N"] for r in annual)),
    "seasonal_pair_sum_fusion": int(sum(r["N"] for r in season_rows)),
}
with open(os.path.join(WORK, "stratification_counts.json"), "w", encoding="utf-8") as fh:
    json.dump(counts, fh, indent=2)
print(json.dumps(counts, indent=2))
