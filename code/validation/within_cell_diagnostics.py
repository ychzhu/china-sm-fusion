# -*- coding: utf-8 -*-
"""Within-cell station-pair diagnostics without an equal-sigma correlation inversion.

For pair (i,j), report RMS(x_i-x_j) relative to the independent-TCA prediction
sqrt(sigma_i^2+sigma_j^2). A positive residual can reflect correlated station
errors and/or point-scale representativeness; these terms are not identifiable
from this archive alone. Uncertainty is bootstrapped by grid cell, not by pair.
"""
from pathlib import Path
import argparse
import csv
import json
import os

import netCDF4
import numpy as np


def open_nc(path):
    with open(path, "rb") as fh:
        return netCDF4.Dataset(os.path.basename(path), memory=fh.read())


def default_release():
    if os.environ.get("SM_RELEASE_DIR"):
        return Path(os.environ["SM_RELEASE_DIR"])
    for p in (Path.cwd(), *Path(__file__).resolve().parents):
        if (p / "data" / "fusion").is_dir() and (p / "data" / "station").is_dir():
            return p
        candidate = p / "datarelease-scienceDB"
        if (candidate / "data" / "fusion").is_dir():
            return candidate
    raise FileNotFoundError("Use --release-root or SM_RELEASE_DIR")


def metrics(obs, pred):
    d = pred - obs
    return {
        "N": int(d.size),
        "MAE": float(np.mean(np.abs(d))),
        "RMSE": float(np.sqrt(np.mean(d ** 2))),
        "BIAS": float(np.mean(d)),
        "CORR": float(np.corrcoef(pred, obs)[0, 1]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-root", type=Path, default=default_release())
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()
    release = args.release_root.resolve()
    out = (args.output_dir or (release / "validation")).resolve()
    out.mkdir(parents=True, exist_ok=True)

    fusion = release / "data" / "fusion" / "fusion_sm_36km_3596_20150331_20261001.nc"
    station = release / "data" / "station" / "staD10_qc_3596_20150101_20261001.nc"
    with open_nc(fusion) as ds:
        glat = np.asarray(ds.variables["lat"][:], dtype=float)
        glon = np.asarray(ds.variables["lon"][:], dtype=float)
    with open_nc(station) as ds:
        x = np.ma.filled(ds.variables["sm"][:], np.nan).astype(float)
        slat = np.asarray(ds.variables["lat"][:], dtype=float)
        slon = np.asarray(ds.variables["lon"][:], dtype=float)
    if x.shape[0] != len(slat):
        x = x.T

    sigma = np.full(len(slat), 0.047, dtype=float)
    with (release / "data" / "station" / "station_tca_error_3596.csv").open(newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            value = float(row["station_sigma"])
            if np.isfinite(value) and value > 0:
                sigma[int(row["station_index"])] = value

    cell_of = np.full(len(slat), -1, dtype=int)
    for s in range(len(slat)):
        hit = np.flatnonzero((np.abs(glat - slat[s]) <= 0.18) &
                             (np.abs(glon - slon[s]) <= 0.18))
        if hit.size:
            cell_of[s] = int(hit[0])

    records = []
    for cell in np.unique(cell_of[cell_of >= 0]):
        stations = np.flatnonzero(cell_of == cell)
        if len(stations) < 2:
            continue
        for ai in range(len(stations) - 1):
            i = stations[ai]
            for j in stations[ai + 1:]:
                valid = np.isfinite(x[i]) & np.isfinite(x[j])
                n = int(valid.sum())
                if n < 30:
                    continue
                rms = float(np.sqrt(np.mean((x[i, valid] - x[j, valid]) ** 2)))
                p1, p2 = np.radians(slat[i]), np.radians(slat[j])
                dp = np.radians(slat[j] - slat[i])
                dl = np.radians(slon[j] - slon[i])
                aa = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
                sep = float(np.degrees(2 * np.arcsin(np.sqrt(aa))) * 111.2)
                independent = float(np.hypot(sigma[i], sigma[j]))
                rho0 = float((sigma[i] ** 2 + sigma[j] ** 2 - rms ** 2) /
                             (2 * sigma[i] * sigma[j]))
                records.append((int(cell), int(i), int(j), n, sep, rms,
                                float(sigma[i]), float(sigma[j]), independent,
                                rms / independent, rho0,
                                float(np.sqrt(max(0.0, rms ** 2 - independent ** 2)))))

    fields = ["grid_index0", "station_i0", "station_j0", "shared_days",
              "separation_km", "rms_difference", "sigma_i", "sigma_j",
              "independent_tca_rms", "standardized_rms", "rho_if_no_representativeness",
              "representativeness_residual_if_rho0"]
    with (out / "within_cell_pairs.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(fields)
        writer.writerows(records)
    arr = np.asarray(records, dtype=float)
    cell = arr[:, 0].astype(int)
    sep, rms = arr[:, 4], arr[:, 5]
    ratio = arr[:, 9]

    def summary(mask):
        return {
            "n_pairs": int(mask.sum()),
            "median_rms_difference": float(np.median(rms[mask])),
            "median_independent_tca_rms": float(np.median(arr[mask, 8])),
            "median_standardized_rms": float(np.median(ratio[mask])),
            "median_representativeness_residual_if_rho0": float(np.median(arr[mask, 11])),
            "median_rho_if_no_representativeness": float(np.median(arr[mask, 10])),
        }

    distance_rows = []
    edges = [0, 5, 10, 20, 30, 50, 80, 200]
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (sep >= lo) & (sep < hi)
        if mask.sum() < 20:
            continue
        row = {"distance_lo_km": lo, "distance_hi_km": hi, **summary(mask)}
        distance_rows.append(row)
    with (out / "within_cell_distance.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(distance_rows[0]))
        writer.writeheader()
        writer.writerows(distance_rows)

    unique_cells = np.unique(cell)
    groups = [np.flatnonzero(cell == c) for c in unique_cells]
    rng = np.random.RandomState(2026)
    boot_ratio = np.empty(args.bootstrap)
    boot_rms = np.empty(args.bootstrap)
    for b in range(args.bootstrap):
        selected = rng.randint(0, len(groups), len(groups))
        idx = np.concatenate([groups[k] for k in selected])
        boot_ratio[b] = np.median(ratio[idx])
        boot_rms[b] = np.median(rms[idx])
    closest = np.argsort(sep)[:min(500, len(sep))]
    closest_mask = np.zeros(len(sep), dtype=bool)
    closest_mask[closest] = True
    result = {
        "model": "RMS pair difference relative to sqrt(sigma_i^2+sigma_j^2)",
        "interpretation": "Station-error correlation and point-scale representativeness are not separately identifiable; rho_if_no_representativeness is diagnostic only.",
        "n_pairs": len(records),
        "n_cells": len(unique_cells),
        "median_shared_days": float(np.median(arr[:, 3])),
        "median_separation_km": float(np.median(sep)),
        **summary(np.ones(len(sep), dtype=bool)),
        "closest_500_pairs": summary(closest_mask),
        "closest_500_median_separation_km": float(np.median(sep[closest_mask])),
        "fraction_rho_if_no_representativeness_outside_bounds": float(np.mean((arr[:, 10] < -1) | (arr[:, 10] > 1))),
        "cluster_bootstrap": {
            "unit": "grid cell", "B": args.bootstrap, "interval": "95% percentile",
            "median_standardized_rms_ci": np.percentile(boot_ratio, [2.5, 97.5]).tolist(),
            "median_rms_difference_ci": np.percentile(boot_rms, [2.5, 97.5]).tolist(),
        },
    }
    with (out / "within_cell_summary.json").open("w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
