#!/usr/bin/env python3
"""Reproducible H-L diagnostic on first-differenced cell innovations.

This preserves the estimator documented in SI S7.7: remove each cell's mean
first-differenced innovation, normalise by its available temporal variability,
and aggregate pair covariances/correlations by great-circle distance.
"""
from pathlib import Path
import argparse
import csv
import json
import os

import numpy as np
from scipy.optimize import curve_fit, minimize_scalar


def haversine(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = np.radians(lat2 - lat1), np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(a)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path,
                        default=Path(os.environ.get("SM_WORK_DIR", "")))
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not (work / "grid_prep_full.npz").is_file():
        raise FileNotFoundError("Set SM_WORK_DIR to the V7 work directory")
    out = (args.output_dir or work).resolve()
    out.mkdir(parents=True, exist_ok=True)

    with np.load(work / "grid_prep_full.npz", allow_pickle=True) as g:
        obs = g["grid_obs"].astype(float)
        bearing = g["bearing"]
        smap = g["smap_bg"][:, bearing].T.astype(float)
        lat = g["glat"][bearing].astype(float)
        lon = g["glon"][bearing].astype(float)
    n = len(bearing)
    x = np.diff(obs - smap, axis=1)
    valid = np.isfinite(x)
    count = valid.sum(axis=1)
    row_mean = np.divide(np.nansum(x, axis=1), count,
                         out=np.zeros(n), where=count > 0)
    xc = np.where(valid, x - row_mean[:, None], 0.0)
    sd = np.sqrt(np.sum(xc ** 2, axis=1))
    corr = (xc @ xc.T) / np.maximum(np.outer(sd, sd), 1e-12)
    common = valid.astype(float) @ valid.astype(float).T
    cov = (xc @ xc.T) / np.maximum(common, 1.0)
    dist = haversine(lat[:, None], lon[:, None], lat[None, :], lon[None, :])
    iu = np.triu_indices(n, 1)
    dd, cc, vv = dist[iu], corr[iu], cov[iu]

    bins = np.arange(0.2, 2.5, 0.1)
    idx = np.digitize(dd, bins)
    records = []
    for b in range(1, len(bins)):
        mask = idx == b
        if mask.sum() < 20:
            continue
        centre = 0.5 * (bins[b - 1] + bins[b])
        records.append({
            "dist_deg": float(centre), "n_pairs": int(mask.sum()),
            "emp_cov": float(np.mean(vv[mask])),
            "emp_corr": float(np.mean(cc[mask])),
            "current_a_model": float(10.0 ** (-centre)),
        })
    fields = ["dist_deg", "n_pairs", "emp_cov", "emp_corr", "current_a_model"]
    with (out / "hl_background_correlation.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)

    d = np.array([r["dist_deg"] for r in records])
    c = np.array([r["emp_corr"] for r in records])
    weights = np.array([r["n_pairs"] for r in records], dtype=float)

    def three_parameter(x, amplitude, length, floor):
        return amplitude * np.exp(-x / length) + floor

    popt, _ = curve_fit(three_parameter, d, c, p0=[0.35, 0.7, 0.08], maxfev=20000)
    amplitude, length, floor = map(float, popt)
    fit_rmse = float(np.sqrt(np.mean((c - three_parameter(d, *popt)) ** 2)))
    use = d <= 1.08

    def weighted_rmse(a):
        return float(np.sqrt(np.sum(weights[use] * (a ** (-d[use]) - c[use]) ** 2) /
                             np.sum(weights[use])))

    a_within = float(minimize_scalar(weighted_rmse, bounds=(2.0, 20.0), method="bounded",
                                     options={"xatol": 1e-9}).x)
    summary = {
        "method": "first-differenced, per-cell-baseline-removed station-minus-raw-SMAP innovations",
        "n_station_bearing_cells": n, "n_analysis_days": int(obs.shape[1]),
        "bin_width_degrees": 0.1,
        "pair_count_definition": "all unique spatial cell pairs in each distance bin",
        "three_parameter_fit": {
            "formula": "A*exp(-d/L)+floor", "amplitude": amplitude,
            "e_folding_length_degrees": length, "floor": floor,
            "equivalent_a_for_exponential_component": float(np.exp(1.0 / length)),
            "unweighted_bin_rmse": fit_rmse,
        },
        "within_OI_radius_unit_amplitude_fit": {
            "radius_degrees": 1.08, "n_bins": int(use.sum()), "a": a_within,
            "e_folding_length_degrees": float(1.0 / np.log(a_within)),
            "pair_weighted_rmse": weighted_rmse(a_within),
        },
        "caveat": "Innovations also contain station and representativeness errors; this diagnostic does not uniquely identify the SMAP background-error covariance.",
    }
    with (out / "hl_background_correlation_summary.json").open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
