#!/usr/bin/env python3
"""Posterior-error reliability using strict station-grouped CV predictions."""
from pathlib import Path
import argparse
import csv
import json
import os

import numpy as np
from scipy.stats import rankdata


def corr(x, y):
    return float(np.corrcoef(x, y)[0, 1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path,
                        default=Path(os.environ.get("SM_WORK_DIR", "")))
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--min-days", type=int, default=30)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not str(work):
        raise ValueError("Set SM_WORK_DIR or --work-dir")
    out = (args.output_dir or work).resolve()
    out.mkdir(parents=True, exist_ok=True)

    with np.load(work / "cv_results_full.npz") as cv:
        required = {"allF", "allO", "allAerr", "sample_grid"}
        if not required.issubset(cv.files):
            raise KeyError("Rerun s9_cv.py: cv_results_full.npz lacks allAerr")
        pred = cv["allF"].astype(float)
        obs = cv["allO"].astype(float)
        aerr = cv["allAerr"].astype(float)
        grid = cv["sample_grid"].astype(int)
    if not (len(pred) == len(obs) == len(aerr) == len(grid)):
        raise ValueError("CV sample arrays have inconsistent lengths")

    order = np.argsort(grid, kind="stable")
    grid, pred, obs, aerr = grid[order], pred[order], obs[order], aerr[order]
    unique, starts = np.unique(grid, return_index=True)
    ends = np.r_[starts[1:], len(grid)]
    cell_rows = []
    for g, lo, hi in zip(unique, starts, ends):
        count = hi - lo
        if count < args.min_days:
            continue
        aa = aerr[lo:hi]
        error = pred[lo:hi] - obs[lo:hi]
        if not np.isfinite(aa).all() or not np.isfinite(error).all():
            continue
        cell_rows.append({
            "grid_index0": int(g), "n_days": count,
            "median_Aerr": float(np.median(aa)),
            "mean_abs_error": float(np.mean(np.abs(error))),
            "rms_error": float(np.sqrt(np.mean(error ** 2))),
            "mean_bias": float(np.mean(error)),
        })
    if not cell_rows:
        raise ValueError("No cells satisfy the minimum-day rule")
    fields = list(cell_rows[0])
    with (out / "aerr_cell_metrics.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(cell_rows)

    a = np.array([r["median_Aerr"] for r in cell_rows])
    e = np.array([r["mean_abs_error"] for r in cell_rows])
    rms = np.array([r["rms_error"] for r in cell_rows])
    bias = np.array([r["mean_bias"] for r in cell_rows])
    ndays = np.array([r["n_days"] for r in cell_rows])
    edges = np.percentile(a, np.arange(0, 101, 10))
    deciles = []
    for k in range(10):
        mask = (a >= edges[k]) & (a <= edges[k + 1] if k == 9 else a < edges[k + 1])
        if mask.sum() < 5:
            continue
        deciles.append({
            "decile": k + 1, "n_cells": int(mask.sum()),
            "median_Aerr": float(np.median(a[mask])),
            "mean_abs_error": float(np.mean(e[mask])),
            "rms_error": float(np.mean(rms[mask])),
            "ratio": float(np.mean(e[mask]) / np.median(a[mask])),
            "mean_bias": float(np.mean(bias[mask])),
            "lo_Aerr": float(edges[k]), "hi_Aerr": float(edges[k + 1]),
        })
    decile_fields = ["decile", "n_cells", "median_Aerr", "mean_abs_error",
                     "rms_error", "ratio", "mean_bias", "lo_Aerr", "hi_Aerr"]
    with (out / "aerr_deciles.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=decile_fields)
        writer.writeheader()
        writer.writerows(deciles)

    decile_mae = np.array([r["mean_abs_error"] for r in deciles])
    reversals = []
    for k in np.flatnonzero(np.diff(decile_mae) < 0):
        reversals.append({
            "from_decile": int(deciles[k]["decile"]),
            "to_decile": int(deciles[k + 1]["decile"]),
            "from_MAE": float(decile_mae[k]),
            "to_MAE": float(decile_mae[k + 1]),
        })
    summary = {
        "source": "strict station-grouped 10-fold CV with fold-specific CDF and dynamic R",
        "min_held_out_days_per_cell": args.min_days,
        "n_cells": len(cell_rows),
        "median_held_out_days_per_cell": float(np.median(ndays)),
        "overall_cell_summary": {
            "median_Aerr": float(np.median(a)),
            "median_cell_mean_abs_error": float(np.median(e)),
            "median_cell_rms_error": float(np.median(rms)),
            "median_cell_MAE_to_Aerr_ratio": float(np.median(e / a)),
            "mean_cell_bias": float(np.mean(bias)),
        },
        "pearson_Aerr_vs_cell_MAE": corr(a, e),
        "pearson_Aerr_vs_cell_RMSE": corr(a, rms),
        "spearman_Aerr_vs_cell_MAE": corr(rankdata(a), rankdata(e)),
        "deciles": deciles,
        "decile_MAE_reversals": reversals,
        "highest_vs_lowest_decile_MAE_increase_percent":
            float((decile_mae[-1] / decile_mae[0] - 1) * 100),
        "interpretation": "Aerr is a model-conditional posterior standard deviation, not calibrated total error. Decile discrimination can coexist with magnitude underestimation.",
    }
    with (out / "aerr_reliability_summary.json").open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    print(json.dumps({k: v for k, v in summary.items() if k != "deciles"}, indent=2))
    print(json.dumps(deciles, indent=2))


if __name__ == "__main__":
    main()
