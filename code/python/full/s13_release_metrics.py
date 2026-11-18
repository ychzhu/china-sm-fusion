# -*- coding: utf-8 -*-
"""Collect the authoritative V2 numbers used by the paper and release notes."""
import csv
import json
import os
import numpy as np
from paths_full import WORK


def read_csv(name):
    with open(os.path.join(WORK, name), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def read_json(name):
    with open(os.path.join(WORK, name), encoding="utf-8") as fh:
        return json.load(fh)


def improvement(old, new):
    return 100.0 * (1.0 - new / old)


def by_key(rows, scope, system):
    return next(r for r in rows if r["scope"] == scope and r["system"] == system)


station = read_csv("station_scope_metrics.csv")
joint = read_csv("reference_station_metrics_joint.csv")
direct = read_csv("direct_gridcell_metrics.csv")
cv = read_json("cv_station_grouped_summary.json")
aerr = read_json("aerr_summary.json")
tca = read_json("tca_sensitivity.json")
release = read_json("release_fusion_summary.json")
counts = read_json("station_scope_counts.json")
roles = read_json("reference_roles.json")

samples = np.load(os.path.join(WORK, "cv_results_full.npz"))
b, f, o = samples["allB"], samples["allF"], samples["allO"]
valid = np.isfinite(b) & np.isfinite(f) & np.isfinite(o)
pair_improved = 100.0 * float(np.mean(np.abs(f[valid] - o[valid]) < np.abs(b[valid] - o[valid])))

bg, fu = direct[0], direct[1]
inside_bg, inside_fu = by_key(station, "inside", "SMAP_CDF"), by_key(station, "inside", "Fusion")
result = {
    "release": release,
    "algorithm": {
        "observation_error": "daily valid-station aggregation: sigma_g(t)=sqrt(sum sigma_i^2)/n_g(t)",
        "analysis_error": "standard OI posterior standard deviation: sigma_b*sqrt(1-c^T w)",
        "production": {"a": 10.0, "radius_degrees": 1.08, "max_near": 9},
        "station_geometry": "inside-cell native pairing is primary; outside nearest-grid matching is sensitivity only",
    },
    "direct_gridcell": {
        "background": bg, "fusion": fu,
        "MAE_reduction_percent": improvement(float(bg["MAE"]), float(fu["MAE"])),
        "RMSE_reduction_percent": improvement(float(bg["RMSE"]), float(fu["RMSE"])),
        "CORR_increase": float(fu["CORR"]) - float(bg["CORR"]),
    },
    "station_scopes": station,
    "station_joint_five_system": joint,
    "inside_station_improvement": {
        "MAE_reduction_percent": improvement(float(inside_bg["MAE"]), float(inside_fu["MAE"])),
        "RMSE_reduction_percent": improvement(float(inside_bg["RMSE"]), float(inside_fu["RMSE"])),
    },
    "strict_cv": cv,
    "strict_cv_improvement": {
        "MAE_reduction_percent": improvement(cv["overall"]["background"]["MAE"], cv["overall"]["fusion"]["MAE"]),
        "RMSE_reduction_percent": improvement(cv["overall"]["background"]["RMSE"], cv["overall"]["fusion"]["RMSE"]),
        "pairs_improved_percent": pair_improved,
    },
    "aerr": aerr,
    "tca_sensitivity": tca,
    "station_counts": counts,
    "reference_roles": roles,
}
with open(os.path.join(WORK, "v2_key_metrics.json"), "w", encoding="utf-8") as fh:
    json.dump(result, fh, indent=2, ensure_ascii=False)
print(json.dumps({
    "direct_MAE_reduction": result["direct_gridcell"]["MAE_reduction_percent"],
    "direct_RMSE_reduction": result["direct_gridcell"]["RMSE_reduction_percent"],
    "inside_MAE_reduction": result["inside_station_improvement"]["MAE_reduction_percent"],
    "cv_MAE_reduction": result["strict_cv_improvement"]["MAE_reduction_percent"],
    "cv_RMSE_reduction": result["strict_cv_improvement"]["RMSE_reduction_percent"],
    "cv_pairs_improved": pair_improved,
}, indent=2))
