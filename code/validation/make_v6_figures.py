# -*- coding: utf-8 -*-
"""Rebuild the V6 quantitative figures from released machine-readable inputs.

Figure contracts:
- Fig. 7/8: quantitative grids comparing stratified native in-grid skill.
- Fig. S1: image plate + quant, separating eligibility from actual daily coverage.
- Fig. S2: quantitative grid testing the independent within-cell error model.
- Fig. S3: single-panel background-error correlation diagnostic.
- Fig. S4: quantitative grid showing weak, uncalibrated Aerr guidance.
- Fig. S5: quantitative grid on one common five-system mask.

All plotting/export is Python/Matplotlib. No source arrays are modified.
"""
from pathlib import Path
import argparse
import csv
import datetime as dt
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4
import numpy as np

_QA_SCRIPTS = os.environ.get("NATURE_FIGURE_SCRIPTS")
if _QA_SCRIPTS:
    sys.path.insert(0, _QA_SCRIPTS)
from audit_panel_alignment import require_matplotlib_panel_alignment

ROOT = Path(__file__).resolve().parents[3]
RELEASE = ROOT / "datarelease-scienceDB"
VALID = RELEASE / "validation"
FUSION_NC = RELEASE / "data" / "fusion" / "fusion_sm_36km_3596_20150331_20261001.nc"
OUT = ROOT / "paper" / "delivery_v6_final" / "figures"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "sans-serif"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.dpi": 600,
})
BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
VERMILION = "#CC79A7"
PURPLE = "#7A5195"
GREY = "#7A7A7A"
LIGHT = "#D9D9D9"


def rows(name):
    with (VALID / name).open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def open_nc(path):
    # Memory-backed reads work with the non-ASCII project path on Windows.
    with path.open("rb") as fh:
        return netCDF4.Dataset(path.name, memory=fh.read())


def panel_label(ax, label):
    ax.annotate(label, xy=(0, 1), xycoords="axes fraction", xytext=(-13, 7),
                textcoords="offset points", fontsize=8, fontweight="bold",
                ha="left", va="bottom")


def save_figure(fig, stem, exclude_axes=()):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    print("alignment", stem, flush=True)
    require_matplotlib_panel_alignment(
        fig, json_out=str(OUT / (stem + ".alignment.json")),
        overlay_svg=str(OUT / (stem + ".alignment.svg")),
        tolerance_pt=1.5, gutter_tolerance_pt=1.5,
        exclude_axes=list(exclude_axes), strict=True)
    print("export PNG", stem, flush=True)
    fig.savefig(OUT / (stem + ".png"), dpi=600, bbox_inches="tight")
    print("export PDF/SVG", stem, flush=True)
    fig.savefig(OUT / (stem + ".pdf"), bbox_inches="tight")
    fig.savefig(OUT / (stem + ".svg"), bbox_inches="tight")
    print("export TIFF", stem, flush=True)
    fig.savefig(OUT / (stem + ".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def stratified_figure(source, stem, xlabel):
    data = rows(source)
    labels = [r.get("class") or r.get("Density class") or r.get("Vegetation class") for r in data]
    mae = np.array([float(r["MAE"]) for r in data])
    corr = np.array([float(r["CORR"]) for r in data])
    x = np.arange(len(data))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.75), constrained_layout=True)
    for ax, values, ylabel, letter in zip(axes, (mae, corr),
                                          ("MAE (m³ m⁻³)", "Pearson correlation"), "ab"):
        bars = ax.bar(x, values, width=0.66, color=[SKY, BLUE, GREEN, PURPLE][:len(x)],
                      edgecolor="white", linewidth=0.5)
        ax.set_xticks(x, labels)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        panel_label(ax, letter)
        if ylabel.startswith("MAE"):
            ax.set_ylim(0, max(values) * 1.22)
            for bar, value in zip(bars, values):
                ax.text(bar.get_x() + bar.get_width()/2, value + max(values)*0.025,
                        f"{value:.4f}", ha="center", va="bottom", fontsize=6)
        else:
            ax.set_ylim(0, 1.0)
    save_figure(fig, stem)


def coverage_figure():
    with open_nc(FUSION_NC) as ds:
        lat = ds.variables["lat"][:].astype(float)
        lon = ds.variables["lon"][:].astype(float)
        sm = np.ma.filled(ds.variables["sm_oi"][:], np.nan)
        oi = np.ma.filled(ds.variables["do_oi"][:], 0).astype(bool)
        ymd = ds.variables["time_ymd"][:].astype(int)
    finite = np.isfinite(sm)
    finite_frac = finite.mean(axis=0)
    oi_frac = oi.mean(axis=0)
    daily = oi.sum(axis=1)
    eligible = oi.any(axis=0)
    dates = [dt.datetime.strptime(str(v), "%Y%m%d") for v in ymd]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), constrained_layout=True)
    map_specs = [
        (axes[0, 0], finite_frac, "Finite-value fraction", "viridis", (0, 1)),
        (axes[0, 1], oi_frac, "Actual OI-update fraction", "viridis", (0, 1)),
    ]
    cbars = []
    for ax, values, title, cmap, limits in map_specs:
        sc = ax.scatter(lon, lat, c=values, s=1.2, cmap=cmap, vmin=limits[0], vmax=limits[1],
                        marker="s", linewidths=0, rasterized=True)
        ax.set_title(title)
        ax.set_xlabel("Longitude (°E)")
        ax.set_ylabel("Latitude (°N)")
        cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.025)
        cbars.append(cb)
    axes[1, 0].plot(dates, daily, color=BLUE, linewidth=0.65)
    axes[1, 0].set_ylabel("Cells with OI update")
    axes[1, 0].set_xlabel("Date")
    axes[1, 0].grid(color=LIGHT, linewidth=0.5, alpha=0.7)
    axes[1, 0].set_axisbelow(True)
    axes[1, 1].scatter(lon[~eligible], lat[~eligible], s=1.2, c=GREY, marker="s",
                       linewidths=0, rasterized=True, label="Not update-eligible")
    axes[1, 1].scatter(lon[eligible], lat[eligible], s=1.2, c=BLUE, marker="s",
                       linewidths=0, rasterized=True, label="Update-eligible")
    axes[1, 1].set_xlabel("Longitude (°E)")
    axes[1, 1].set_ylabel("Latitude (°N)")
    axes[1, 1].legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.38))
    for ax, letter in zip(axes.flat, "abcd"):
        panel_label(ax, letter)
    save_figure(fig, "figS1_coverage",
                exclude_axes=[cb.ax for cb in cbars])


def within_cell_figure():
    pairs = rows("within_cell_pairs.csv")
    distance = rows("within_cell_distance.csv")
    with (VALID / "within_cell_summary.json").open(encoding="utf-8") as fh:
        summary = json.load(fh)
    standard = np.array([float(r["standardized_rms"]) for r in pairs])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True)
    ax = axes[0]
    ax.hist(standard, bins=np.arange(0, 3.01, 0.075), color=SKY, edgecolor="white", linewidth=0.35)
    ax.axvline(1, color="black", linestyle="--", linewidth=1, label="Independent-model expectation")
    med = summary["median_standardized_rms"]
    ci = summary["cluster_bootstrap"]["median_standardized_rms_ci"]
    ax.axvline(med, color=VERMILION, linewidth=1.2, label=f"Median = {med:.3f}")
    ax.axvspan(ci[0], ci[1], color=VERMILION, alpha=0.12, linewidth=0,
               label="Cell-clustered 95% interval")
    ax.set_xlabel("RMS / √(σᵢ² + σⱼ²)")
    ax.set_ylabel("Station pairs")
    ax.set_xlim(0, 3)
    ax.legend(loc="upper right")
    panel_label(ax, "a")
    ax = axes[1]
    centres = np.array([(float(r["distance_lo_km"]) + float(r["distance_hi_km"]))/2 for r in distance])
    rms = np.array([float(r["median_rms_difference"]) for r in distance])
    expected = np.array([float(r["median_independent_tca_rms"]) for r in distance])
    ns = [int(r["n_pairs"]) for r in distance]
    ax.plot(centres, rms, color=BLUE, marker="o", markersize=4, label="Observed pair RMS")
    ax.plot(centres, expected, color=GREY, marker="s", markersize=3.5,
            linestyle="--", label="Zero-correlation TCA prediction")
    close = summary["closest_500_pairs"]
    ax.scatter([summary["closest_500_median_separation_km"]], [close["median_rms_difference"]],
               marker="D", s=32, facecolor="white", edgecolor=VERMILION, linewidth=1.1,
               label="Closest 500 pairs", zorder=4)
    ax.set_xlabel("Pair separation (km)")
    ax.set_ylabel("RMS difference (m³ m⁻³)")
    ax.set_ylim(0.055, 0.108)
    ax.legend(loc="lower right")
    panel_label(ax, "b")
    save_figure(fig, "figS2_within_cell")


def hl_figure():
    data = rows("hl_background_correlation.csv")
    d = np.array([float(r["dist_deg"]) for r in data])
    corr = np.array([float(r["emp_corr"]) for r in data])
    good = np.isfinite(corr)
    d, corr = d[good], corr[good]
    xx = np.linspace(0.05, 2.45, 400)
    fig, ax = plt.subplots(figsize=(4.4, 3.15), constrained_layout=True)
    ax.scatter(d, corr, s=18, facecolor="white", edgecolor=GREEN, linewidth=0.8,
               label="Empirical first-difference correlation", zorder=4)
    for a in (3, 5, 20, 30):
        ax.plot(xx, a ** (-xx), color=LIGHT, linewidth=0.75, zorder=1)
    ax.plot(xx, 10 ** (-xx), color=BLUE, linewidth=1.35, label="Operational a = 10")
    ax.plot(xx, 0.406*np.exp(-xx/0.8869) + 0.067, color=VERMILION,
            linewidth=1.35, linestyle="--", label="Three-parameter fit")
    ax.axvline(1.08, color="black", linestyle=":", linewidth=1, label="Search radius 1.08°")
    ax.axhline(0, color="black", linewidth=0.6)
    ax.set_xlim(0, 2.45)
    ax.set_ylim(-0.04, 0.66)
    ax.set_xlabel("Great-circle distance (degrees)")
    ax.set_ylabel("Normalised background-error correlation")
    ax.legend(loc="upper right")
    save_figure(fig, "figS3_hl_correlation")


def aerr_figure():
    data = rows("aerr_deciles.csv")
    dec = np.array([int(r["decile"]) for r in data])
    aerr = np.array([float(r["median_Aerr"]) for r in data])
    mae = np.array([float(r["mean_abs_error"]) for r in data])
    rmse = np.array([float(r["rms_error"]) for r in data])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True)
    ax = axes[0]
    ax.plot(dec, aerr, marker="o", color=GREY, label="Median Aerr")
    ax.plot(dec, mae, marker="o", color=BLUE, label="Realised MAE")
    ax.plot(dec, rmse, marker="o", color=VERMILION, label="Realised RMSE")
    ax.set_xticks(dec)
    ax.set_xlabel("Aerr decile")
    ax.set_ylabel("Error (m³ m⁻³)")
    ax.set_ylim(0.03, 0.11)
    ax.legend(loc="upper left")
    panel_label(ax, "a")
    ax = axes[1]
    line = np.linspace(0.03, 0.095, 100)
    ax.plot(line, line, color="black", linestyle="--", linewidth=1, label="1:1")
    ax.plot(line, 1.6*line, color=GREY, linestyle=":", linewidth=1.1, label="1.6:1")
    ax.scatter(aerr, mae, s=28, color=BLUE, edgecolor="white", linewidth=0.5, zorder=4)
    ax.annotate("Decile 7→8 reversal",
                xy=((aerr[6] + aerr[7]) / 2, (mae[6] + mae[7]) / 2),
                xytext=(0.97, 0.45), textcoords="axes fraction",
                ha="right", va="center", fontsize=6,
                arrowprops=dict(arrowstyle="-|>", color=GREY, linewidth=0.7,
                                shrinkA=3, shrinkB=4))
    ax.set_xlim(0.03, 0.095)
    ax.set_ylim(0.03, 0.095)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Median Aerr (m³ m⁻³)")
    ax.set_ylabel("Realised MAE (m³ m⁻³)")
    ax.legend(loc="lower right")
    panel_label(ax, "b")
    save_figure(fig, "figS4_aerr_reliability")


def reference_figure():
    data = [r for r in rows("reference_station_metrics_joint.csv") if r["scope"] == "inside"]
    order = ["SMAP_raw", "SMAP_CDF", "Fusion", "GLDAS", "ERA5"]
    labels = ["Raw\nSMAP", "CDF\nbackground", "Fusion", "GLDAS", "ERA5"]
    by = {r["system"]: r for r in data}
    mae = [float(by[k]["MAE"]) for k in order]
    corr = [float(by[k]["CORR"]) for k in order]
    x = np.arange(len(order))
    colors = [GREY, SKY, BLUE, GREEN, PURPLE]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.85), constrained_layout=True)
    for ax, values, ylabel, letter in zip(axes, (mae, corr),
                                          ("MAE (m³ m⁻³)", "Pearson correlation"), "ab"):
        ax.bar(x, values, color=colors, width=0.68, edgecolor="white", linewidth=0.5)
        ax.set_xticks(x, labels)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color=LIGHT, linewidth=0.5)
        ax.set_axisbelow(True)
        panel_label(ax, letter)
    axes[0].set_ylim(0, 0.12)
    axes[1].set_ylim(0, 0.9)
    fig.suptitle("Common mask: N = 6,973,869 native in-grid pairs", fontsize=7.5)
    save_figure(fig, "figS5_reference")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-coverage", action="store_true")
    args = parser.parse_args()
    stratified_figure("station_density_metrics_inside.csv", "fig07_density",
                      "Neighbouring station-bearing cells")
    stratified_figure("station_ndvi_metrics_inside.csv", "fig08_ndvi",
                      "Growing-season NDVI class")
    if not args.skip_coverage:
        coverage_figure()
    within_cell_figure()
    hl_figure()
    aerr_figure()
    reference_figure()
    print("V6 quantitative figures exported to", OUT)


if __name__ == "__main__":
    main()
