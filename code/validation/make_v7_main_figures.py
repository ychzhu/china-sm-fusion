#!/usr/bin/env python3
"""Regenerate V7 main-manuscript quantitative figures from authoritative outputs."""
from pathlib import Path
import argparse
import json
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

_QA = os.environ.get("NATURE_FIGURE_SCRIPTS")
if _QA:
    sys.path.insert(0, _QA)
try:
    from audit_panel_alignment import require_matplotlib_panel_alignment
except ImportError:
    require_matplotlib_panel_alignment = None

ROOT = Path(__file__).resolve().parents[3]
WORK = Path(os.environ.get("SM_WORK_DIR", ROOT / "analysis" / "v7_work"))
RELEASE = Path(os.environ.get("SM_RELEASE_DIR", ROOT / "datarelease-scienceDB-v2.1"))
OUT = Path(os.environ.get("SM_FIGURE_DIR", ROOT / "paper" / "delivery_v7_final" / "figures"))
GEO_PATH = Path(__file__).resolve().parents[1] / "resources" / "ne_10m_admin_0_countries.geojson"
GEO = json.loads(GEO_PATH.read_text(encoding="utf-8"))

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8.0, "axes.titlesize": 8.5, "axes.labelsize": 8.0,
    "xtick.labelsize": 7.2, "ytick.labelsize": 7.2, "legend.fontsize": 7.2,
    "axes.linewidth": 0.7, "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
    "savefig.facecolor": "white",
})
COLORS = {"Fusion": "#c84c4c", "SMAP_CDF": "#168c88", "SMAP_raw": "#3478a8",
          "ERA5": "#65607f", "GLDAS": "#d1a143"}
LABELS = {"Fusion": "Fusion", "SMAP_CDF": "CDF", "SMAP_raw": "SMAP",
          "ERA5": "ERA5", "GLDAS": "GLDAS"}
SYSTEMS = ["SMAP_raw", "SMAP_CDF", "Fusion", "ERA5", "GLDAS"]


def panel_label(ax, label):
    ax.text(0, 1, label, transform=ax.transAxes, ha="left", va="bottom",
            fontweight="bold", fontsize=9, clip_on=False)


def save(fig, stem, exclude_axes=None, require_labels=True):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    if require_matplotlib_panel_alignment is not None:
        require_matplotlib_panel_alignment(
            fig, json_out=str(OUT / f"{stem}.alignment.json"),
            overlay_svg=str(OUT / f"{stem}.alignment.svg"), tolerance_pt=1.5,
            gutter_tolerance_pt=1.5, require_panel_labels=require_labels,
            exclude_axes=exclude_axes or [], strict=True)
    fig.savefig(OUT / f"{stem}.png", dpi=600, facecolor="white")
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white")
    fig.savefig(OUT / f"{stem}.svg", facecolor="white")
    plt.close(fig)
    print(stem, flush=True)


def draw_boundary(ax):
    for feature in GEO["features"]:
        if feature["properties"].get("ADMIN") not in {"China", "Taiwan"}:
            continue
        geom = feature["geometry"]
        polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
        for poly in polys:
            ring = np.asarray(poly[0], dtype=float)
            ax.plot(ring[:, 0], ring[:, 1], color="#3d3d3d", lw=0.45, zorder=5)
    ax.set_xlim(73, 136); ax.set_ylim(17, 54); ax.set_aspect(1.18)
    ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")


def metrics(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    e = x - y
    return len(x), np.mean(np.abs(e)), np.sqrt(np.mean(e * e)), np.mean(e), np.corrcoef(x, y)[0, 1]


def fig3_maps(g, cdf, fusion):
    dates = ["20230715", "20240115", "20250410"]
    titles = ["SMAP L3 (V09) raw", "CDF-corrected background", "OI fusion"]
    arrays = [g["smap_bg"], cdf, fusion]
    ymd = g["common_ymd"].astype(str)
    fig = plt.figure(figsize=(7.2, 6.4), constrained_layout=False)
    gs = fig.add_gridspec(3, 3, left=0.075, right=0.865, bottom=0.07, top=0.92,
                          wspace=0.34, hspace=0.78)
    axes, image = [], None
    for r, date in enumerate(dates):
        ti = int(np.where(ymd == date)[0][0])
        for c in range(3):
            ax = fig.add_subplot(gs[r, c]); axes.append(ax)
            values = arrays[c][ti]
            ok = np.isfinite(values)
            image = ax.scatter(g["glon"][ok], g["glat"][ok], c=values[ok], s=1.1,
                               cmap="viridis", vmin=0.05, vmax=0.45, linewidths=0,
                               rasterized=True)
            draw_boundary(ax)
            date_label = f"{date[:4]}-{date[4:6]}-{date[6:]}"
            ax.set_title(f"{titles[c]}\n{date_label}" if r == 0 else date_label, pad=8)
            ax.set_xticks([80, 100, 120]); ax.set_yticks([20, 30, 40, 50])
            if r < 2:
                ax.set_xlabel("")
            if c > 0:
                ax.set_ylabel("")
            panel_label(ax, chr(97 + r * 3 + c))
    cax = fig.add_axes([0.895, 0.16, 0.016, 0.68])
    cb = fig.colorbar(image, cax=cax)
    cb.set_label("Soil moisture (m³ m⁻³)")
    save(fig, "fig03_fusion_maps", exclude_axes=[cax])


def fig4_density(g, cdf, fusion):
    obs = g["grid_obs"].T
    ids = g["bearing"]
    series = [g["smap_bg"][:, ids], cdf[:, ids], fusion[:, ids]]
    names = ["SMAP raw", "CDF background", "OI fusion"]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.55), sharex=True, sharey=True,
                             constrained_layout=True)
    bins = np.linspace(0, 0.6, 121)
    for i, (ax, pred, name) in enumerate(zip(axes, series, names)):
        x, y = pred.ravel(), obs.ravel()
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        h, xe, ye = np.histogram2d(y, x, bins=[bins, bins])
        ax.pcolormesh(xe, ye, h.T + 1, cmap="viridis",
                      norm=LogNorm(vmin=1, vmax=max(10, h.max())), rasterized=True)
        ax.plot([0, 0.6], [0, 0.6], color="white", lw=0.7, ls="--")
        n, mae, rmse, bias, corr = metrics(x, y)
        ax.text(0.03, 0.97, f"N = {n:,}\nMAE = {mae:.4f}\nRMSE = {rmse:.4f}\nr = {corr:.3f}",
                transform=ax.transAxes, va="top", ha="left", color="white", fontsize=7.0,
                bbox=dict(facecolor="black", alpha=0.45, edgecolor="none", pad=2))
        ax.set_title(name); ax.set_xlabel("In situ observation (m³ m⁻³)")
        if i == 0:
            ax.set_ylabel("Product (m³ m⁻³)")
        ax.set_xlim(0, 0.6); ax.set_ylim(0, 0.6); ax.set_aspect("equal")
        panel_label(ax, chr(97 + i))
    save(fig, "fig04_scatter_density")


def temporal_figure(filename, stem, period_title):
    df = pd.read_csv(RELEASE / "validation" / filename)
    periods = list(dict.fromkeys(df["Period"].tolist()))
    x = np.arange(len(periods))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.45), constrained_layout=False)
    for ax, metric, title, ylabel in ((axes[0], "MAE", f"{period_title} MAE", "MAE (m³ m⁻³)"),
                                      (axes[1], "CORR", f"{period_title} correlation", "Pearson r")):
        for system in SYSTEMS:
            sub = df[df["system"] == system].set_index("Period").loc[periods]
            ax.plot(x, sub[metric], "o-", ms=3.5, lw=1.25, color=COLORS[system], label=LABELS[system])
        ax.set_xticks(x); ax.set_xticklabels(periods)
        ax.set_title(title); ax.set_ylabel(ylabel)
        ax.set_xlabel("Month" if period_title == "Monthly" else "Season")
        ax.grid(axis="y", color="#dddddd", lw=0.5); ax.set_axisbelow(True)
    handles, labels = axes[1].get_legend_handles_labels()
    fig.subplots_adjust(left=0.08, right=0.95, bottom=0.25, top=0.88, wspace=0.28)
    fig.legend(handles, labels, loc="lower center", ncol=5,
               bbox_to_anchor=(0.5, 0.005), fontsize=7.0)
    save(fig, stem, require_labels=False)


def fig9_annual():
    df = pd.read_csv(RELEASE / "validation" / "station_annual_metrics_inside.csv")
    years = df["Year"].to_numpy(dtype=int)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.9), sharex=True, constrained_layout=True)
    axes[0, 0].plot(years, df["Raw stations"], "o-", ms=3, color=COLORS["SMAP_raw"], label="Raw stations")
    axes[0, 0].plot(years, df["QC-passing"], "s-", ms=3, color=COLORS["Fusion"], label="QC-passing station-years")
    axes[0, 0].set_ylabel("Count")
    for ax, metric, ylabel in [(axes[0, 1], "MAE", "MAE (m³ m⁻³)"),
                               (axes[1, 0], "CORR", "Pearson r"),
                               (axes[1, 1], "BIAS", "BIAS (m³ m⁻³)")]:
        ax.plot(years, df[f"Fusion_{metric}"], "o-", ms=3, color=COLORS["Fusion"], label="Fusion")
        if metric != "BIAS":
            ax.plot(years, df[f"ERA5_{metric}"], "s-", ms=3, color=COLORS["ERA5"], label="ERA5-Land")
        ax.set_ylabel(ylabel); ax.grid(axis="y", color="#e2e2e2", lw=0.45)
    axes[1, 1].axhline(0, color="#555555", lw=0.7)
    handles = [Line2D([0], [0], color=COLORS["SMAP_raw"], marker="o", lw=1.1),
               Line2D([0], [0], color=COLORS["Fusion"], marker="s", lw=1.1),
               Line2D([0], [0], color=COLORS["Fusion"], marker="o", lw=1.1),
               Line2D([0], [0], color=COLORS["ERA5"], marker="s", lw=1.1)]
    fig.get_layout_engine().set(rect=[0, 0, 1, 0.92])
    fig.legend(handles, ["Raw stations", "QC-passing station-years", "Fusion", "ERA5-Land"],
               loc="upper center", ncol=4, bbox_to_anchor=(0.5, 0.995), fontsize=6.8)
    for ax in axes.flat:
        ax.set_xticks(years); ax.set_xticklabels(years, rotation=45, ha="right", rotation_mode="anchor")
    axes[1, 0].set_xlabel("Year"); axes[1, 1].set_xlabel("Year")
    for ax, lab in zip(axes.flat, "abcd"):
        panel_label(ax, lab)
    save(fig, "fig09_yearly")


def main():
    global WORK, RELEASE, OUT
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, default=WORK)
    parser.add_argument("--release-root", type=Path, default=RELEASE)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    args = parser.parse_args()
    WORK, RELEASE, OUT = args.work_dir.resolve(), args.release_root.resolve(), args.output_dir.resolve()
    g = np.load(WORK / "grid_prep_full.npz", allow_pickle=True)
    cdf = np.load(WORK / "smap_cdf_full.npy", mmap_mode="r")
    fusion = np.load(WORK / "fusion_full.npy", mmap_mode="r")
    fig3_maps(g, cdf, fusion)
    fig4_density(g, cdf, fusion)
    temporal_figure("reference_station_metrics_monthly.csv", "fig05_monthly", "Monthly")
    temporal_figure("reference_station_metrics_seasonal.csv", "fig06_seasonal", "Seasonal")
    fig9_annual()


if __name__ == "__main__":
    main()
