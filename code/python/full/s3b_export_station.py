# -*- coding: utf-8 -*-
"""Export the complete V2.1 station product and its reproducible QC tables.

This stage repairs the V2.0 provenance break: the raw/cleaned/flag arrays,
station NetCDF, CF time coordinate, QC decisions and TCA table are now written
from the same continuous 4292-day production axis.
"""
from paths_full import *

import datetime as dt
import json
import os
from pathlib import Path
import shutil
import tempfile

import netCDF4
import numpy as np
import pandas as pd


PROVINCE_EN = {
    "<undefined>": "undefined", "北京市": "Beijing", "天津市": "Tianjin",
    "上海市": "Shanghai", "重庆市": "Chongqing", "河北省": "Hebei",
    "山西省": "Shanxi", "辽宁省": "Liaoning", "吉林省": "Jilin",
    "黑龙江省": "Heilongjiang", "江苏省": "Jiangsu", "浙江省": "Zhejiang",
    "安徽省": "Anhui", "福建省": "Fujian", "江西省": "Jiangxi",
    "山东省": "Shandong", "河南省": "Henan", "湖北省": "Hubei",
    "湖南省": "Hunan", "广东省": "Guangdong", "海南省": "Hainan",
    "四川省": "Sichuan", "贵州省": "Guizhou", "云南省": "Yunnan",
    "陕西省": "Shaanxi", "甘肃省": "Gansu", "青海省": "Qinghai",
    "内蒙古自治区": "Inner Mongolia", "广西壮族自治区": "Guangxi",
    "西藏自治区": "Tibet", "宁夏回族自治区": "Ningxia",
    "新疆维吾尔自治区": "Xinjiang",
}


def _daily_stats(dates, counts):
    counts = np.asarray(counts, dtype=int)
    maximum = int(counts.max())
    return {
        "median_reporting_stations_per_day": float(np.median(counts)),
        "maximum_reporting_stations_on_one_day": maximum,
        "maximum_dates": [dates[i] for i in np.flatnonzero(counts == maximum)],
        "total_reporting_station_days": int(counts.sum()),
    }


def main():
    days = station_axis_ymd().astype(str)
    n_days = len(days)
    raw = np.load(os.path.join(WORK, "raw3596.npy"))
    obj_clean = np.load(os.path.join(WORK, "obj_clean3596.npy"))
    final_flag = np.load(os.path.join(WORK, "final_flag3596.npy"))
    if raw.shape != obj_clean.shape or raw.shape != final_flag.shape:
        raise ValueError(f"Station array shape mismatch: {raw.shape}, {obj_clean.shape}, {final_flag.shape}")
    if raw.shape[1] != n_days:
        raise ValueError(f"Expected {n_days} station days, found {raw.shape[1]}")

    meta = pd.read_csv(os.path.join(WORK, "stations3596_meta.csv"))
    province_source = meta["province"].fillna("<undefined>").astype(str)
    english_provinces = province_source.map(PROVINCE_EN)
    if english_provinces.isna().any():
        unknown = sorted(province_source[english_provinces.isna()].unique())
        raise ValueError(f"Unmapped province labels: {unknown}")
    qc = pd.read_csv(os.path.join(WORK, "final_station_year.csv"))
    tca = pd.read_csv(os.path.join(WORK, "tca3596.csv"))
    cleaned = np.where(final_flag == 0, raw, np.nan).astype(np.float32)
    if not np.array_equal(np.isfinite(cleaned), np.isfinite(obj_clean)):
        raise AssertionError("Final QC mask and Stage-1 cleaned matrix disagree")
    if not np.allclose(cleaned[np.isfinite(cleaned)], raw[np.isfinite(cleaned)], rtol=0, atol=0):
        raise AssertionError("QC export changed a retained raw observation")

    out_station = Path(STATION)
    out_station.mkdir(parents=True, exist_ok=True)
    np.save(out_station / "raw3596.npy", raw)
    np.save(out_station / "cleaned3596.npy", cleaned)
    np.save(out_station / "qc_flag3596.npy", final_flag.astype(np.int8))
    qc.to_csv(out_station / "qc_station_year_3596.csv", index=False)
    tca.to_csv(out_station / "station_tca_error_3596.csv", index=False)

    reasons = qc["reasons"].fillna("").astype(str)
    yearly = qc.groupby("year", as_index=False).agg(
        n_raw=("cma_id", "size"),
        n_keep=("decision", lambda x: int((x == "keep").sum())),
        masked_days=("masked_days", "sum"),
    )
    yearly["n_drop"] = yearly["n_raw"] - yearly["n_keep"]
    yearly["drop_pct"] = (100 * yearly["n_drop"] / yearly["n_raw"]).round(2)
    yearly = yearly[["year", "n_raw", "n_keep", "n_drop", "drop_pct", "masked_days"]]
    yearly.to_csv(out_station / "qc_yearly_summary_3596.csv", index=False)

    raw_daily = np.isfinite(raw).sum(axis=0)
    qc_daily = np.isfinite(cleaned).sum(axis=0)
    pd.DataFrame({
        "date": days,
        "raw_reporting_stations": raw_daily,
        "qc_valid_stations": qc_daily,
    }).to_csv(os.path.join(WORK, "station_daily_counts.csv"), index=False)

    annual_network = []
    for year in range(2015, 2027):
        sel = np.array([x.startswith(str(year)) for x in days])
        qrow = yearly.loc[yearly["year"] == year].iloc[0]
        annual_network.append({
            "year": year,
            "source_station_identifiers": int(qrow["n_raw"]),
            "qc_retained_station_identifiers": int(qrow["n_keep"]),
            "stations_with_at_least_one_raw_value": int((np.isfinite(raw[:, sel]).sum(axis=1) > 0).sum()),
            "stations_with_at_least_one_qc_value": int((np.isfinite(cleaned[:, sel]).sum(axis=1) > 0).sum()),
            "raw_daily_median": float(np.median(raw_daily[sel])),
            "qc_daily_median": float(np.median(qc_daily[sel])),
            "raw_daily_max": int(raw_daily[sel].max()),
            "qc_daily_max": int(qc_daily[sel].max()),
        })
    network = {
        "n_station_identifiers": int(raw.shape[0]),
        "n_days": n_days,
        "axis_start": days[0], "axis_end": days[-1],
        "raw": _daily_stats(days.tolist(), raw_daily),
        "qc": _daily_stats(days.tolist(), qc_daily),
        "stations_with_at_least_one_qc_retained_year": int((qc["decision"] == "keep").groupby(qc["cma_id"]).any().sum()),
        "annual": annual_network,
        "restored_dates": {
            x: {
                "raw_reporting_stations": int(raw_daily[np.flatnonzero(days == x)[0]]),
                "qc_valid_stations": int(qc_daily[np.flatnonzero(days == x)[0]]),
            } for x in ("20230101", "20240101")
        },
    }
    with open(os.path.join(WORK, "station_network_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(network, fh, ensure_ascii=False, indent=2)

    release_summary = {
        "version": RELEASE_VERSION, "shape_station_time": list(raw.shape),
        "n_days": n_days, "finite_raw": int(np.isfinite(raw).sum()),
        "finite_cleaned": int(np.isfinite(cleaned).sum()),
        "flag_0_cells": int((final_flag == 0).sum()),
        "flag_1_cells": int((final_flag == 1).sum()),
        "flag_2_cells": int((final_flag == 2).sum()),
        "masked_station_days": int(qc["masked_days"].sum()),
        "station_years": int(len(qc)),
        "kept_station_years": int((qc["decision"] == "keep").sum()),
        "dropped_station_years": int((qc["decision"] != "keep").sum()),
        "station_tca_valid": int(tca["station_sigma"].notna().sum()),
        "station_tca_median": float(tca["station_sigma"].median()),
        "station_tca_mean": float(tca["station_sigma"].mean()),
        "station_fallback": STA_FALLBACK_ERR,
        "smap_background_error": SMAP_TCA_ERR,
        "time_matches_time_ymd": True,
    }
    with open(os.path.join(WORK, "station_release_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(release_summary, fh, ensure_ascii=False, indent=2)

    coverage = meta.copy()
    coverage["province"] = english_provinces
    coverage = coverage.rename(columns={"cma_id": "station_id"})
    coverage["nvalid_days"] = np.isfinite(cleaned).sum(axis=1)
    coverage["missing_rate_pct"] = (100 * (1 - coverage["nvalid_days"] / n_days)).round(2)
    coverage = coverage[["station_id", "lat", "lon", "province", "nvalid_days", "missing_rate_pct"]]
    coverage_dir = Path(DATA) / "coverage"
    coverage_dir.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(coverage_dir / "coverage_station.csv", index=False)

    fd, temp_path = tempfile.mkstemp(suffix=".nc", prefix="china_station_v21_")
    os.close(fd)
    os.unlink(temp_path)
    ids = meta["cma_id"].astype(int).astype(str).to_numpy()
    provinces = english_provinces.to_numpy()
    id_len = max(6, max(len(x) for x in ids))
    pr_len = max(14, max(len(x) for x in provinces))
    out = None
    try:
        out = netCDF4.Dataset(temp_path, "w", format="NETCDF4")
        out.createDimension("station", raw.shape[0])
        out.createDimension("time", n_days)
        out.createDimension("id_len", id_len)
        out.createDimension("pr_len", pr_len)
        v = out.createVariable("station_id", "S1", ("station", "id_len"))
        v.long_name = "CMA station identifier"
        v[:] = netCDF4.stringtochar(ids.astype(f"S{id_len}"))
        v = out.createVariable("province", "S1", ("station", "pr_len"))
        v.long_name = "administrative province in English"
        v[:] = netCDF4.stringtochar(provinces.astype(f"S{pr_len}"))
        v = out.createVariable("lat", "f8", ("station",))
        v.units = "degrees_north"
        v[:] = meta["lat"].to_numpy(float)
        v = out.createVariable("lon", "f8", ("station",))
        v.units = "degrees_east"
        v[:] = meta["lon"].to_numpy(float)
        v = out.createVariable("time", "f8", ("time",))
        v.units = "days since 2015-01-01 00:00:00"
        v.calendar = "standard"
        v[:] = np.arange(n_days, dtype=np.float64)
        v = out.createVariable("time_ymd", "i4", ("time",))
        v.long_name = "calendar date"
        v[:] = days.astype(np.int32)
        v = out.createVariable("sm", "f4", ("time", "station"), zlib=True,
                               complevel=6, shuffle=True, fill_value=np.float32(-9999),
                               chunksizes=(30, raw.shape[0]))
        v.units = "m3 m-3"
        v.long_name = "three-stage QC daily mean volumetric soil moisture at 10 cm depth"
        v[:] = np.ma.masked_invalid(cleaned.T)
        v = out.createVariable("qc_flag", "i1", ("time", "station"), zlib=True,
                               complevel=6, shuffle=True, chunksizes=(30, raw.shape[0]))
        v.flag_values = np.array([0, 1, 2], dtype=np.int8)
        v.flag_meanings = "valid masked_non_informative hard_dropped"
        v[:] = final_flag.T.astype(np.int8)
        out.title = "China automatic soil moisture network (3596 stations), QC daily 10 cm soil moisture"
        out.Conventions = "CF-1.8"
        out.license = "CC-BY-4.0"
        out.version = RELEASE_VERSION
        out.temporal_coverage = "2015-01-01 to 2026-10-01 inclusive; continuous 4292-day calendar"
        out.station_background_error_fallback = STA_FALLBACK_ERR
        out.station_fallback_basis = "three-decimal median of valid station TCA error estimates; equals SMAP fixed background error"
        out.date_created = dt.datetime.now(dt.timezone.utc).isoformat()
        out.history = "V2.1 restores 2023-01-01 and 2024-01-01 from annual source CSVs and rebuilds CF time consistently with time_ymd."
        out.close()
        out = None
        shutil.copyfile(temp_path, STATION3596_NC)
        with open_nc(STATION3596_NC) as check:
            cf_dates = netCDF4.num2date(
                check.variables["time"][:], check.variables["time"].units,
                check.variables["time"].calendar, only_use_cftime_datetimes=False,
            )
            cf_ymd = np.array([int(x.strftime("%Y%m%d")) for x in cf_dates], dtype=np.int32)
            stored_ymd = np.asarray(check.variables["time_ymd"][:], dtype=np.int32)
            if not np.array_equal(cf_ymd, stored_ymd) or not np.array_equal(stored_ymd, days.astype(np.int32)):
                raise AssertionError("CF time and time_ymd do not describe the same continuous date axis")
    finally:
        if out is not None:
            out.close()
        if os.path.exists(temp_path):
            os.unlink(temp_path)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.2, 5.8), constrained_layout=True)
    sc = ax.scatter(coverage["lon"], coverage["lat"], c=coverage["missing_rate_pct"],
                    s=5, cmap="viridis_r", vmin=0, vmax=100, linewidths=0)
    ax.set(xlabel="Longitude (°E)", ylabel="Latitude (°N)", xlim=(73, 136), ylim=(17, 54),
           title="V2.1 station missing-data rate (continuous 2015-01-01 to 2026-10-01)")
    ax.set_aspect("equal", adjustable="box")
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("Missing days (%)")
    fig.savefig(coverage_dir / "coverage_station.png", dpi=300)
    plt.close(fig)

    print(json.dumps({"station_release": release_summary, "network": network}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
