# download_era5land_swvl1.py
"""Download ERA5-Land volumetric soil water layer 1 (0-7 cm) over China.

Product : ERA5-Land (CDS dataset 'reanalysis-era5-land'), 6-hourly
Variable: volumetric_soil_water_layer_1 (swvl1, layer 1 = 0-7 cm)
Area    : China box [54, 73, 18, 135] (N, W, S, E)
Period  : 2023-01 .. 2026-10 (46 calendar months)
Sampling: 00/06/12/18 UTC; daily mean is computed downstream. 6-hourly keeps
          the CDS request volume (and download time) ~6x smaller than hourly;
          soil moisture is a slow variable so the daily-mean error is negligible.
Output  : one NetCDF per calendar month: era5land_swvl1_YYYYMM.nc

The CDS service enforces per-request cost limits, so requests are split by
month. Up to MAX_WORKERS requests run in parallel; failures retry once.

Companion note: no separate extension script is needed; extend PERIOD above
for later months.

Run:  python download_era5land_swvl1.py
"""
import cdsapi
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

OUT = os.environ.get("ERA5_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "era5land"))
AREA = [54, 73, 18, 135]  # north, west, south, east
MAX_WORKERS = 3
TIME_STEPS = ["00:00", "06:00", "12:00", "18:00"]

# Base months covered by the released dataset: 2023-01 .. 2026-10. The CDS
# archive keeps publishing up to ~5 days before today; extend the last range
# when you update the product.
PERIOD = [(y, m) for y in range(2023, 2026) for m in range(1, 13)] + \
         [(2026, m) for m in range(1, 11)]


def build_req(year, month):
    return {
        "variable": "volumetric_soil_water_layer_1",
        "year": f"{year:04d}",
        "month": f"{month:02d}",
        "day": [f"{d:02d}" for d in range(1, 32)],
        "time": TIME_STEPS,
        "area": AREA,
        "format": "netcdf",
    }


def fetch_one(year, month):
    name = f"era5land_swvl1_{year}{month:02d}.nc"
    target = os.path.join(OUT, name)
    if os.path.exists(target) and os.path.getsize(target) > 1e6:
        return f"skip (exists) {name}"
    c = cdsapi.Client()
    for attempt in (1, 2):
        try:
            c.retrieve("reanalysis-era5-land", build_req(year, month), target)
            sz = os.path.getsize(target) / 1e6
            return f"DONE {name} ({sz:.1f} MB)"
        except Exception as e:
            if attempt == 2:
                return f"FAIL {name}: {e}"
            time.sleep(30)
    return f"FAIL {name}"


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] tasks={len(PERIOD)} workers={MAX_WORKERS}", flush=True)
    done = 0
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = {ex.submit(fetch_one, y, m): (y, m) for y, m in PERIOD}
        for fut in as_completed(futures):
            done += 1
            print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {done}/{len(PERIOD)} {fut.result()}", flush=True)
    print("ALL DOWNLOADS FINISHED", flush=True)


if __name__ == "__main__":
    main()
