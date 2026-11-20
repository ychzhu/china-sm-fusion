# Processing code for the China 36 km daily soil-moisture fusion dataset

This folder contains the Python production implementation and a stage-by-stage MATLAB port of the full-domain (3,596-station) three-stage quality-control and fusion pipeline:

- a **Python** implementation in `python/full/`, and
- a **MATLAB** implementation in `matlab/full/` that mirrors it stage by stage.

The Python implementation is authoritative for the V2.1 release and every number in the revised paper. The MATLAB implementation mirrors the same production stages; its calendar and station-error fallback have been updated for V2.1, but an end-to-end MATLAB rerun is not claimed in this revision.

## Repository layout

```
code/
├── python/full/              authoritative Python pipeline (s0 -> s10)
│   ├── paths_full.py         shared locations, constants, in-memory netCDF helper
│   ├── s0_build_raw.py       build the 3596 x 4292 raw matrix on a continuous calendar from yearly CSVs
│   ├── s1_stage1.py          objective hard-fault screening (no error estimate)
│   ├── s2_stage2_tca.py      Triple Collocation error estimation
│   ├── s3_stage3.py          sigma-driven decisions (extreme threshold 0.1074)
│   ├── s3b_export_station.py export the 4292-day station NetCDF, arrays and QC tables
│   ├── s4_grid_agg.py        station -> grid aggregation, rho=0 cell errors
│   ├── s5_neighbors.py       deterministic capped neighbour selection
│   ├── s6_cdf.py             rank-based CDF adjustment of the SMAP background
│   ├── s7_fusion.py          daily dynamic-error OI and standard posterior Aerr
│   ├── s8_eval_direct.py     in-sample pooled evaluation
│   ├── s9_cv.py              10-fold spatial cross-validation
│   ├── s10_compare_refs.py   five-system reference comparison
│   └── run_pipeline_full.py  run s0 -> s10 in order
├── matlab/full/              faithful MATLAB port of the same stages
│   ├── paths_full.m         shared locations and constants (returns struct P)
│   ├── read_ymd.m, read_tgrid.m   netCDF reading helpers
│   ├── m0_build_raw.m ... m10_compare_refs.m   ordered pipeline stages
│   ├── cv_fold_perm.csv      fixed 10-fold split (matches numpy RandomState 2026)
│   └── run_pipeline_full.m   run m0 -> m10 in order
├── lib/
│   └── cdf_adjusting2.m      rank-based CDF matching (shared by the MATLAB chain)
├── resources/                small geometry/shape .mat files used by the code
├── reference/download/       downloaders for the public reference data
│   ├── download_era5land_swvl1.py   ERA5-Land monthly CDS downloader
│   ├── download_gldas.py            GLDAS-2.1 Noah 3-hourly subset
│   └── build_gldas_daily.py         3-hourly -> daily mean
└── ndvi/
    └── download_mod13c2.py   CMR API downloader for MOD13C2 C6.1
```

Regenerated intermediates are written to `python/work_full/` and
`matlab/full/work_full/`; both are excluded from version control and rebuilt by
the pipeline.

## Requirements

- **Python** 3.9+ with numpy, scipy, pandas, netCDF4 and h5py (see
  `requirements.txt`). No MATLAB is required. Non-ASCII paths are handled through
  in-memory netCDF reads.
- **MATLAB** R2022b or later (Mapping Toolbox optional; the great-circle
  distance is computed with a built-in haversine routine). The MATLAB chain only
  relies on base functions plus `lib/cdf_adjusting2.m`.

## Data location

Both `paths_full` modules resolve the input data relative to the repository. By
default they look for the released data under the repository's data root; set the
environment variable `SM_DATA_DIR` to point at the directory that holds the
`station/`, `smap_v09/`, `reference/` and related folders if the layout differs.

## Running the pipelines

Python, from `code/python/full/`:

```
python run_pipeline_full.py
```

MATLAB, from `code/matlab/full/`:

```
run_pipeline_full
```

## Key reproduced metrics

| Stage | Metric | Value |
|---|---|---|
| s8 / m8 direct grid-cell diagnostic | Fusion MAE / CORR | 0.033851 / 0.865114 (N = 5,003,186; in-sample) |
| s9 / m9 strict station-grouped 10-fold CV | Fusion MAE / CORR | 0.068362 / 0.576726 (background MAE 0.075254; 9.16% MAE reduction; N = 5,006,022) |
| s10 / m10 native in-grid comparison | Fusion MAE / CORR | 0.045562 / 0.798390 (N = 7,232,510) |
