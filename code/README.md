# Processing code for the China 36 km daily soil-moisture fusion dataset

This folder contains the two equivalent, self-contained implementations of the
full-domain (3,596-station) three-stage quality-control and fusion pipeline that
produced the dataset and every number reported in the accompanying paper:

- a **Python** implementation in `python/full/`, and
- a **MATLAB** implementation in `matlab/full/` that mirrors it stage by stage.

Both chains run the same eleven stages and produce the same products and metrics
(any differences are limited to last-digit floating-point accumulation).

## Repository layout

```
code/
├── python/full/              authoritative Python pipeline (s0 -> s10)
│   ├── paths_full.py         shared locations, constants, in-memory netCDF helper
│   ├── s0_build_raw.py       build the 3596 x 4290 raw matrix from yearly CSVs
│   ├── s1_stage1.py          objective hard-fault screening (no error estimate)
│   ├── s2_stage2_tca.py      Triple Collocation error estimation
│   ├── s3_stage3.py          sigma-driven decisions (extreme threshold 0.1074)
│   ├── s4_grid_agg.py        station -> grid aggregation, rho=0 cell errors
│   ├── s5_neighbors.py       deterministic capped neighbour selection
│   ├── s6_cdf.py             rank-based CDF adjustment of the SMAP background
│   ├── s7_fusion.py          static-weight OI fusion
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
| s8 / m8 in-sample | Fusion MAE / CORR | 0.0323 / 0.8767 |
| s9 / m9 10-fold CV | OI MAE (background 0.0714) | 0.0644 (9.8% improvement) |
| s10 / m10 comparison | Fusion MAE / CORR (best of five systems) | 0.0461 / 0.7917 |
