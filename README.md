# China 36 km Daily Soil Moisture Fusion — Code

Reproducible code for the **China 36 km daily surface soil-moisture (10 cm)
fusion dataset** described in the accompanying *Earth System Science Data*
paper. It combines NASA SMAP L3 SPL3SMP V009 retrievals with the China
Meteorological Administration (CMA) in-situ network over mainland China for
2015–2026, using a three-stage quality control and a reference-free Triple
Collocation error calibration. Version 2 uses day-specific observation errors,
standard OI posterior uncertainty, and station-grouped spatial cross-validation.

This repository contains **code and small static resources only**. The large
data products are downloaded separately (see [Data availability](#data-availability)).

## Quick start

```bash
# 1. clone and enter the repository
git clone https://github.com/ychzhu/china-sm-fusion.git
cd china-sm-fusion

# 2. install dependencies (Python 3.9+; a virtual environment is recommended)
pip install -r requirements.txt

# 3. download the released data package into ./data (see data/README.md)
#    expected: data/station, data/smap_v09, data/fusion, data/reference, ...

# 4. run the authoritative full-domain (3596) pipeline
python code/python/full/run_pipeline_full.py
```

All paths are resolved relative to the scripts and NetCDFs are read through
in-memory buffers, so the pipeline runs from any working directory and on
non-ASCII paths.

## Repository layout

```
.
├── code/
│   ├── config.m                 central MATLAB configuration (paths, parameters)
│   ├── run_all.m                MATLAB master entry (preprocess -> oi_fusion)
│   ├── core/                    canonical MATLAB processing stages
│   │   ├── preprocess.m         station QC -> TCA aggregation -> SMAP CDF
│   │   └── oi_fusion.m          optimal interpolation (analysis + error)
│   ├── lib/                     helper functions (station QC, CDF, logging)
│   ├── validation/              scripts behind every number in the paper
│   │   ├── tca_station_error.py Triple Collocation error calibration
│   │   ├── qc_seasonal.py       seasonal active-window QC (full 3596 domain)
│   │   ├── cv_10fold.m/.py      10-fold spatial cross-validation
│   │   └── hl_background_correlation.m  Hollingsworth–Lonnberg correlation
│   ├── reference/               independent ERA5 / GLDAS validation + downloaders
│   ├── ndvi/                    MODIS NDVI download and stratification
│   ├── resources/               small (<1 MB) geometry/error files (.mat)
│   ├── python/
│   │   ├── full/                authoritative V2 full-domain pipeline (s0–s14)
│   │   ├── p1–p4 + run_pipeline.py   dependency-free Python port (2909 domain)
│   │   ├── eval_direct.py / cv_10fold.py / compare_refs.py
│   │   └── work, work_full/     regenerated intermediates (git-ignored)
│   └── verify_release.py        self-contained read-back verifier
├── data/                        downloaded data (git-ignored; see data/README.md)
├── requirements.txt
├── CITATION.cff
└── LICENSE
```

## Methodology overview

1. **Three-stage quality control**, applied independently per station-year so a
   later sensor drift cannot remove good early data:
   - **Stage 1 — objective hard-fault screening** (physical range, flatline/stuck,
     jumps, active-window completeness, spatial consistency; no error estimate).
   - **Stage 2 — Triple Collocation** on the cleaned station + SMAP + ERA5 triplet
     (minimum 100 common days, positive-variance gate), giving per-station error
     standard deviations without treating any product as truth.
   - **Stage 3 — sigma-driven decisions**: extreme errors are down-weighted and
     only hard-dropped with independent objective evidence, avoiding circularity.
2. **Grid aggregation** of stations to the 36 km SMAP EASE-Grid 2.0 cells
   (0.18° search radius), using the stations actually valid on each day and
   propagating their errors under the empirically verified rho = 0 assumption.
3. **CDF bias correction** of the SMAP background where at least three
   neighbouring grids are available; each cross-validation fold rebuilds its
   transfer functions from training-grid observations only.
4. **Optimal interpolation** with an exponential background-error correlation
   (`a=10`, `R=1.08°`, at most nine neighbours), daily observation-error weights,
   and the standard OI posterior error `Aerr=sqrt(max(A0²-B²/A0²,0))`.

## Reproducing the paper results

- **Full-domain V2 pipeline (recommended):** `python code/python/full/run_pipeline_full.py`
  runs stages s0–s14, including TCA sensitivity, direct inside/outside-grid validation,
  leakage-controlled station-grouped 10-fold CV, native and common-mask five-system
  comparisons, stratified metrics, NetCDF export, and release packaging.
- **Path overrides:** set `SM_DATA_DIR` to the extracted archive root,
  `SM_WORK_DIR` to a writable intermediate directory, and `SM_RELEASE_DIR` to the
  release-package root when they differ from the repository defaults.
- **Read-back verification:** stage s14 independently reads the exported NetCDF and
  writes `validation/verify_report.txt` plus `validation/SHA256SUMS.csv`; it exits
  nonzero if the independent checks fail.
- **MATLAB:** the canonical product code (R2022b + Mapping Toolbox) is in
  `code/core/`; run `code/run_all.m`. `config.m` supports both the operational
  project layout and this standalone repository (auxiliary files in
  `code/resources/`).

## Data availability

The released NetCDF/CSV data (~596 MB) are archived separately and cited by DOI;
they are intentionally not committed here. Download instructions and the expected
folder layout are in [`data/README.md`](data/README.md). The raw CMA station data
are restricted; the quality-controlled station products used in the paper are
included in the data archive.

## Citation

Please cite the ESSD paper and the data archive when using this code. See
[`CITATION.cff`](CITATION.cff) (the ESSD paper DOI will be added after assignment).

## License

CC-BY-4.0 — see [`LICENSE`](LICENSE). No operational credentials (FTP, mail,
CDS/Earthdata tokens) are stored in this repository.
