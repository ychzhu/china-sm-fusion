# China 36 km Daily Surface Soil Moisture Fusion — Code

Reproducible code for the **China 36 km daily surface soil moisture (0–10 cm) fusion
dataset** described in the companion *Earth System Science Data* paper. The pipeline
combines NASA SMAP L3 (SPL3SMP V009) retrievals with the China Meteorological
Administration (CMA) in-situ network over mainland China for 2015–2026, using a
three-stage quality control, a reference-free triple-collocation error calibration,
per-cell CDF bias correction and optimal interpolation.

The dataset itself is archived at **https://doi.org/10.57760/sciencedb.018rl** (CC BY 4.0).

## Quick start

```bash
git clone https://github.com/ychzhu/china-sm-fusion.git
cd china-sm-fusion
pip install -r requirements.txt
# place the extracted data archive so that data/station, data/smap_v09, data/fusion exist
python code/python/full/run_pipeline_full.py
```

All paths are resolved relative to the scripts, and NetCDF files are read through
in-memory buffers, so the pipeline runs from any working directory and survives
non-ASCII paths.

## Repository layout

```
code/
├── python/full/          authoritative full-domain (3,596-station) pipeline, stages s0–s10
├── matlab/full/          stage-by-stage MATLAB port producing the same products
├── validation/           parameter sensitivity, fold-specific CV and error-field reliability
├── lib/                  rank-based CDF matching shared by the MATLAB chain
├── reference/download/   downloaders and preprocessing for ERA5 and GLDAS
├── ndvi/                 MOD13C2 download and stratification
└── resources/            small geometry files (<1 MB)
data/                     extracted release (git-ignored; see data/README.md)
```

## Pipeline stages

| Stage | What it does |
|---|---|
| s0 | Build the 3,596 × 4,290 raw station matrix from the yearly files |
| s1 | Stage 1 objective hard-fault screening (no error estimate) |
| s2 | Stage 2 triple collocation on the station–SMAP–ERA5 triplet |
| s3 | Stage 3 error-driven decisions (extreme threshold 0.1074) |
| s4 | Station-to-grid aggregation and cell error |
| s5 | Deterministic capped neighbour selection |
| s6 | Rank-based CDF adjustment of the SMAP background |
| s7 | Optimal-interpolation fusion with daily observation error |
| s8 | In-sample pooled evaluation |
| s9 | Ten-fold spatial cross-validation |
| s10 | Five-system reference comparison |

## Reproducing the paper

| Result | Script |
|---|---|
| Table 1, native-grid direct evaluation, and Table 3, annual metrics | `code/python/full/s8_eval_direct.py` |
| Table 2 and Tables S5–S6, strict station-grouped CV | `code/python/full/s9_cv.py` |
| Table 4, seasonal metrics, and Table 5/Fig. 7, neighbour-density strata | `code/python/full/s12_stratified_metrics.py` |
| Table 6 and Table S7, BASE plus 22 sensitivity variants | `code/validation/s15_sensitivity.py` |
| Tables S9 and Fig. S4, posterior error reliability | Released Aerr decile tables and `code/validation/make_v6_figures.py` |
| Fig. S1, coverage | `code/validation/make_v6_figures.py` and `validation/coverage_summary.csv` |
| Fig. S2, within-cell station differences | `code/validation/within_cell_diagnostics.py` and `make_v6_figures.py` |
| Fig. S3, H–L background-error structure | Released H–L tables and `code/validation/make_v6_figures.py` |
| Final release read-back and full headline-metric recomputation | `code/validation/verify_release.py` |

The sensitivity scripts need the intermediate products produced by stages s0–s5 and a separate writable work directory. See `code/validation/README.md` for the exact conventions, common-mask scoring, variant tags, and the status of legacy MATLAB scripts.

## Requirements

- **Python** 3.9+ with numpy, scipy, pandas, netCDF4 and h5py. No MATLAB is required for the production pipeline or V6 validation.
- **MATLAB** R2022b or later is optional for the historical MATLAB port. Legacy MATLAB sensitivity scripts are retained for traceability but do not define the V6 results.

## Data availability

The released NetCDF and CSV products are archived at Science Data Bank under
https://doi.org/10.57760/sciencedb.018rl and are cited by DOI. They are not committed to
this repository. `data/README.md` gives the expected folder layout. The raw CMA station
data are restricted; the quality-controlled station products used in the paper are part
of the archive.

## Notes on the error model

Production OI treats co-located station errors as independent; this is an explicit approximation, not a conclusion established by the pair-difference test. The V6 diagnostic evaluates 2,906 co-located pairs with station-specific TCA errors and a 2,000-replicate grid-cluster bootstrap. The median standardized pair RMS is 1.235182 (95% interval 1.193845–1.264823), and the nearest 500 pairs give 1.213708. Observed pair differences therefore exceed the independent-TCA expectation, including at short separation.

Station-error correlation and point-scale representativeness differences cannot be separately identified from these data. An equal-variance diagnostic inversion produces many physically invalid correlations and is not used in production. The earlier V5 claim that a pooled inversion near zero demonstrates independence is withdrawn. The posterior error field describes interpolation error against the box mean under the stated covariance approximation; comparison with individual stations also contains point-scale representativeness differences. See the companion paper, Sect. 3.3 and Supplement S7.4, and `code/validation/within_cell_diagnostics.py` for the algebra and diagnostics.

## Citation

Cite the data archive DOI and the companion ESSD paper. See `CITATION.cff`.

## Licence

Code: CC BY 4.0 (see `LICENSE`). No operational credentials (FTP, mail, CDS or Earthdata
tokens) are stored in this repository.
