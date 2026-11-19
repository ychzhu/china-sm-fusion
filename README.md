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
| Main tables 1–4, 6 and the panels of Figs. 3–9 | `code/python/full/run_pipeline_full.py` (s8–s10) |
| Table 5, parameter sensitivity (21 variants) | `code/validation/sens_all.m` |
| Table 6 and Tables S5–S6, strict CV | `code/validation/t_cv_official3.m` |
| Tables S9 and Fig. S4, posterior error reliability | `code/validation/t_aerr.m` |
| Fig. S1, coverage | `code/validation/make_si_figs.m` |
| Fig. S2, within-cell station differences | `code/validation/make_si_fig3.m` |

The sensitivity scripts need the intermediate matrices produced by stages s0–s5; see
`code/validation/README.md` for the exact prerequisites and for the variant tags.

## Requirements

- **Python** 3.9+ with numpy, scipy, pandas, netCDF4 and h5py. No MATLAB is required.
- **MATLAB** R2022b or later for the MATLAB port and the sensitivity scripts. Only base
  functions are used apart from `lib/cdf_adjusting2.m`.

## Data availability

The released NetCDF and CSV products are archived at Science Data Bank under
https://doi.org/10.57760/sciencedb.018rl and are cited by DOI. They are not committed to
this repository. `data/README.md` gives the expected folder layout. The raw CMA station
data are restricted; the quality-controlled station products used in the paper are part
of the archive.

## Notes on the error model

Co-located stations have independent errors. A pair-difference test on 2,906 co-located
pairs gives an inter-station error correlation of -0.02 (bootstrap 5-95 %: -0.07 to
+0.02), so the 1/sqrt(n) aggregation used for the cell observation error is sound.

The same test shows that the point-scale repeatability of a station, about 0.066
m^3 m^-3, is larger than the triple-collocation error of 0.046 m^3 m^-3 that enters the
observation-error covariance. The collocation quantity is the right one when
interpolating against a 36-km box mean, but it is not the uncertainty of an individual
station value, and it does not depend on the separation of the pair, so the difference is
a property of the measurements rather than of sub-grid variability. The consequence is
that the posterior error field describes interpolation error against the box mean;
comparisons of the analysis with single stations will show a larger discrepancy than the
posterior error suggests. See the companion paper, Sect. 3.3 and Supplement S7.4, for the
algebra and the distance-stratified test.

## Citation

Cite the data archive DOI and the companion ESSD paper. See `CITATION.cff`.

## Licence

Code: CC BY 4.0 (see `LICENSE`). No operational credentials (FTP, mail, CDS or Earthdata
tokens) are stored in this repository.
