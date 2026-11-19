# Sensitivity and reliability scripts

These scripts sit on top of the authoritative pipeline in `code/matlab/full` and
reproduce the sensitivity analysis, the cross-validation re-run and the posterior
error reliability test reported in the companion paper.

They were run with **MATLAB R2022b** on the released intermediate products. Every
script expects `paths_full()` on the path and reads `work_full/` from
`code/matlab/full`, so either run them from a copy placed next to that folder or
add both directories to the path.

## Preparatory step

`sens_*` works on the intermediate matrices produced by stages m0–m5 of the
authoritative chain (`raw3596.mat`, `obj_clean3596.mat`, `tca3596.csv`,
`grid_prep_full.mat`, `neighbors_full.mat`). Run those stages once:

```
run_pipeline_full      % or call m0_build_raw ... m5_neighbors individually
```

If you start from the released products instead, `grid_prep_full.mat` must be
rebuilt by `m4_grid_agg` because the sensitivity sweep needs the daily cell
aggregates `gridObs`, not only the exported analysis.

## Files

| Script | Purpose | Output |
|---|---|---|
| `sens_geo_build.m` | Neighbour superset over a ±1.62° box, sorted by distance | `sens_geo_R162K48.mat` |
| `sens_prep.m` | Defines the fixed scoring domain (4,689 cells) | `sens_base_domain.mat` |
| `sens_cdf.m` | Stage A: neighbour selection and CDF correction for one variant | `senscdf_<tag>.mat` |
| `sens_oi.m` | Stage B: static OI weights, daily fusion, in-sample diagnostics | `sensoi_<tag>.mat` |
| `sens_cv.m` | Stage C: strict station-grouped CV for one variant | `senscv_<tag>.mat` |
| `sens_all.m` | Driver: the 21 variants of Table 5 | `sens_summary.csv` |
| `sens_table.m` | Collects the sweep into one CSV | `sens_summary.csv` |
| `cdf_cell.m` | CDF transfer function with a settable polynomial order | — |
| `t_cv_official3.m` | Strict CV with fold-specific CDF reconstruction over all grid cells | `cv3_<mode>.mat` |
| `t_aerr.m` | Reliability of the posterior error against held-out truth | `aerr_reliability.mat` |
| `t_verify.m` | Compares the re-run analysis with the released `fusion_full.mat` | console |
| `t_diag.m` | Coverage provenance: OI, background-only and missing cell-days | console |
| `make_aerr_fig.m` | Supplementary figure of the reliability test | PNG and CSV |
| `make_si_figs.m` | Coverage figure | PNG and CSV |
| `make_si_fig3.m` | Within-cell station-difference figure | PNG and CSV |

## Variant tags

`BASE` (published settings) · `R070` `R090` `R126` `R144` (search box, degrees) ·
`K03` `K06` `K15` `K24` (neighbour cap) · `A03` `A05` `A20` `A30` (background-error
decay) · `SE040` `SE052` `SE060` (background-error standard deviation) ·
`H090` `H180` (minimum history, days) · `P1` `P3` (CDF polynomial order) ·
`M05` (minimum neighbours).

## Reproduction notes

`t_cv_official3.m` reads `validation/cv_fold_assignments.csv` from the released
data archive. Point `SM_DATA_DIR` at the extracted archive and place that one CSV
at `<repo>/datarelease-scienceDB/validation/cv_fold_assignments.csv`, or edit the
path in the script.

The fold-level CDF application counts printed by `t_cv_official3.m` reproduce the
ten values in `validation/cv_station_grouped_summary.json` exactly
(4,560 / 4,537 / 4,519 / 4,467 / 4,515 / 4,547 / 4,364 / 4,544 / 4,558 / 4,526).
The pooled metrics agree to about 1 %, the residual difference coming from the
handling of missing neighbour observations.
