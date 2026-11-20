# V7 validation, sensitivity, and delivery scripts

This directory contains the reproducible validation used for the V7 manuscript and the local V2.1 data-release candidate. The production configuration is: 1.08° search radius, at most 9 neighbours, at least 3 neighbours, `a=10`, `sigma_b=0.046`, station-error fallback `0.046`, minimum history 30 days, quadratic CDF matching, and `rho=0` in the OI covariance model.

## Current V7 scripts

| Script | Purpose | Principal outputs |
|---|---|---|
| `verify_release.py` | NetCDF structure, continuous time axes, finite-value partitions, update flags, restored-date OI checks, SHA-256 manifest, and optional full metric recomputation | `release_verification.json`, `verify_report.txt` |
| `rebuild_manifest.py` | Hash every release file except the manifest itself | `validation/SHA256SUMS.csv` |
| `build_release_zip.py` | Build the release ZIP and independently verify CRC, membership, and packed-file hashes | Release ZIP and `PACKAGE_QA.json` |
| `within_cell_diagnostics.py` | Standardized paired differences, grid-cluster bootstrap, distance strata, and diagnostic equal-variance inversion | `within_cell_pairs.csv`, `within_cell_distance.csv`, `within_cell_summary.json` |
| `hl_background_correlation.py` | Estimate empirical H–L bins and all fitted curves from one source | `hl_background_correlation.csv`, `hl_background_correlation_summary.json` |
| `aerr_reliability.py` | Strict-CV cell-level posterior-error reliability and decile diagnostics | `aerr_cell_metrics.csv`, `aerr_deciles.csv`, `aerr_reliability_summary.json` |
| `s15_sensitivity.py` | Fold-specific CDFs, dynamic daily errors, `dynamic_oi_day`, and common-mask scoring for BASE plus 22 variants | Variant checkpoints, `sensitivity_common_metrics.csv`, `sensitivity_summary.json` |
| `make_v7_main_figures.py` | Main Figs. 3–6 and 9 from authoritative arrays and release tables | PDF/PNG/SVG/TIFF figures and alignment manifests |
| `make_v7_figures.py` | Figs. 7–8 and S1–S5 from authoritative release tables and diagnostics | PDF/PNG/SVG/TIFF figures and alignment manifests |
| `prepare_v7_word.py` | Prepare citation-keyed Word sources and insert the 9 main and 5 SI figures | `manuscript_word.md`, `supplement_word.md` |
| `build_v7_documents.py` | Export DOCX with `md-to-word`, keep each figure with its own caption, and render PDF through Word on Windows | Final DOCX/PDF and build logs |
| `audit_v7_documents.py` | Audit citations, DOCX table/figure counts, unresolved placeholders, and PDF structure | `WORD_PDF_QA.json` |

The quantitative scripts use released data and authoritative validation tables rather than silently substituting an old production run. Set `SM_DATA_DIR` to the extracted candidate release's `data` directory and `SM_WORK_DIR` to a separate writable analysis directory when running the full Python pipeline. The V7 station NetCDF is required for station-correlation sensitivity variants.

## Release verification

Quick verification checks the release structure and every manifest entry, including files not listed in the manifest:

```powershell
python code/validation/verify_release.py --release-root /path/to/datarelease-scienceDB-v2.1
```

For a staged rebuild, full numerical verification can be run before final manifest generation:

```powershell
python code/validation/verify_release.py --release-root /path/to/datarelease-scienceDB-v2.1 --work-dir /path/to/writable-work --full --full-output /path/to/full-metrics --skip-checksum
python code/validation/rebuild_manifest.py --release-root /path/to/datarelease-scienceDB-v2.1
python code/validation/verify_release.py --release-root /path/to/datarelease-scienceDB-v2.1
```

`--skip-checksum` is an explicit staged-build option, not a final SHA-256 pass. The final quick run must be repeated after manifest regeneration. After writing its JSON and text reports, the verifier refreshes only those two report rows in the manifest, preventing self-generated reports from invalidating their own hashes.

Create and read back the distributable ZIP after the final quick verification:

```powershell
python code/validation/build_release_zip.py --release-root /path/to/datarelease-scienceDB-v2.1 --output /path/to/china_sm_36km_v2.1_20261121.zip --audit-output /path/to/PACKAGE_QA.json
```

The package builder rechecks every source hash, requires saved quick and full verification reports to pass, tests ZIP CRCs and membership, and recomputes every listed file hash from the ZIP itself.

## Within-cell station diagnostics

The V7 diagnostic does not infer a physical correlation from a single pooled RMS value under an equal-variance assumption. For each station pair it calculates `RMS(difference) / sqrt(sigma_i^2 + sigma_j^2)` using station-specific TCA errors, then bootstraps grid cells rather than treating pairs as independent.

Across 2,906 pairs in 641 grid cells, the median standardized RMS is 1.236959 (grid-cluster bootstrap 95% interval 1.195952–1.266620, 2,000 replicates). The nearest 500 pairs give 1.213530. Thus observed differences exceed the independent-TCA expectation even at short separation.

Station correlation and point-scale representativeness differences are not separately identifiable from these data. The diagnostic equal-variance inversion yields 38.37% of estimates outside [−1,1]; it is therefore not used for production or interpreted as a physical estimate. The earlier claim that a pooled inversion near zero demonstrates independence is withdrawn.

## H–L and posterior-error diagnostics

`hl_background_correlation.py` writes the empirical points and fitted curves from the same estimator, preventing the former mismatch between plotted points and the SI fit. The V7 three-parameter fit uses 4,116 days and 1,932 station-bearing cells: amplitude 0.406426, e-folding length 0.886579°, floor 0.067340 and residual RMSE 0.010086. A separate pair-weighted unit-amplitude fit within 1.08° gives `a=7.859919`; it is diagnostic and does not replace the declared production value.

`aerr_reliability.py` reads the per-held-out-day posterior errors saved by `s9_cv.py`. Across 1,826 valid cells, median `Aerr` has Spearman correlation 0.341 with realised cell MAE. The only decile reversal is 4→5, the highest-to-lowest decile MAE increase is 66.47%, and the median cell-wise MAE/Aerr ratio is 1.563. These results do not calibrate `Aerr` as total uncertainty.

## Formal V7 sensitivity

`s15_sensitivity.py` reuses the production Python implementation, including fold-specific CDF matching, dynamic daily observation errors and `dynamic_oi_day`. It evaluates BASE and 22 variants over all 4,116 common analysis days. Scoring uses the intersection of finite observation/background/analysis masks across every listed variant so that changing sample availability is not mistaken for a parameter effect.

The final numerical results are the generated `sensitivity_common_metrics.csv` and `sensitivity_summary.json`; manuscript Table 6 is derived from those files. K03 has no eligible OI target in strict CV. Its common CSV background and fusion values are equal because of background fallback and must not be described as a validated 0% OI improvement. Sensitivity results alone do not justify replacing the documented production configuration.

## Legacy MATLAB scripts

The `sens_*`, `t_cv_official3.m`, `t_aerr.m`, and related MATLAB scripts are retained for historical traceability. Their three-day stepping, static-error, full-sample-CDF, or zero-filled-neighbour conventions do not define the V7 sensitivity or posterior-error results. `t_within_cell_v2.m` and `make_si_fig3.m` are superseded for the within-cell interpretation by `within_cell_diagnostics.py`.

Use the current Python scripts and released common-mask CSV/JSON files when checking V7 manuscript claims. Do not mix legacy outputs with Table 6 or SI S7.4 statistics. The MATLAB mirror source has been updated for the 4,292-day calendar and 0.046 station-error fallback, but an end-to-end MATLAB numerical rerun is not claimed for this candidate.

## Figures and documents

`make_v7_main_figures.py` and `make_v7_figures.py` write the 12 regenerated quantitative figures. The static station map and editable Fig. 2 workflow are supplied separately. Figure QA includes source preflight, render-time panel alignment, minimum text size, collision checks and rendered-image review. QA overlay files are diagnostic only and are not submission figures.

`prepare_v7_word.py` reads the authoritative `paper/delivery_v7_final/source/manuscript.md` and `supplement.md`. It must not be run in reverse: Word-specific Markdown files are generated derivatives, not the scientific source of truth. On Windows, `build_v7_documents.py` runs `md-to-word` with `source/format.txt`, the verified 30-entry DOI-backed BibTeX file and the Copernicus CSL, then renders DOCX to PDF through Word while keeping every figure with its own caption. `audit_v7_documents.py` checks the generated DOCX/PDF deliverables after actual PDF rendering.
