# V6 validation, sensitivity, and delivery scripts

This directory contains the reproducible validation used for the V6 manuscript and the V2 data-release read-back. The production configuration is unchanged: 1.08° search radius, at most 9 neighbours, at least 3 neighbours, `a=10`, `sigma_b=0.046`, minimum history 30 days, quadratic CDF matching, and `rho=0` in the OI covariance model.

## Current V6 scripts

| Script | Purpose | Principal outputs |
|---|---|---|
| `verify_release.py` | NetCDF structure, finite-value partitions, update flags, fallback dates, SHA-256 manifest, and optional full metric recomputation | `release_verification.json`, `verify_report.txt` |
| `rebuild_manifest.py` | Hash every release file except the manifest itself | `validation/SHA256SUMS.csv` |
| `build_release_zip.py` | Build the release ZIP and independently verify CRC, membership, and packed-file hashes | Release ZIP and `PACKAGE_QA.json` |
| `within_cell_diagnostics.py` | Standardized paired differences, grid-cluster bootstrap, distance strata, and diagnostic equal-variance inversion | `within_cell_pairs.csv`, `within_cell_distance.csv`, `within_cell_summary.json` |
| `s15_sensitivity.py` | V2 sensitivity using fold-specific CDFs, dynamic daily errors, and `dynamic_oi_day` | Variant checkpoints and common-mask summary tables |
| `make_v6_figures.py` | Seven regenerated quantitative figures from authoritative release tables | PDF/PNG/SVG/TIFF figures and QA inputs |
| `prepare_v6_word.py` | Deterministically prepare citation-keyed Word sources and insert the 9 main and 5 SI figures | `manuscript_word.md`, `supplement_word.md` |
| `audit_v6_documents.py` | Audit citations, DOCX table/figure counts, unresolved placeholders, and PDF structure | `WORD_PDF_QA.json` |

The quantitative scripts use the released data and authoritative validation tables rather than recalculating production fusion without evidence. Set `SM_DATA_DIR` to the extracted release's `data` directory and `SM_WORK_DIR` to a separate writable analysis directory when running the full Python pipeline.

## Release verification

Quick verification checks the release structure and every manifest entry, including files not listed in the manifest:

```bash
python code/validation/verify_release.py --release-root /path/to/datarelease-scienceDB
```

For a staged rebuild, full numerical verification can be run before final manifest generation:

```bash
python code/validation/verify_release.py \
  --release-root /path/to/datarelease-scienceDB \
  --work-dir /path/to/writable-work \
  --full --full-output /path/to/full-metrics \
  --skip-checksum
python code/validation/rebuild_manifest.py --release-root /path/to/datarelease-scienceDB
python code/validation/verify_release.py --release-root /path/to/datarelease-scienceDB
```

`--skip-checksum` is an explicit staged-build option. A report produced with that option is not a final SHA-256 pass. The final quick run must be repeated after manifest regeneration. After writing its JSON and text reports, the verifier refreshes only those two report rows in the manifest, preventing self-generated reports from invalidating their own hashes.

Create and read back the distributable ZIP after the final quick verification:

```bash
python code/validation/build_release_zip.py \
  --release-root /path/to/datarelease-scienceDB \
  --output /path/to/china_sm_36km_v2.0.zip \
  --audit-output /path/to/PACKAGE_QA.json
```

The package builder rechecks every source hash, requires saved quick and full verification reports to be PASS, tests ZIP CRCs and membership, and recomputes every listed file hash from the ZIP itself.

## Within-cell station diagnostics

The V6 diagnostic does not infer a physical correlation from a single pooled RMS value under an equal-variance assumption. For each station pair it calculates

`RMS(difference) / sqrt(sigma_i^2 + sigma_j^2)` using station-specific TCA errors, then bootstraps grid cells rather than treating pairs as independent.

Across 2,906 pairs in 641 grid cells, the median standardized RMS is 1.235182 (grid-cluster bootstrap 95% interval 1.193845–1.264823, 2,000 replicates). The nearest 500 pairs give 1.213708. Thus observed differences exceed the independent-TCA expectation even at short separation.

Station correlation and point-scale representativeness differences are not separately identifiable from these data. The diagnostic equal-variance inversion yields a median implied correlation of −0.581083 and 38.1624% of estimates outside [−1,1]; it is therefore not used for production or interpreted as a physical estimate. The V5 claim that `rho=-0.02` demonstrates independence is withdrawn.

## Formal V2 sensitivity

`s15_sensitivity.py` reuses the production Python V2 implementation, including fold-specific CDF matching, dynamic daily observation errors, and `dynamic_oi_day`. It evaluates BASE and 22 variants over all 4,114 analysis days. Scoring uses a common finite mask so that changes in sample availability are not mistaken for parameter effects.

On the common CV mask, BASE reduces MAE by 9.1967799%. K15 and K24 give 10.30200% and 10.77461%, while RHO020 and RHO050 give 9.09555% and 8.94722%. These diagnostics do not by themselves justify replacing the documented production configuration or regenerating the released product.

K03 has no eligible OI targets in strict CV. Its common CSV background and fusion values are equal because of background fallback; this must not be described as a validated 0% OI improvement.

## Legacy MATLAB scripts

The `sens_*`, `t_cv_official3.m`, `t_aerr.m`, and related MATLAB scripts are retained for historical traceability. Their three-day stepping, static-error, full-sample-CDF, or zero-filled-neighbour conventions do not define the V6 sensitivity results. `t_within_cell_v2.m` and `make_si_fig3.m` are superseded for the within-cell interpretation by `within_cell_diagnostics.py`.

Use the current Python scripts and the released common-mask CSV/JSON files when checking V6 manuscript claims. Do not mix legacy outputs with the V6 Table 6 or SI S7.4 statistics.

## Figures and documents

`make_v6_figures.py` writes final quantitative figures; the editable Fig. 2 workflow is supplied separately as SVG. Figure QA includes source preflight, minimum text size, collision and alignment checks, and rendered-image review.

`prepare_v6_word.py` reads the authoritative `paper/delivery_v6_final/source/manuscript.md` and `supplement.md`. It must not be run in reverse: the Word-specific Markdown files are generated derivatives, not the scientific source of truth. DOCX export uses `md-to-word`, `source/format.txt`, the 19-entry DOI-backed BibTeX file, and the Copernicus CSL. `audit_v6_documents.py` checks the generated DOCX/PDF deliverables after actual PDF rendering.
