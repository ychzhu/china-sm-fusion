# Download the published data here

The data files are **not stored in this Git repository** (they total about 596 MB
and GitHub is intended for code). Download the released data package from the
data archive and place its contents in this `data/` folder, keeping the folder
names below.

> **Data archive (SCIDB):** https://doi.org/10.57760/sciencedb.018rl
> The archive is the same package described in the top-level README and is
> licensed CC-BY-4.0.

## Expected layout after downloading

```
data/
├── station/      QC station NetCDF, raw/cleaned/flag arrays, TCA + QC tables
├── smap_v09/     processed SMAP V09 daily / 2-day / error NetCDFs
├── fusion/       full-domain fusion NetCDF
├── reference/    GCS ERA5 and GLDAS-2.1 regridded NetCDFs + comparison tables
├── ndvi/         monthly NDVI NetCDF and per-station NDVI table
└── coverage/     valid-data counts and missing-rate statistics
```

## Verify the download

After copying the files, run the self-contained read-back verifier from the
repository root:

```
python code/verify_release.py
```

It locates the data from its own path and exits 0 only when every check passes.
The verifier needs only numpy and netCDF4.
