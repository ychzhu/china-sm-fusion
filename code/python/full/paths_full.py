# -*- coding: utf-8 -*-
"""Shared locations for the full-domain (3596) three-stage QC + fusion pipeline.

All paths are resolved relative to this file so the pipeline runs from any
working directory and handles non-ASCII paths. NetCDFs are opened through
in-memory buffers because the HDF5/netCDF C libraries cannot create or open
files on this host's non-ASCII release path directly.

Data directory: by default the released data are expected under the repository
``data/`` folder. Set the environment variable ``SM_DATA_DIR`` to the root of
an extracted data archive (the folder containing ``station/``, ``smap_v09/``,
``fusion/`` ...) to run from a different location.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.dirname(HERE)
CODE_DIR = os.path.dirname(PYDIR)
ROOT = os.path.dirname(CODE_DIR)

# Released data root: SM_DATA_DIR override wins, otherwise repo/data.
_data_env = os.environ.get("SM_DATA_DIR")
DATA = _data_env if _data_env else os.path.join(ROOT, "data")
STATION = os.path.join(DATA, "station")
SMAP = os.path.join(DATA, "smap_v09")
FUSION = os.path.join(DATA, "fusion")
REFERENCE = os.path.join(DATA, "reference")
WORK = os.environ.get("SM_WORK_DIR", os.path.join(PYDIR, "work_full"))
os.makedirs(WORK, exist_ok=True)

# Authoritative released NetCDFs (3596 domain).
STATION3596_NC = os.path.join(
    STATION, "staD10_qc_3596_20150101_20261001.nc")
FUSION3596_NC = os.path.join(
    FUSION, "fusion_sm_36km_3596_20150331_20261001.nc")

# key constants of the full-domain pipeline
# V2.1: both fixed fallbacks use the three-decimal median of their
# respective valid TCA distributions, as approved for this revision.
SMAP_TCA_ERR = 0.046
STA_FALLBACK_ERR = 0.046
SIGMA_EXTREME = 0.1074
GRID_RES = 0.18
RELEASE_VERSION = "2.1"


def open_nc(path):
    """Open a netCDF file in memory, supporting non-ASCII paths."""
    import netCDF4
    with open(path, "rb") as fh:
        buf = fh.read()
    return netCDF4.Dataset(os.path.basename(path), memory=buf)


def station_axis_ymd():
    """Return the complete station calendar as YYYYmmdd strings.

    V2.0 incorrectly treated its already-incomplete station NetCDF as the
    date authority.  The declared calendar is now independent of that file;
    s0_build_raw verifies every annual source column against it.
    """
    import datetime as dt
    import numpy as np

    start = dt.date(2015, 1, 1)
    end = dt.date(2026, 10, 1)
    return np.array([
        (start + dt.timedelta(days=i)).strftime("%Y%m%d")
        for i in range((end - start).days + 1)
    ], dtype="U8")


def grid_geometry():
    """36-km analysis grid coordinates (7313 cells) from the 3596 fusion product."""
    ds = open_nc(FUSION3596_NC)
    glat = ds.variables["lat"][:].astype(float)
    glon = ds.variables["lon"][:].astype(float)
    ds.close()
    return glat, glon
