# Build daily-mean GLDAS soil moisture (0-10 cm) from the 3-hourly subsets
# downloaded by download_gldas.py.
#
# Input : <rawdir>/YYYYMMDD/HHMM.nc   (8 steps per day, netCDF classic)
# Output: gldas_daily.mat             (lat, lon, dates, sm, nstep)
#         sm : 145 x 249 x N float32, daily mean volumetric soil moisture
#              [m3/m3] = SoilMoi0_10cm_inst [kg/m2] / 100
#         A day is kept only when >= 4 of the 8 steps are present (GLDAS
#         2.1 final stream has gaps after ~2026-07 due to publication lag).
# Usage: python build_gldas_daily.py --raw <rawdir> --out <file.mat>
import os, sys, glob, argparse
import numpy as np
import netCDF4
from scipy.io import savemat

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw', default=os.path.join(os.path.dirname(__file__), 'raw'))
    ap.add_argument('--out', default=os.path.join(os.path.dirname(__file__), 'gldas_daily.mat'))
    a = ap.parse_args()
    days = sorted(os.listdir(a.raw))
    lat = lon = None
    sm_list = []
    dates = []
    nstep = []
    miss_days = []
    for d in days:
        files = sorted(glob.glob(os.path.join(a.raw, d, '*.nc')))
        if not files:
            miss_days.append(d)
            continue
        acc = []
        for fp in files:
            try:
                f = netCDF4.Dataset(fp)
                v = f['SoilMoi0_10cm_inst'][0] / 100.0  # kg/m2 -> m3/m3
                v[v < -1.0] = np.nan                    # GLDAS fill -9999 -> NaN
                lat = f['lat'][:].astype(np.float32) if lat is None else lat
                lon = f['lon'][:].astype(np.float32) if lon is None else lon
                acc.append(v)
                f.close()
            except Exception as e:
                print('ERR', fp, e)
        ns = len(acc)
        nstep.append(ns)
        if ns >= 4:
            sm_list.append(np.nanmean(np.stack(acc, 0), axis=0).astype(np.float32))
            dates.append(int(d))
        else:
            miss_days.append(f'{d} (nstep={ns})')
    sm = np.stack(sm_list, -1)  # 145 x 249 x N
    dates = np.array(dates, dtype=np.int64)
    nstep = np.array(nstep, dtype=np.int64)
    print('grid:', sm.shape, 'days:', len(dates), 'range:', dates[0], dates[-1])
    print('nstep: min %d max %d mean %.1f' % (nstep.min(), nstep.max(), nstep.mean()))
    print('missing days:', miss_days[:10], '... total', len(miss_days))
    print('SM mean %.4f min %.4f max %.4f' % (np.nanmean(sm), np.nanmin(sm), np.nanmax(sm)))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    savemat(a.out, {'lat': lat, 'lon': lon, 'dates': dates, 'sm': sm, 'nstep': nstep})
    print('saved:', a.out, os.path.getsize(a.out))

if __name__ == '__main__':
    main()
