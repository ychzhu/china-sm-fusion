# Download GLDAS-2.1 Noah 3-hourly surface soil moisture (SoilMoi0_10cm_inst)
#
# China subset (lat 18.125-54.125N, lon 73.125-135.125E), one NetCDF per
# 3-hourly step: GLDAS_NOAH025_3H.A{YYYYMMDD}.{HHMM}.021.nc4 via OPeNDAP.
#
# Coverage: GLDAS-2.1 final stream currently publishes up to ~2026-07-31
# (requesting later dates returns HTTP 404, recorded in missing.txt).
#
# Authentication (Earthdata Login / GES DISC):
#   Set the environment variable EDL_TOKEN to a valid Earthdata Login bearer
#   token (profile -> "Generate Token" at urs.earthdata.nasa.gov), or store
#   machine credentials in your ~/.netrc under urs.earthdata.nasa.gov.
#   A GES DISC EULA must be accepted once in a browser for the account.
#
# Output layout (override with --outdir):
#   <outdir>/YYYYMMDD/HHMM.nc      (netCDF classic, China soil-moisture subset)
# Usage: python download_gldas.py [--start 20230101] [--end 20261031]
#                                [--outdir ./raw] [--workers 8]
import os, sys, time, datetime, threading, random, argparse
import requests, netrc

HOURS = ['0000', '0300', '0600', '0900', '1200', '1500', '1800', '2100']
COLL = 'C1342986035-GES_DISC'
# lat index 312..456 (18.125-54.125 N), lon index 1012..1260 (73.125-135.125 E)
SUB = 'SoilMoi0_10cm_inst[0:1:0][312:1:456][1012:1:1260]'

def auth_headers():
    """Bearer token from env EDL_TOKEN; falls back to ~/.netrc credentials."""
    tok = os.environ.get('EDL_TOKEN')
    if tok:
        return {'Authorization': 'Bearer ' + tok}
    n = netrc.netrc()
    _, _, pwd = n.authenticators('urs.earthdata.nasa.gov')
    return {'Authorization': 'Basic ' + __import__('base64').b64encode(
        (n.hosts['urs.earthdata.nasa.gov'][0] + ':' + pwd).encode()).decode()}

def daterange(d0, d1):
    d = datetime.date(*d0)
    end = datetime.date(*d1)
    while d <= end:
        yield d
        d += datetime.timedelta(days=1)

def granule_url(yyyymmdd, hhmm):
    d = datetime.datetime.strptime(yyyymmdd, '%Y%m%d').date()
    doy = (d - datetime.date(d.year, 1, 1)).days + 1
    return (f'https://opendap.earthdata.nasa.gov/collections/{COLL}/granules/'
            f'GLDAS_NOAH025_3H.2.1%3AGLDAS_NOAH025_3H.A{yyyymmdd}.{hhmm}.021.nc4.nc?{SUB}')

def fetch(yyyymmdd, hhmm, outdir, headers):
    od = os.path.join(outdir, yyyymmdd)
    os.makedirs(od, exist_ok=True)
    out = os.path.join(od, hhmm + '.nc')
    if os.path.exists(out) and os.path.getsize(out) > 100000:
        return True, None
    url = granule_url(yyyymmdd, hhmm)
    for attempt in range(4):
        try:
            r = requests.get(url, headers=headers, stream=True, timeout=90)
            if r.status_code == 200:
                with open(out, 'wb') as f:
                    for c in r.iter_content(65536):
                        f.write(c)
                r.close()
                if os.path.getsize(out) > 100000:
                    return True, None
                return False, f'{yyyymmdd}{hhmm} small file {os.path.getsize(out)}'
            elif r.status_code in (403, 429, 500, 502, 503, 504):
                r.close()
                time.sleep(5 * (attempt + 1))
                continue
            else:
                r.close()
                return False, f'{yyyymmdd}{hhmm} HTTP {r.status_code}'
        except Exception:
            time.sleep(4 * (attempt + 1))
    return False, f'{yyyymmdd}{hhmm} failed after retries'

class Worker(threading.Thread):
    def __init__(self, jobs, lock, loglock, outdir, headers):
        super().__init__()
        self.jobs, self.lock, self.loglock = jobs, lock, loglock
        self.outdir, self.headers = outdir, headers
        self.ok = self.fail = 0
    def run(self):
        while True:
            with self.lock:
                if not self.jobs:
                    return
                y, h = self.jobs.pop()
            ok, err = fetch(y, h, self.outdir, self.headers)
            with self.loglock:
                if ok:
                    self.ok += 1
                    with open(os.path.join(self.outdir, 'done.txt'), 'a') as f:
                        f.write(y + h + '\n')
                else:
                    self.fail += 1
                    with open(os.path.join(self.outdir, 'missing.txt'), 'a') as f:
                        f.write(err + '\n')
                    print('MISS:', err, flush=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--start', default='20230101')
    ap.add_argument('--end', default='20261031')
    ap.add_argument('--outdir', default=os.path.join(os.path.dirname(__file__), 'raw'))
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    for p in ('done.txt', 'missing.txt'):
        fp = os.path.join(a.outdir, p)
        if not os.path.exists(fp):
            open(fp, 'a').close()
    headers = auth_headers()
    jobs = []
    for d in daterange((int(a.start[:4]), int(a.start[4:6]), int(a.start[6:])),
                       (int(a.end[:4]), int(a.end[4:6]), int(a.end[6:]))):
        y = d.strftime('%Y%m%d')
        for h in HOURS:
            jobs.append((y, h))
    random.shuffle(jobs)
    print('total jobs:', len(jobs), flush=True)
    lock = threading.Lock(); loglock = threading.Lock()
    ws = [Worker(jobs, lock, loglock, a.outdir, headers) for _ in range(max(1, a.workers))]
    t0 = time.time()
    for w in ws: w.start()
    for w in ws: w.join()
    dt = time.time() - t0
    tot_ok = sum(w.ok for w in ws); tot_fail = sum(w.fail for w in ws)
    print(f'DONE ok={tot_ok} fail={tot_fail} elapsed={dt:.0f}s rate={tot_ok/max(dt,1):.2f}/s', flush=True)

if __name__ == '__main__':
    main()
