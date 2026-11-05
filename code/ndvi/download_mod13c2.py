"""Download MOD13C2 (Collection 6.1, monthly CMG 0.05 deg) NDVI via CMR.

The LP DAAC cloud data pool redirects through NASA Earthdata Login; requests
drops the Authorization header on cross-host redirects, so the chain is
followed manually, re-attaching basic auth at every step.
"""
import os
import time
import requests

OUT = r"F:\E盘备份\sm_data_products_autoupdate_project\scripts\casmos_chain_v3\ndvi\raw"
CMR = ("https://cmr.earthdata.nasa.gov/search/granules.json"
       "?short_name=MOD13C2&version=061"
       "&temporal=2023-01-01T00:00:00Z,2026-12-31T23:59:59Z"
       "&page_size=300")

parts = open(os.path.expanduser("~/.netrc")).read().split()
AUTH = (parts[parts.index("login") + 1], parts[parts.index("password") + 1])


def collect_urls():
    r = requests.get(CMR, timeout=60)
    r.raise_for_status()
    out = []
    for e in r.json().get("feed", {}).get("entry", []):
        for link in e.get("links", []):
            if link.get("rel") == "http://esipfed.org/ns/fedsearch/1.1/data#" \
                    and link.get("href", "").endswith(".hdf"):
                out.append((e.get("producer_granule_id"), link["href"]))
    return out


def authenticated_get(sess, url, stream=False):
    """Follow the Earthdata redirect chain, re-attaching auth at each hop."""
    for _ in range(8):
        r = sess.get(url, auth=AUTH, allow_redirects=False, timeout=300, stream=stream)
        if r.status_code in (301, 302, 303, 307, 308) and r.headers.get("Location"):
            loc = r.headers["Location"]
            if loc.startswith("/"):
                loc = "https://data.lpdaac.earthdatacloud.nasa.gov" + loc
            r.close()
            url = loc
            continue
        if r.status_code == 200:
            return r
        r.close()
        raise RuntimeError(f"HTTP {r.status_code}")
    raise RuntimeError("too many redirects")


def main():
    os.makedirs(OUT, exist_ok=True)
    sess = requests.Session()
    for gid, url in collect_urls():
        target = os.path.join(OUT, gid + ".hdf")
        if os.path.exists(target) and os.path.getsize(target) > 1000:
            print(f"[have] {gid}", flush=True)
            continue
        for attempt in range(3):
            try:
                r = authenticated_get(sess, url, stream=True)
                with open(target, "wb") as f:
                    for chunk in r.iter_content(1 << 18):
                        f.write(chunk)
                r.close()
                print(f"[ok] {gid} {os.path.getsize(target)} bytes", flush=True)
                break
            except Exception as e:
                print(f"[retry {attempt}] {gid}: {e}", flush=True)
                time.sleep(3)
        time.sleep(0.4)
    print("done", flush=True)


if __name__ == "__main__":
    main()
