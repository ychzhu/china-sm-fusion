# -*- coding: utf-8 -*-
"""V6 parameter sensitivity using the formal V2 fold-CDF and dynamic-R kernel.

All variants use the released station-component folds and all 4,114 common days.
Only station-bearing targets are evaluated: non-bearing CDF grids cannot affect
these scores, so omitting them changes runtime, not the validation definition.
Two passes are used: variant predictions are checkpointed, then every variant is
scored again on the intersection of finite observation/background/analysis masks.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python" / "full"))
from paths_full import WORK, open_nc, STATION3596_NC
from cdf_utils import cdf_adjusting2
from oi_utils import build_geometry, dynamic_oi_day

import argparse
import csv
import json
import os
import warnings

import numpy as np


VARIANTS = {
    "BASE": {},
    "R070": {"radius": 0.70}, "R090": {"radius": 0.90},
    "R126": {"radius": 1.26}, "R144": {"radius": 1.44},
    "K03": {"cap": 3}, "K06": {"cap": 6}, "K15": {"cap": 15}, "K24": {"cap": 24},
    "A03": {"a": 3.0}, "A05": {"a": 5.0}, "A20": {"a": 20.0}, "A30": {"a": 30.0},
    "SE040": {"sigma_b": 0.040}, "SE052": {"sigma_b": 0.052}, "SE060": {"sigma_b": 0.060},
    "H090": {"min_hist": 90}, "H180": {"min_hist": 180},
    "P1": {"poly_order": 1}, "P3": {"poly_order": 3}, "M05": {"min_nbr": 5},
    "RHO020": {"rho": 0.20}, "RHO050": {"rho": 0.50},
}
DEFAULT = dict(radius=1.08, cap=9, a=10.0, sigma_b=0.046,
               min_hist=30, min_nbr=3, poly_order=2, rho=0.0)


def parameters(tag):
    p = DEFAULT.copy()
    p.update(VARIANTS[tag])
    p["tag"] = tag
    return p


def cdf_adjust(smap, reference, order):
    if order == 2:
        return cdf_adjusting2(smap, reference)
    s = np.asarray(smap, float).copy()
    t = np.asarray(reference, float).copy()
    s[(s > 0.6) | (s < 0.04)] = np.nan
    t[(t > 0.6) | (t < 0.04)] = np.nan
    sv, tv = s[np.isfinite(s)], t[np.isfinite(t)]
    if len(sv) < order + 1 or len(tv) < order + 1:
        return s
    ps = np.arange(1, len(sv) + 1) / (len(sv) + 1)
    pt = np.arange(1, len(tv) + 1) / (len(tv) + 1)
    # Same linear-extrapolating quantile rule as cdf_utils.interp1_extrap.
    x, y = ps, np.sort(sv)
    sint = np.interp(pt, x, y)
    left, right = pt < x[0], pt > x[-1]
    sint[left] = y[0] + (y[1] - y[0]) / (x[1] - x[0]) * (pt[left] - x[0])
    sint[right] = y[-1] + (y[-1] - y[-2]) / (x[-1] - x[-2]) * (pt[right] - x[-1])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", np.exceptions.RankWarning)
        coeff = np.polyfit(sint, np.sort(tv) - sint, order)
    out = s + np.polyval(coeff, s)
    out[(out > 1) | (out < 0)] = np.nan
    return out


def haversine(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = np.radians(lat2 - lat1), np.radians(lon2 - lon1)
    aa = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(aa)))


def metrics(obs, background, fusion):
    result = {"N": int(obs.size)}
    for name, values in (("background", background), ("fusion", fusion)):
        d = values - obs
        result[name] = {
            "MAE": float(np.mean(np.abs(d))),
            "RMSE": float(np.sqrt(np.mean(d ** 2))),
            "BIAS": float(np.mean(d)),
            "CORR": float(np.corrcoef(values, obs)[0, 1]),
        }
    result["mae_reduction_pct"] = float(100 * (1 - result["fusion"]["MAE"] / result["background"]["MAE"]))
    result["rmse_reduction_pct"] = float(100 * (1 - result["fusion"]["RMSE"] / result["background"]["RMSE"]))
    return result


class Sensitivity:
    def __init__(self, out):
        self.out = Path(out)
        self.out.mkdir(parents=True, exist_ok=True)
        g = np.load(os.path.join(WORK, "grid_prep_full.npz"), allow_pickle=True)
        nb = np.load(os.path.join(WORK, "neighbors_full.npz"), allow_pickle=True)
        self.bearing = g["bearing"].astype(int)
        self.glat, self.glon = g["glat"], g["glon"]
        self.grid_obs = g["grid_obs"].astype(float)
        self.grid_err_base = g["grid_err_daily"].astype(float)
        self.smap_bg = g["smap_bg"].astype(float)
        self.common_ymd = g["common_ymd"].astype(int)
        self.members = g["members"]
        self.smap_err = float(g["smap_err"][0])
        self.n_eff = len(self.bearing)
        self.n_days = self.smap_bg.shape[0]
        self.fold_of = np.full(self.n_eff, -1, dtype=int)
        va_path = Path(WORK) / "cv_fold_assignments.csv"
        if not va_path.is_file():
            raise FileNotFoundError(va_path)
        with va_path.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                self.fold_of[int(row["bearing_position"])] = int(row["fold"]) - 1
        if np.any(self.fold_of < 0):
            raise AssertionError("Incomplete released fold assignments")
        self._build_superset()
        # The superset stores at most 24 candidates. Counts can saturate at 24,
        # but capped identities/distances and every eligibility test (3 or 5)
        # remain identical to production.
        base_cap, base_dist, base_num = self.geometry(parameters("BASE"))
        if not (np.array_equal(base_cap, nb["cap_idx"][self.bearing]) and
                np.allclose(base_dist, nb["cap_dist"][self.bearing], equal_nan=True) and
                np.array_equal(base_num >= 3, nb["num_around"][self.bearing] >= 3) and
                np.array_equal(base_num >= 5, nb["num_around"][self.bearing] >= 5)):
            cap_diff = int(np.any(base_cap != nb["cap_idx"][self.bearing], axis=1).sum())
            num_diff = int((base_num != nb["num_around"][self.bearing]).sum())
            dist_diff = float(np.nanmax(np.abs(base_dist - nb["cap_dist"][self.bearing])))
            raise AssertionError("BASE geometry mismatch: cap_rows=%d num_rows=%d max_dist=%g" %
                                 (cap_diff, num_diff, dist_diff))
        self.error_cache = {0.0: self.grid_err_base}

    def _build_superset(self):
        blat, blon = self.glat[self.bearing], self.glon[self.bearing]
        self.super_idx = np.full((self.n_eff, 24), -1, dtype=int)
        self.super_dist = np.full((self.n_eff, 24), np.nan)
        for v, gi in enumerate(self.bearing):
            cand = np.flatnonzero((np.abs(blat - self.glat[gi]) <= 1.44) &
                                  (np.abs(blon - self.glon[gi]) <= 1.44))
            if not len(cand):
                continue
            dist = haversine(self.glat[gi], self.glon[gi], blat[cand], blon[cand])
            keep = np.lexsort((cand, dist))[:24]
            self.super_idx[v, :len(keep)] = cand[keep]
            self.super_dist[v, :len(keep)] = dist[keep]

    def geometry(self, p):
        cap = np.full((self.n_eff, p["cap"]), -1, dtype=int)
        dist = np.full((self.n_eff, p["cap"]), np.nan)
        num = np.zeros(self.n_eff, dtype=int)
        blat, blon = self.glat[self.bearing], self.glon[self.bearing]
        for v, gi in enumerate(self.bearing):
            ci = self.super_idx[v]
            ci = ci[ci >= 0]
            cd = self.super_dist[v, :len(ci)]
            inside = (np.abs(blat[ci] - self.glat[gi]) <= p["radius"]) & \
                     (np.abs(blon[ci] - self.glon[gi]) <= p["radius"])
            ci, cd = ci[inside], cd[inside]
            num[v] = len(ci)
            m = min(p["cap"], len(ci))
            cap[v, :m], dist[v, :m] = ci[:m], cd[:m]
        return cap, dist, num

    def grid_error(self, rho):
        rho = float(rho)
        if rho in self.error_cache:
            return self.error_cache[rho]
        cache = self.out / ("grid_err_rho%03d.npy" % round(rho * 100))
        if cache.is_file():
            result = np.load(cache)
            self.error_cache[rho] = result
            return result
        with open_nc(STATION3596_NC) as ds:
            x = np.ma.filled(ds.variables["sm"][:], np.nan).astype(float)
            station_ymd = ds.variables["time_ymd"][:].astype(int)
            nstation = len(ds.dimensions["station"])
        if x.shape[0] != nstation:
            x = x.T
        day_index = {int(day): i for i, day in enumerate(station_ymd)}
        missing_days = [int(day) for day in self.common_ymd if int(day) not in day_index]
        if missing_days:
            raise AssertionError("Station archive lacks analysis dates: %s" % missing_days)
        x = x[:, [day_index[int(day)] for day in self.common_ymd]]
        sigma = np.full(nstation, 0.047)
        tca_path = Path(STATION3596_NC).parent / "station_tca_error_3596.csv"
        with tca_path.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                value = float(row["station_sigma"])
                if np.isfinite(value) and value > 0:
                    sigma[int(row["station_index"])] = value
        result = np.full((self.n_eff, self.n_days), np.nan)
        for v, gi in enumerate(self.bearing):
            mem = np.asarray(self.members[gi], dtype=int)
            valid = np.isfinite(x[mem]) & np.isfinite(sigma[mem, None]) & (sigma[mem, None] > 0)
            n = valid.sum(axis=0)
            sum_sigma = np.where(valid, sigma[mem, None], 0.0).sum(axis=0)
            sum_square = np.where(valid, sigma[mem, None] ** 2, 0.0).sum(axis=0)
            good = n > 0
            variance = ((1 - rho) * sum_square[good] + rho * sum_sigma[good] ** 2) / n[good] ** 2
            result[v, good] = np.sqrt(np.maximum(variance, 0))
        if rho == 0 and not np.allclose(result, self.grid_err_base, equal_nan=True, atol=1e-10, rtol=1e-8):
            raise AssertionError("rho=0 grid error reconstruction differs from production")
        np.save(cache, result)
        self.error_cache[rho] = result
        return result

    def build_cdf(self, p, cap, fold=None):
        raw = self.smap_bg[:, self.bearing]
        cdf = raw.copy()
        for v in range(self.n_eff):
            if np.isfinite(raw[:, v]).sum() <= p["min_hist"]:
                continue
            ci = cap[v]
            ci = ci[ci >= 0]
            if fold is not None:
                ci = ci[self.fold_of[ci] != fold]
            if len(ci) < p["min_nbr"]:
                continue
            cdf[:, v] = cdf_adjust(raw[:, v], self.grid_obs[ci, :].flatten(), p["poly_order"])
        return cdf

    def apply_oi(self, p, cdf, target, cap, dist, grid_err, strict_fold=False):
        background = cdf[:, target].T.copy()
        obs = self.grid_obs[target]
        fusion = background.copy()
        active, neighbour_lists, distance_lists = [], [], []
        for row, v in enumerate(target):
            ci = cap[v]
            ci = ci[ci >= 0]
            cd = dist[v, :len(ci)]
            if strict_fold:
                ci, cd = ci[self.fold_of[ci] != self.fold_of[v]], cd[self.fold_of[ci] != self.fold_of[v]]
            if len(ci) >= p["min_nbr"] and np.isfinite(cdf[:, v]).sum() > p["min_hist"]:
                active.append(row)
                neighbour_lists.append(ci)
                distance_lists.append(cd)
        active = np.asarray(active, dtype=int)
        if len(active):
            bp, bmat, cvec = build_geometry(
                self.glat, self.glon, self.bearing, self.bearing[target[active]],
                neighbour_lists, distance_lists, p["a"])
            safe = np.where(bp >= 0, bp, 0)
            for di in range(self.n_days):
                target_bg = cdf[di, target[active]]
                inc, _, _, _ = dynamic_oi_day(
                    target_bg, self.grid_obs[safe, di], cdf[di, safe],
                    grid_err[safe, di], bp, bmat, cvec, p["sigma_b"])
                fusion[active, di] = target_bg + inc
        return background, fusion, obs

    def run(self, tag):
        p = parameters(tag)
        cap, dist, _ = self.geometry(p)
        grid_err = self.grid_error(p["rho"])
        print("\n=== %s %s ===" % (tag, p), flush=True)
        paths = {kind: self.out / ("%s_%s.npz" % (tag, kind)) for kind in ("ins", "cv")}
        if all(path.is_file() for path in paths.values()):
            print("checkpoint exists; skipping", flush=True)
            return
        # In-sample production-domain diagnostic.
        if not paths["ins"].is_file():
            cdf = self.build_cdf(p, cap)
            target = np.arange(self.n_eff)
            b, f, o = self.apply_oi(p, cdf, target, cap, dist, grid_err)
            valid = np.isfinite(o) & np.isfinite(b) & np.isfinite(f)
            full = metrics(o[valid], b[valid], f[valid])
            np.savez(paths["ins"], B=b.astype(np.float32), F=f.astype(np.float32), valid=valid)
            with (self.out / ("%s_ins_full.json" % tag)).open("w", encoding="utf-8") as fh:
                json.dump(full, fh, indent=2)
            print("in-sample", json.dumps(full), flush=True)
            del cdf, b, f, o, valid
        # Strict station-component CV; every fold has its own CDF and dynamic R.
        if not paths["cv"].is_file():
            ball = np.full((self.n_eff, self.n_days), np.nan)
            fall = np.full_like(ball, np.nan)
            for fold in range(10):
                target = np.flatnonzero(self.fold_of == fold)
                cdf = self.build_cdf(p, cap, fold)
                b, f, o = self.apply_oi(p, cdf, target, cap, dist, grid_err, strict_fold=True)
                ball[target], fall[target] = b, f
                print("%s fold %d targets=%d" % (tag, fold + 1, len(target)), flush=True)
            valid = np.isfinite(self.grid_obs) & np.isfinite(ball) & np.isfinite(fall)
            full = metrics(self.grid_obs[valid], ball[valid], fall[valid])
            np.savez(paths["cv"], B=ball.astype(np.float32), F=fall.astype(np.float32), valid=valid)
            with (self.out / ("%s_cv_full.json" % tag)).open("w", encoding="utf-8") as fh:
                json.dump(full, fh, indent=2)
            print("strict CV", json.dumps(full), flush=True)

    def verify_base(self):
        checks = [
            (self.out / "BASE_cv_full.json", Path(WORK) / "cv_station_grouped_summary.json", "overall"),
            (self.out / "BASE_ins_full.json", Path(WORK) / "direct_gridcell_metrics.csv", None),
        ]
        for generated, reference, key in checks:
            if not generated.is_file() or not reference.is_file():
                continue
            with generated.open(encoding="utf-8") as fh:
                got = json.load(fh)
            if key:
                with reference.open(encoding="utf-8") as fh:
                    want = json.load(fh)[key]
                assert got["N"] == want["fusion"]["N"]
                for system in ("background", "fusion"):
                    for metric in ("MAE", "RMSE", "CORR"):
                        assert abs(got[system][metric] - want[system][metric]) < 1e-7, (system, metric)
            else:
                with reference.open(newline="", encoding="utf-8-sig") as fh:
                    rows = {r["system"]: r for r in csv.DictReader(fh)}
                assert got["N"] == int(rows["Fusion"]["N"])
                for system, name in (("background", "SMAP_CDF"), ("fusion", "Fusion")):
                    for metric in ("MAE", "RMSE", "CORR"):
                        assert abs(got[system][metric] - float(rows[name][metric])) < 1e-7, (system, metric)
        print("BASE matches both released direct and strict-CV references", flush=True)

    def finalize(self, tags):
        common = {}
        for kind in ("ins", "cv"):
            mask = None
            for tag in tags:
                with np.load(self.out / ("%s_%s.npz" % (tag, kind))) as z:
                    valid = z["valid"]
                mask = valid.copy() if mask is None else (mask & valid)
            common[kind] = mask
        rows = []
        for tag in tags:
            p = parameters(tag)
            row = dict(tag=tag, radius=p["radius"], cap=p["cap"], a=p["a"],
                       sigma_b=p["sigma_b"], min_hist=p["min_hist"],
                       min_nbr=p["min_nbr"], poly_order=p["poly_order"], rho=p["rho"])
            for kind in ("ins", "cv"):
                with np.load(self.out / ("%s_%s.npz" % (tag, kind))) as z:
                    b, f = z["B"], z["F"]
                mask = common[kind]
                m = metrics(self.grid_obs[mask], b[mask], f[mask])
                for prefix, source in (("maeB", ("background", "MAE")), ("maeF", ("fusion", "MAE")),
                                       ("rmseB", ("background", "RMSE")), ("rmseF", ("fusion", "RMSE")),
                                       ("rB", ("background", "CORR")), ("rF", ("fusion", "CORR"))):
                    row[kind + "_" + prefix] = m[source[0]][source[1]]
                row[kind + "_N"] = m["N"]
                row[kind + "_improve_mae_pct"] = m["mae_reduction_pct"]
            rows.append(row)
            print("common", tag, row["cv_N"], row["cv_maeB"], row["cv_maeF"], flush=True)
        fieldnames = list(rows[0])
        with (self.out / "sensitivity_common_metrics.csv").open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        summary = {
            "formal_v2_kernel": True, "fold_specific_cdf": True,
            "dynamic_daily_R": True, "day_stride": 1, "n_days": self.n_days,
            "common_mask": "intersection of finite observation/background/analysis masks across all listed variants",
            "common_insample_N": int(common["ins"].sum()), "common_cv_N": int(common["cv"].sum()),
            "variants": tags,
        }
        with (self.out / "sensitivity_summary.json").open("w", encoding="utf-8") as fh:
            json.dump(summary, fh, indent=2)
        print(json.dumps(summary, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path(WORK) / "sensitivity_v6")
    parser.add_argument("--tags", default=",".join(VARIANTS))
    parser.add_argument("--finalize", action="store_true")
    args = parser.parse_args()
    tags = args.tags.split(",")
    unknown = set(tags) - set(VARIANTS)
    if unknown:
        raise ValueError("Unknown tags: %s" % sorted(unknown))
    engine = Sensitivity(args.output_dir)
    if args.finalize:
        engine.finalize(tags)
        return
    for tag in tags:
        engine.run(tag)
    if "BASE" in tags:
        engine.verify_base()


if __name__ == "__main__":
    main()
