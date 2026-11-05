# -*- coding: utf-8 -*-
"""Stage 3: robust sigma threshold + error-driven QC (no circularity).

High sigma alone never deletes: extreme sigma hard-drops only with an
independent Stage-1 objective fault; otherwise the station is kept and the OI
step down-weights it (1/sigma^2).
"""
from paths_full import *

import os
import numpy as np
import pandas as pd


tca = pd.read_csv(os.path.join(WORK, "tca3596.csv"))
obj = pd.read_csv(os.path.join(WORK, "obj_station_year.csv"))
obj_flag = np.load(os.path.join(WORK, "obj_flag3596.npy"))
sids = np.load(os.path.join(WORK, "ids3596.npy"))

# robust sigma statistics from the new station TCA distribution
sig = tca.station_sigma.to_numpy()
med = np.median(sig)
mad = np.median(np.abs(sig - med))
robust_sd = 1.4826 * mad
sigma_extreme = med + 3 * robust_sd
print("station sigma robust stats: n=%d median=%.4f MAD=%.4f robustSD=%.4f" %
      (len(sig), med, mad, robust_sd))
print("p90=%.4f p95=%.4f p99=%.4f max=%.4f" %
      (np.percentile(sig, 90), np.percentile(sig, 95),
       np.percentile(sig, 99), sig.max()))
print("NEW SIGMA_EXTREME=%.4f (old 0.1129)" % sigma_extreme)

# map per-station sigma
sig_by_station = dict(zip(tca.station_index, tca.station_sigma))
pos = {int(s): k for k, s in enumerate(sids)}

final_rows = []
for _, r in obj.iterrows():
    k = pos[int(r.cma_id)]
    s = sig_by_station.get(k, np.nan)
    objective_drop = r.decision == "drop"
    extreme = np.isfinite(s) and s > sigma_extreme
    l3 = []
    if extreme and objective_drop:
        l3.append("L3_extreme_sigma")
    elif extreme:
        l3.append("HIGH_SIGMA_kept")
    reasons = str(r.reasons)
    if l3 and reasons and reasons != "nan":
        reasons = reasons + ";" + l3[0]
    elif l3:
        reasons = l3[0]
    final_rows.append({
        "cma_id": int(r.cma_id), "year": int(r.year),
        "n_reported": int(r.n_reported), "expected_days": int(r.expected_days),
        "active_completeness": r.active_completeness,
        "masked_days": int(r.masked_days),
        "sigma": round(s, 4) if np.isfinite(s) else np.nan,
        "decision": r.decision, "reasons": reasons})

final = pd.DataFrame(final_rows)
final.to_csv(os.path.join(WORK, "final_station_year.csv"), index=False)

# final flag = objective flag (drop/mask unchanged); keep as authoritative
np.save(os.path.join(WORK, "final_flag3596.npy"), obj_flag)

n_extreme = int(((final.sigma > sigma_extreme)).sum())
n_extreme_kept = int(final.reasons.str.contains("HIGH_SIGMA_kept").sum())
n_extreme_dropped = int(final.reasons.str.contains("L3_extreme_sigma").sum())
print("\nstation-years with extreme sigma: %d (kept/down-weighted=%d, dropped w/ corroboration=%d)" %
      (n_extreme, n_extreme_kept, n_extreme_dropped))
print("FINAL station-years: keep=%d drop=%d (%.1f%%)" %
      ((final.decision == 'keep').sum(), (final.decision == 'drop').sum(),
       100 * (final.decision == 'drop').mean()))
print("unique stations with >=1 keep year:", final[final.decision=='keep'].cma_id.nunique())
