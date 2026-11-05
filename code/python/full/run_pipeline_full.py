# -*- coding: utf-8 -*-
"""Run the full-domain (3596) three-stage QC + fusion pipeline in order.

s0 builds the raw station matrix from the yearly CSVs; s1-s3 are the three
quality-control stages (objective hard-fault screening, TCA, sigma-driven
decisions); s4-s7 rebuild the grid aggregation, neighbours, CDF and OI fusion;
s8-s10 run the direct, cross-validation and reference comparisons.
"""
import os
import runpy

HERE = os.path.dirname(os.path.abspath(__file__))
steps = [
    "s0_build_raw.py", "s1_stage1.py", "s2_stage2_tca.py", "s3_stage3.py",
    "s4_grid_agg.py", "s5_neighbors.py", "s6_cdf.py", "s7_fusion.py",
    "s8_eval_direct.py", "s9_cv.py", "s10_compare_refs.py",
]
for step in steps:
    print("\n########## running %s ##########" % step)
    runpy.run_path(os.path.join(HERE, step), run_name="__main__")
print("\nFULL PIPELINE DONE")
