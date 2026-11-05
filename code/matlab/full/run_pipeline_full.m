function run_pipeline_full()
%RUN_PIPELINE_FULL  Run the full MATLAB 3596 three-stage QC + fusion chain.
%
%   Stages mirror the Python authoritative pipeline:
%     m0  build raw 3596 matrix
%     m1  Stage 1 objective hard-fault screening
%     m2  Stage 2 Triple Collocation error estimation
%     m3  Stage 3 error-driven QC
%     m4  station->grid aggregation and cell error
%     m5  neighbour search
%     m6  CDF bias correction
%     m7  OI fusion
%     m8  direct validation
%     m9  10-fold cross-validation
%     m10 independent reference comparison

P = paths_full();
addpath(P.here);
addpath(fullfile(P.codeDir,'lib'));

stages = {'m0_build_raw','m1_stage1','m2_stage2_tca','m3_stage3', ...
          'm4_grid_agg','m5_neighbors','m6_cdf','m7_fusion', ...
          'm8_eval_direct','m9_cv','m10_compare_refs'};
t0 = tic;
for k = 1:numel(stages)
    fprintf('\n########## running %s ##########\n', stages{k});
    feval(stages{k});
end
fprintf('\nFULL MATLAB PIPELINE DONE in %.1f min\n', toc(t0)/60);
end
