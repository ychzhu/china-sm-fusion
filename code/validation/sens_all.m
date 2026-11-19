function sens_all()
%SENS_ALL  Full parameter sweep: stage A (CDF) + stage B (OI) for every
%   variant, then stage C (strict 10-fold CV, day stride 3) for every variant.
tags = {'BASE', ...
        'R070','R090','R126','R144', ...
        'K03','K06','K15','K24', ...
        'A03','A05','A20','A30', ...
        'SE040','SE052','SE060', ...
        'H090','H180', ...
        'P1','P3', ...
        'M05'};
t0 = tic;
for k = 1:numel(tags)
    tag = tags{k};
    fprintf('\n===== %s (%d/%d) =====\n', tag, k, numel(tags));
    t = tic; sens_cdf(tag,1.62,48); fprintf('  A %.1f s\n', toc(t));
    t = tic; sens_oi(tag, false);   fprintf('  B %.1f s\n', toc(t));
end
for k = 1:numel(tags)
    tag = tags{k};
    fprintf('\n----- CV %s (%d/%d) -----\n', tag, k, numel(tags));
    t = tic; sens_cv(tag,3); fprintf('  C %.1f s\n', toc(t));
end
fprintf('\nSWEEP DONE in %.1f min\n', toc(t0)/60);
sens_table();
end
