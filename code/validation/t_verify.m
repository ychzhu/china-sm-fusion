function t_verify()
% Compare this harness's BASE fusion against the released fusion_full.mat.
P = paths_full();
So = load(fullfile(P.work,'fusion_full.mat'));      official = So.fusion;      % nDays x nGrid
Sn = load(fullfile(P.work,'sensoi_BASE.mat'));     mine = Sn.fusion;
doOI = Sn.doOI;
fprintf('official size=%s  mine=%s\n', mat2str(size(official)), mat2str(size(mine)));
upd = find(doOI);
D = double(official(:,upd)) - double(mine(:,upd));
ok = isfinite(D);
fprintf('updated cells=%d, compared values=%d\n', numel(upd), sum(ok(:)));
fprintf('max |diff| = %.3e, mean |diff| = %.3e, median |diff| = %.3e\n', ...
    max(abs(D(ok))), mean(abs(D(ok))), median(abs(D(ok))));
fprintf('99th pct |diff| = %.3e\n', prctile(abs(D(ok)),99));
fprintf('share |diff| > 1e-4 = %.4f\n', mean(abs(D(ok))>1e-4));
% non-updated cells should also agree (both = background)
d2 = double(official(:,~doOI)) - double(mine(:,~doOI));
ok2 = isfinite(d2);
fprintf('non-updated cells max |diff| = %.3e\n', max(abs(d2(ok2))));
end
