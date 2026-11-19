function t_diag()
% (1) CV sample check, (2) coverage breakdown of the released product.
P = paths_full();
Sg = load(fullfile(P.work,'grid_prep_full.mat'));
gridObs = double(Sg.gridObs); nEff = size(gridObs,1); nDays = size(gridObs,2);
fprintf('bearing cells=%d, days=%d, station-days finite=%d\n', nEff, nDays, sum(isfinite(gridObs(:))));

Si = load(fullfile(P.work,'sensoi_BASE.mat'));
fprintf('CV scored cells=%d\n', Si.nCv);
fprintf('expected pairs at stride 3 on scored cells: %d\n', ...
    sum(sum(isfinite(gridObs(Si.cvBearing, 1:3:nDays)))) );

% coverage of the released product
Sn = load(fullfile(P.work,'smap_cdf_full.mat')); smapCdf = Sn.smapCdf;
Sf = load(fullfile(P.work,'fusion_full.mat')); fusion = Sf.fusion;
doOI = Si.doOI;
nGrid = numel(doOI);
src = zeros(nGrid,1);   % 0 = no value, 1 = raw SMAP (no CDF), 2 = CDF only, 3 = OI-updated
finBg  = sum(isfinite(smapCdf),1);
cdfApplied = abs(smapCdf - Sg.smapBg) > 0;      % sparse test
for gi = 1:nGrid
    if finBg(gi) == 0, src(gi) = 0; continue; end
    if doOI(gi), src(gi) = 3; continue; end
    if cdfApplied(gi), src(gi) = 2; else, src(gi) = 1; end
end
fprintf('\n--- grid-cell provenance (7313 cells) ---\n');
fprintf('OI-updated (station-informed) : %5d (%.1f%%)\n', sum(src==3), 100*mean(src==3));
fprintf('CDF-corrected only           : %5d (%.1f%%)\n', sum(src==2), 100*mean(src==2));
fprintf('raw SMAP, untouched          : %5d (%.1f%%)\n', sum(src==1), 100*mean(src==1));
fprintf('no background value          : %5d (%.1f%%)\n', sum(src==0), 100*mean(src==0));

% grid-day provenance
cells3 = find(src==3); cells2 = find(src==2); cells1 = find(src==1); cells0 = find(src==0);
nDayFin = sum(isfinite(fusion),1);
fprintf('\n--- grid-cell-days ---\n');
tot = numel(fusion);
fprintf('total cells x days = %d\n', tot);
fprintf('OI-updated days    = %d (%.1f%% of all)\n', sum(nDayFin(cells3)), 100*sum(nDayFin(cells3))/tot);
fprintf('CDF-only days      = %d (%.1f%%)\n', sum(nDayFin(cells2)), 100*sum(nDayFin(cells2))/tot);
fprintf('raw SMAP days      = %d (%.1f%%)\n', sum(nDayFin(cells1)), 100*sum(nDayFin(cells1))/tot);
fprintf('missing days       = %d (%.1f%%)\n', tot - sum(isfinite(fusion(:))), 100*(tot-sum(isfinite(fusion(:))))/tot);
fprintf('valid values in released fusion = %d\n', sum(isfinite(fusion(:))));
end
