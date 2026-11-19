function make_si_fig3(outDir)
%MAKE_SI_FIG3  Figure S3: within-cell station error independence.
%   The pair-difference statistic is the RMS of co-located station differences.
%   Under the model xi = theta + ei with Var(ei) = sigma^2 and inter-station
%   error correlation rho, the expected value is sigma*sqrt(2(1-rho)).  Comparing
%   the observed value with the TCA single-station sigma therefore tests rho,
%   and the expected values for rho = 0 and rho = 1 are drawn for reference.
if nargin < 1 || isempty(outDir), outDir = 'C:\Users\zhuyo\AppData\Local\Temp\smsens\out'; end
P = paths_full(); w0 = warning('off','all');
Sg = load(fullfile(P.work,'grid_prep_full.mat'));
glat = double(Sg.glat(:)); glon = double(Sg.glon(:));

sta  = double(ncread(P.station3596Nc,'sm'));       % station x day
slat = ncread(P.station3596Nc,'lat');
slon = ncread(P.station3596Nc,'lon');
nSta = size(sta,1);
SIG  = 0.046;                                      % TCA single-station error

cellOf = zeros(nSta,1);
for s = 1:nSta
    d = abs(glat - slat(s)) <= 0.18 & abs(glon - slon(s)) <= 0.18;
    k = find(d,1);
    if ~isempty(k), cellOf(s) = k; end
end
uc = unique(cellOf(cellOf>0));
nPerCell = histcounts(cellOf(cellOf>0), [uc; max(uc)+1]);

pairRMS = []; pairN = []; cellPairRMS = []; cellPairN = [];
for k = 1:numel(uc)
    g = uc(k);
    s = find(cellOf == g);
    if numel(s) < 2, continue; end
    X = sta(s,:);
    v = []; vn = [];
    for a = 1:numel(s)-1
        for b = a+1:numel(s)
            both = isfinite(X(a,:)) & isfinite(X(b,:));
            if sum(both) < 30, continue; end
            dv = X(a,both) - X(b,both);
            v(end+1,1) = sqrt(mean(dv.^2));      % RMS pair difference
            vn(end+1,1) = sum(both);
        end
    end
    if isempty(v), continue; end
    pairRMS = [pairRMS; v]; pairN = [pairN; vn];
    cellPairRMS(end+1,1) = median(v);
    cellPairN(end+1,1) = numel(s);
end

rho_hat = 1 - (median(pairRMS)/SIG)^2 / 2;
fprintf('stations assigned %d, cells %d (1:%d 2:%d >=3:%d, max %d)\n', ...
    sum(cellOf>0), numel(uc), sum(nPerCell==1), sum(nPerCell==2), sum(nPerCell>=3), max(nPerCell));
fprintf('pairs %d, median shared days %.0f\n', numel(pairRMS), median(pairN));
fprintf('RMS pair difference: median %.4f, p10 %.4f, p90 %.4f\n', ...
    median(pairRMS), prctile(pairRMS,10), prctile(pairRMS,90));
fprintf('expected if rho=0: %.4f ; if rho=1: %.4f (sigma=%.3f)\n', ...
    SIG*sqrt(2), SIG, SIG);
fprintf('implied rho = %.3f\n', rho_hat);

f = figure('Visible','off','Units','centimeters','Position',[1 1 17 6.5]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
histogram(pairRMS, 0:0.01:0.30, 'FaceColor',[0.20 0.45 0.65], 'EdgeColor','none');
hold on;
xl = [0 0.30];
plot([SIG*sqrt(2) SIG*sqrt(2)], [0 1], '-', 'Color',[0.84 0.37 0.00], 'LineWidth',1.4);
plot([SIG SIG], [0 1], '-', 'Color',[0.45 0.45 0.45], 'LineWidth',1.4);
plot([median(pairRMS) median(pairRMS)], [0 1], '--', 'Color',[0.00 0.00 0.00], 'LineWidth',1.4);
hold off; box on; grid on; set(gca,'FontSize',8,'GridAlpha',0.2);
xlim(xl); ylim([0 inf]);
legend({'observed pairs','\rho = 0','\rho = 1','observed median'}, ...
    'FontSize',7,'Box','off','Location','northeast');
xlabel('RMS difference of co-located stations (m^3 m^{-3})','FontSize',8);
ylabel('Pairs','FontSize',8);
title('(a) Co-located pair differences','FontSize',8);

nexttile; hold on;
plot([SIG sqrt(2)*SIG], [0 1], '--', 'Color',[0.55 0.55 0.55],'LineWidth',1.1);
scatter(cellPairN, cellPairRMS, 26, [0.84 0.37 0.00], 'filled', 'MarkerFaceAlpha',0.6);
yline(median(pairRMS), ':', 'Color',[0.20 0.20 0.20],'LineWidth',1.2);
hold off; box on; grid on; set(gca,'FontSize',8,'GridAlpha',0.2);
set(gca,'XScale','log','XLim',[1.8 16],'YLim',[0 0.20]);
xlabel('Stations in the cell','FontSize',8);
ylabel('Median RMS pair difference (m^3 m^{-3})','FontSize',8);
title('(b) Estimate versus stations per cell','FontSize',8);
exportgraphics(f, fullfile(outDir,'figS3_within_cell.png'),'Resolution',400);
close(f);

fid = fopen(fullfile(outDir,'within_cell.csv'),'w');
fprintf(fid,'quantity,value\n');
fprintf(fid,'stations_assigned,%d\n', sum(cellOf>0));
fprintf(fid,'cells_used,%d\n', numel(uc));
fprintf(fid,'cells_1_station,%d\n', sum(nPerCell==1));
fprintf(fid,'cells_2_stations,%d\n', sum(nPerCell==2));
fprintf(fid,'cells_3_or_more,%d\n', sum(nPerCell>=3));
fprintf(fid,'max_stations_in_cell,%d\n', max(nPerCell));
fprintf(fid,'pairs_evaluated,%d\n', numel(pairRMS));
fprintf(fid,'median_shared_days,%.0f\n', median(pairN));
fprintf(fid,'pair_rms_median,%.5f\n', median(pairRMS));
fprintf(fid,'pair_rms_p10,%.5f\n', prctile(pairRMS,10));
fprintf(fid,'pair_rms_p90,%.5f\n', prctile(pairRMS,90));
fprintf(fid,'expected_rho0,%.5f\n', SIG*sqrt(2));
fprintf(fid,'expected_rho1,%.5f\n', SIG);
fprintf(fid,'implied_rho,%.3f\n', rho_hat);
fclose(fid);
warning(w0);
end
