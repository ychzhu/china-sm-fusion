function make_si_figs(outDir)
%MAKE_SI_FIGS  Figure S1 (coverage) and Figure S3 (within-cell independence).
if nargin < 1 || isempty(outDir)
    outDir = 'C:\Users\zhuyo\AppData\Local\Temp\smsens\out';
end
if ~isfolder(outDir), mkdir(outDir); end
P = paths_full(); w0 = warning('off','all');
Sg = load(fullfile(P.work,'grid_prep_full.mat'));
glat = double(Sg.glat(:)); glon = double(Sg.glon(:));
gridObs = double(Sg.gridObs); nEff = size(gridObs,1); nDays = size(gridObs,2);
bearing = Sg.bearing(:); members = Sg.members;
Si = load(fullfile(P.work,'sensoi_BASE.mat')); doOI = Si.doOI;
smapCdf = double(load(fullfile(P.work,'smap_cdf_full.mat'),'smapCdf').smapCdf);
nGrid = numel(glat); bPosMap = -ones(nGrid,1); bPosMap(bearing) = (1:nEff)';

%% ---------------- Figure S1: coverage ----------------
fracValid = sum(isfinite(smapCdf),1).' / nDays;
fracOI = zeros(nGrid,1);
gOI = find(doOI);
for k = 1:numel(gOI)
    gi = gOI(k);
    fracOI(gi) = sum(isfinite(smapCdf(:,gi)))/nDays;
end
updB = find(doOI(bearing)); nUpd = numel(updB);
cntDay = zeros(nDays,1);
for d = 1:nDays, cntDay(d) = sum(isfinite(gridObs(updB,d))); end
xlims = [72 136]; ylims = [17 55];

f = figure('Visible','off','Units','centimeters','Position',[1 1 17 12]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
scatter(glon, glat, 7, fracValid, 'filled');
cb = colorbar; cb.Label.String = 'fraction of days';
colormap(gca, flipud(gray)); clim([0 1]); xlim(xlims); ylim(ylims);
box on; set(gca,'FontSize',7,'Layer','top');
xlabel('Longitude (deg E)','FontSize',7); ylabel('Latitude (deg N)','FontSize',7);
title('(a) Finite soil moisture','FontSize',8);

nexttile;
scatter(glon, glat, 7, fracOI, 'filled');
cb = colorbar; cb.Label.String = 'fraction of days';
colormap(gca, flipud(gray)); clim([0 1]); xlim(xlims); ylim(ylims);
box on; set(gca,'FontSize',7,'Layer','top');
xlabel('Longitude (deg E)','FontSize',7); ylabel('Latitude (deg N)','FontSize',7);
title('(b) Days with an OI update','FontSize',8);

nexttile;
plot(1:nDays, cntDay, 'Color',[0 0.35 0.6], 'LineWidth',0.8); box on; grid on; set(gca,'GridAlpha',0.2);
set(gca,'FontSize',7);
xlabel('Day index (from 31 March 2015)','FontSize',7);
ylabel('Station cells with a valid update','FontSize',7);
title('(c) Updated station cells per day','FontSize',8);

nexttile;
scatter(glon(doOI), glat(doOI), 7, [0.90 0.55 0.10], 'filled');
hold on;
gB = setdiff(1:nGrid, find(doOI));
scatter(glon(gB), glat(gB), 7, [0.70 0.70 0.70], 'filled');
hold off; box on; xlim(xlims); ylim(ylims); set(gca,'FontSize',7,'Layer','top');
xlabel('Longitude (deg E)','FontSize',7); ylabel('Latitude (deg N)','FontSize',7);
title('(d) Cells with (orange) and without (grey) an OI update','FontSize',7);
exportgraphics(f, fullfile(outDir,'figS1_coverage.png'), 'Resolution', 400);
close(f);

fid = fopen(fullfile(outDir,'coverage_summary.csv'),'w');
fprintf(fid,'quantity,value\n');
fprintf(fid,'cells_total,%d\n', nGrid);
fprintf(fid,'cells_with_oi,%d\n', sum(doOI));
fprintf(fid,'cells_background_only,%d\n', sum(~doOI & fracValid>0));
fprintf(fid,'cells_without_background,%d\n', sum(fracValid==0));
fprintf(fid,'gridcell_days_total,%d\n', nGrid*nDays);
fprintf(fid,'gridcell_days_oi,%d\n', sum(sum(isfinite(smapCdf(:,gOI)))));
fprintf(fid,'gridcell_days_all_finite,%d\n', sum(isfinite(smapCdf(:))));
fclose(fid);

%% ---------------- Figure S3: within-cell independence ----------------
medBy = nan(16,1); q1By = nan(16,1); q3By = nan(16,1); nBy = zeros(16,1);
for k = 2:16
    v = pairdev_k(k, nGrid, members, gridObs, bPosMap);
    nBy(k) = numel(v);
    if numel(v) >= 20
        medBy(k) = median(v); q1By(k) = prctile(v,25); q3By(k) = prctile(v,75);
    end
end
nCell = cellfun(@numel, members);

f = figure('Visible','off','Units','centimeters','Position',[1 1 17 6.5]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
for k = 2:16
    if isnan(medBy(k)), continue; end
    plot([k k],[q1By(k) q3By(k)],'-','Color',[0.78 0.78 0.78],'LineWidth',4);
end
plot(2:16, medBy(2:16), 'o-','Color',[0.84 0.37 0.00], ...
    'MarkerFaceColor',[0.84 0.37 0.00],'MarkerSize',5,'LineWidth',1.3);
yline(0.046,'--','single-station TCA error 0.046','Color',[0.25 0.25 0.25],'FontSize',7);
hold off; box on; grid on; set(gca,'FontSize',8,'GridAlpha',0.2);
xlabel('Stations in the cell','FontSize',8);
ylabel('RMS pair difference / sqrt(2)','FontSize',8);
title('(a) Within-cell pair-difference error','FontSize',8);
nexttile;
histogram(nCell, 0:1:16, 'FaceColor',[0.20 0.45 0.65], 'EdgeColor','none');
box on; grid on; set(gca,'FontSize',8,'GridAlpha',0.2);
xlabel('Stations per cell','FontSize',8); ylabel('Cells','FontSize',8);
title('(b) Stations per station-bearing cell','FontSize',8);
exportgraphics(f, fullfile(outDir,'figS3_within_cell.png'), 'Resolution', 400);
close(f);

fid = fopen(fullfile(outDir,'within_cell.csv'),'w');
fprintf(fid,'stations_in_cell,n_pairs,median_implied_error,q25,q75\n');
for k = 1:16
    fprintf(fid,'%d,%d,%.4f,%.4f,%.4f\n',k,nBy(k),medBy(k),q1By(k),q3By(k));
end
fclose(fid);
fprintf('figures and tables written to %s\n', outDir);
warning(w0);
end

function v = pairdev_k(kwant, nGrid, members, gridObs, bPosMap)
v = [];
for gi = 1:nGrid
    m = members{gi};
    if numel(m) ~= kwant, continue; end
    pp = bPosMap(m); pp = pp(pp > 0);
    if numel(pp) ~= numel(m), continue; end
    X = gridObs(pp, :);
    ok = isfinite(X);
    for a = 1:numel(m)-1
        for b = a+1:numel(m)
            both = ok(a,:) & ok(b,:);
            if sum(both) < 30, continue; end
            dv = X(a,both) - X(b,both);
            v(end+1,1) = sqrt(mean(dv.^2)/2);
        end
    end
end
end
