function t_within_cell_v2(outDir)
%T_WITHIN_CELL_V2  Within-cell station error correlation, corrected statistics.
%
%   Model: xi = theta + ei, Var(ei)=sigma^2, corr(ei,ej)=rho.
%   Then RMS(xi-xj) = sigma*sqrt(2*(1-rho)), i.e.  rho = 1 - (RMS/sigma)^2 / 2.
%
%   Three things are checked that the earlier version got wrong:
%     1. rho is inferred with the correct factor of 2;
%     2. pairs are stratified by their actual great-circle separation, because
%        the 0.18 deg box admits separations of up to ~28 km and a pair that is
%        not co-located also contains a sub-grid representativeness difference;
%     3. the implied sigma is compared with the TCA single-station sigma, so the
%        two hypotheses (rho>0 versus larger sigma) can be separated.
if nargin < 1 || isempty(outDir), outDir = 'C:\Users\zhuyo\AppData\Local\Temp\smsens\out'; end
if ~isfolder(outDir), mkdir(outDir); end
P = paths_full(); w0 = warning('off','all');

Sg = load(fullfile(P.work,'grid_prep_full.mat'));
glat = double(Sg.glat(:)); glon = double(Sg.glon(:));
nFull = numel(glat);

sta  = double(ncread(P.station3596Nc,'sm'));        % station x day
slat = ncread(P.station3596Nc,'lat');
slon = ncread(P.station3596Nc,'lon');
nSta = numel(slat);
SIG_TCA = 0.046;

% --- assign each station to its containing cell (production rule: 0.18 deg box) ---
cellOf = zeros(nSta,1);
for s = 1:nSta
    d = abs(glat - slat(s)) <= 0.18 & abs(glon - slon(s)) <= 0.18;
    k = find(d,1);
    if ~isempty(k), cellOf(s) = k; end
end
uc = unique(cellOf(cellOf>0));

% --- collect every co-located pair with its true separation ---
pRMS = []; pSep = []; pN = []; pCell = [];
for k = 1:numel(uc)
    g = uc(k);
    s = find(cellOf == g);
    if numel(s) < 2, continue; end
    X = sta(s,:);
    for a = 1:numel(s)-1
        for b = a+1:numel(s)
            both = isfinite(X(a,:)) & isfinite(X(b,:));
            if sum(both) < 30, continue; end
            dv = X(a,both) - X(b,both);
            pRMS(end+1,1) = sqrt(mean(dv.^2));
            pSep(end+1,1) = hav(slat(s(a)),slon(s(a)),slat(s(b)),slon(s(b)));
            pN(end+1,1)   = sum(both);
            pCell(end+1,1) = numel(s);
        end
    end
end
fprintf('co-located pairs: %d, median shared days %d\n', numel(pRMS), round(median(pN)));
fprintf('pair separation: median %.1f km, p90 %.1f km, max %.1f km\n', ...
    median(pSep)*111.2, prctile(pSep,90)*111.2, max(pSep)*111.2);

% --- rho from the correct formula, using the TCA sigma ---
rhoAll = 1 - (pRMS/SIG_TCA).^2 / 2;
fprintf('\n--- pooled over all pairs (TCA sigma = %.3f) ---\n', SIG_TCA);
fprintf('  RMS pair difference median %.4f -> rho = %.3f\n', median(pRMS), median(rhoAll));
fprintf('  implied sigma if rho = 0: %.4f  (vs TCA %.3f)\n', median(pRMS)/sqrt(2), SIG_TCA);

% --- stratify by separation: close pairs isolate instrument error ---
edges = [0 5 10 20 30 50 80 200] / 111.2;   % convert km to degrees
fprintf('\n--- stratify by separation ---\n');
fprintf('  distance       n     RMS median   implied sigma (rho=0)   rho (TCA sigma)\n');
strat = [];
for e = 1:numel(edges)-1
    m = pSep >= edges(e) & pSep < edges(e+1);
    if sum(m) < 20, continue; end
    r = median(pRMS(m)); sg = r/sqrt(2); rh = 1 - (r/SIG_TCA)^2/2;
    fprintf('  %5.0f-%-5.0f km %5d    %.4f        %.4f                %+.3f\n', ...
        edges(e)*111.2, edges(e+1)*111.2, sum(m), r, sg, rh);
    strat(end+1,:) = [edges(e)*111.2, edges(e+1)*111.2, sum(m), r, sg, rh];
end

% --- nearest-neighbour pairs only: the best available "co-located" test ---
[~, ordIdx] = sort(pSep);
bestN = min(500, numel(ordIdx));
sel = ordIdx(1:bestN);
fprintf('\n--- the %d closest pairs ---\n', bestN);
fprintf('  median separation %.1f km, RMS pair difference %.4f\n', ...
    median(pSep(sel))*111.2, median(pRMS(sel)));
fprintf('  implied sigma if rho = 0 : %.4f  (TCA %.3f)\n', median(pRMS(sel))/sqrt(2), SIG_TCA);
fprintf('  implied rho if sigma = TCA: %+.3f\n', 1 - (median(pRMS(sel))/SIG_TCA)^2/2);

% --- bootstrap CI on rho ---
rng(7); B = 2000; bs = zeros(B,1);
n = numel(pRMS);
for b = 1:B
    ii = randi(n, n, 1);
    bs(b) = 1 - (median(pRMS(ii))/SIG_TCA)^2/2;
end
fprintf('\n  rho = %+.3f  (bootstrap 5-95%%: %+.3f to %+.3f)\n', ...
    median(rhoAll), prctile(bs,5), prctile(bs,95));

fid = fopen(fullfile(outDir,'within_cell_v2.csv'),'w');
fprintf(fid,'distance_lo_km,distance_hi_km,n_pairs,rms_median,implied_sigma_rho0,rho_tca_sigma\n');
for k = 1:size(strat,1), fprintf(fid,'%.0f,%.0f,%d,%.5f,%.5f,%+.3f\n',strat(k,:)); end
fclose(fid);

f = figure('Visible','off','Units','centimeters','Position',[1 1 17 6.5]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
plot([SIG_TCA*1 SIG_TCA*1],[0 1],'-','Color',[0.30 0.30 0.30],'LineWidth',1.4);
plot([SIG_TCA*sqrt(2) SIG_TCA*sqrt(2)],[0 1],'-','Color',[0.84 0.37 0.00],'LineWidth',1.4);
histogram(pRMS, 0:0.01:0.30, 'FaceColor',[0.20 0.45 0.65],'EdgeColor','none');
hold off; box on; grid on; set(gca,'FontSize',8,'GridAlpha',0.2);
xlim([0 0.30]);
legend({'rho = 1 (sigma)','rho = 0 (sigma sqrt 2)','observed pairs'},'FontSize',7,'Box','off','Location','northeast');
xlabel('RMS difference of co-located stations','FontSize',8); ylabel('Pairs','FontSize',8);
title('(a) All pairs','FontSize',8);

nexttile; hold on;
if ~isempty(strat)
    mid = mean(strat(:,[1 2]),2);
    semilogy(mid, strat(:,4), 'o-','Color',[0.20 0.45 0.65],'MarkerFaceColor',[0.20 0.45 0.65],'MarkerSize',5);
end
yline(SIG_TCA*sqrt(2),'--','rho = 0 expectation','Color',[0.84 0.37 0.00],'FontSize',7);
hold off; box on; grid on; set(gca,'FontSize',8,'GridAlpha',0.2,'XScale','log');
xlabel('Station separation (km)','FontSize',8);
ylabel('RMS pair difference','FontSize',8);
title('(b) Against separation','FontSize',8);
exportgraphics(f, fullfile(outDir,'figS2_within_cell_v2.png'),'Resolution',400);
close(f);
fprintf('\nfigure written\n');
warning(w0);
end

function d = hav(lat1,lon1,lat2,lon2)
p1 = lat1*pi/180; p2 = lat2*pi/180;
dphi = (lat2-lat1)*pi/180; dl = (lon2-lon1)*pi/180;
aa = sin(dphi/2).^2 + cos(p1).*cos(p2).*sin(dl/2).^2;
d = (2*asin(sqrt(aa)))*180/pi;
end
