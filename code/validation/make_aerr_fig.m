function make_aerr_fig()
%MAKE_AERR_FIG  Figure S4: reliability of the posterior error field.
P = paths_full(); w0 = warning('off','all');
C = load(fullfile(P.work,'cv3_fold.mat'));
a = C.cellAerr; e = C.cellAbsErr; rms = C.cellRMS; n = C.cellN;
q = isfinite(a) & isfinite(e) & n >= 30;
a = a(q); e = e(q); rms = rms(q);

edges = prctile(a, 0:10:100);
mA = zeros(10,1); mE = zeros(10,1); mR = zeros(10,1); nB = zeros(10,1);
for k = 1:10
    if k == 10, m = a >= edges(k) & a <= edges(k+1);
    else,        m = a >= edges(k) & a <  edges(k+1); end
    nB(k) = sum(m); mA(k) = median(a(m)); mE(k) = mean(e(m)); mR(k) = mean(rms(m));
end

f = figure('Visible','off','Units','centimeters','Position',[1 1 17 7]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
plot(1:10, mA, '-o','LineWidth',1.4,'Color',[0.00 0.45 0.70],'MarkerFaceColor',[0.00 0.45 0.70],'MarkerSize',5);
plot(1:10, mE, '-s','LineWidth',1.4,'Color',[0.84 0.37 0.00],'MarkerFaceColor',[0.84 0.37 0.00],'MarkerSize',5);
plot(1:10, mR, '-^','LineWidth',1.4,'Color',[0.30 0.30 0.30],'MarkerFaceColor',[0.30 0.30 0.30],'MarkerSize',5);
hold off; box on; grid on;
set(gca,'XTick',1:10,'FontSize',8,'GridAlpha',0.25);
xlabel('Aerr decile','FontSize',8);
ylabel('m^3 m^{-3}','FontSize',8);
legend({'Aerr (median)','mean |error|','RMS error'},'Location','northwest','FontSize',7,'Box','off');
title('(a) Held-out error by Aerr decile','FontSize',8);

nexttile; hold on;
plot(mA, mE, '-o','LineWidth',1.4,'Color',[0.84 0.37 0.00],'MarkerFaceColor',[0.84 0.37 0.00],'MarkerSize',5);
xl = [0.0365 0.0470];
plot(xl, xl, '--','Color',[0.5 0.5 0.5],'LineWidth',1);
plot(xl, 1.6*xl, ':','Color',[0.20 0.20 0.20],'LineWidth',1.2);
hold off; box on; grid on;
set(gca,'FontSize',8,'GridAlpha',0.25,'XLim',xl,'YLim',[0.055 0.09]);
xlabel('median Aerr','FontSize',8);
ylabel('mean |error|','FontSize',8);
legend({'decile means','1:1','1.6:1'},'Location','northwest','FontSize',7,'Box','off');
title('(b) Calibration','FontSize',8);

exportgraphics(f, 'F:\E盘备份\sm_data_products_autoupdate_project\scripts\casmos_chain_v3\figs\si_v2\figS4_aerr_reliability.png', 'Resolution', 400);
close(f);

fid = fopen('F:\E盘备份\sm_data_products_autoupdate_project\scripts\casmos_chain_v3\figs\si_v2\aerr_deciles.csv','w');
fprintf(fid, 'decile,n_cells,median_Aerr,mean_abs_error,rms_error,ratio\n');
for k=1:10, fprintf(fid,'%d,%d,%.5f,%.5f,%.5f,%.2f\n',k,nB(k),mA(k),mE(k),mR(k),mE(k)/mA(k)); end
fclose(fid);
fprintf('figure and table written\n');
warning(w0);
end
