function t_aerr
%T_AERR  Reliability of the posterior error field against held-out truth.
%   Uses the fold-specific-CDF CV (cv3_fold.mat): for each validation cell the
%   median posterior Aerr is compared with the mean absolute and the RMS error
%   of the analysis at that cell over all held-out days. Cells are then binned
%   by Aerr decile to test whether the error field discriminates.
P = paths_full(); w0 = warning('off','all');
C = load(fullfile(P.work,'cv3_fold.mat'));
a  = C.cellAerr; e = C.cellAbsErr; rms = C.cellRMS; n = C.cellN; bias = C.cellErr;
q  = isfinite(a) & isfinite(e) & n >= 30;
a = a(q); e = e(q); rms = rms(q); n = n(q); bias = bias(q);
N = numel(a);
fprintf('cells=%d, held-out days per cell: median %.0f\n', N, median(n));

fprintf('\n--- overall ---\n');
fprintf('Aerr median %.4f | |error| median %.4f | RMS error median %.4f\n', ...
    median(a), median(e), median(rms));
fprintf('ratio median(|error|)/Aerr = %.2f\n', median(e./a));
fprintf('Pearson r(Aerr, mean|error|) = %.3f   r(Aerr, RMS error) = %.3f\n', ...
    pcorr(a,e), pcorr(a,rms));
fprintf('Spearman rho(Aerr, mean|error|) = %.3f\n', spear(a,e));

fprintf('\n--- Aerr deciles ---\n');
fprintf(' decile   n    Aerr    mean|err|   RMS err   ratio   BIAS\n');
ed = prctile(a, 0:10:100);
for k = 1:10
    lo = ed(k); hi = ed(k+1);
    if k == 10, m = a >= lo & a <= hi; else, m = a >= lo & a < hi; end
    if sum(m) < 5, continue; end
    fprintf(' %2d/10  %4d   %.4f   %.4f     %.4f    %.2f   %+.4f\n', ...
        k, sum(m), median(a(m)), mean(e(m)), mean(rms(m)), mean(e(m))/median(a(m)), mean(bias(m)));
end

fprintf('\n--- Aerr quintiles of the fraction of pairs within Aerr ---\n');
for k = 1:5
    lo = prctile(a, (k-1)*20); hi = prctile(a, k*20);
    if k == 5, m = a >= lo & a <= hi; else, m = a >= lo & a < hi; end
    if sum(m) < 5, continue; end
    fprintf(' Q%d: Aerr %.4f - %.4f, mean|err| %.4f, mean|err|/Aerr %.2f\n', ...
        k, lo, hi, mean(e(m)), mean(e(m))/mean(a(m)));
end

save(fullfile(P.work,'aerr_reliability.mat'),'a','e','rms','n','bias');
warning(w0);
end

function c = pcorr(x,y)
x = x(:) - mean(x); y = y(:) - mean(y);
c = sum(x.*y)/sqrt(sum(x.^2)*sum(y.^2));
end

function rho = spear(x,y)
rx = tiedrank_(x); ry = tiedrank_(y);
rho = pcorr(rx,ry);
end

function r = tiedrank_(v)
v = v(:); [s,o] = sort(v); r = zeros(size(v)); n = numel(v);
i = 1;
while i <= n
    j = i;
    while j < n && s(j+1) == s(i), j = j + 1; end
    r(o(i:j)) = (i+j)/2;
    i = j + 1;
end
end
