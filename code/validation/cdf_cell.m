function c = cdf_cell(smapRow, staRef, npoly)
% CDF_CDF  Exact reimplementation of lib/cdf_adjusting2.m with a settable
% polynomial order. Column/row broadcasting bug in the earlier draft is gone.
out = nan(size(smapRow));
s = double(smapRow(:));  s(s>0.6 | s<0.04) = nan; s = s(~isnan(s));
t = double(staRef(:));   t(t>0.6 | t<0.04) = nan; t = t(~isnan(t));
if numel(s) < 10 || numel(t) < 10, c = out; return; end
ps = (1:numel(s)).'/(numel(s)+1);
pt = (1:numel(t)).'/(numel(t)+1);
try
    sInt = interp1(ps, sort(s), pt, 'linear', 'extrap');   % numel(t) x 1
catch
    c = out; return;
end
sInt = sInt(:);
DIFF = sort(t) - sInt;                 % same length: both numel(t) x 1
COEFF = polyfit(sInt, DIFF, npoly);
y = polyval(COEFF, double(smapRow)) + double(smapRow);
y(y>1 | y<0) = nan;                    % identical to cdf_adjusting2
c = y;
end
