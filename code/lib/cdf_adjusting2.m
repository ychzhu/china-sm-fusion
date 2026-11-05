function smapCDF2 = cdf_adjusting2(smap_thisday, sta_thisday)
%CDF_ADJUSTING2 CDF-based bias correction of a satellite time series
%   towards the station climatology (per-grid, per-day application).
%
%   Method (quantile-mapping variant):
%     1. restrict both series to the physically plausible range
%        [0.04, 0.6] m3 m-3 and keep valid samples only;
%     2. build empirical CDF ranks for both samples (Weibull plotting
%        position);
%     3. interpolate the SMAP sample onto the station CDF ranks;
%     4. fit a 2nd-order polynomial to the rank-mapped difference
%        DIFF = sorted(sta) - interp(smap), and evaluate it on the full
%        SMAP series: corrected = SMAP + polyval(DIFF-fit).
%     5. clip the corrected values to [0, 1].
%
%   INPUT
%     smap_thisday vector of SMAP background values (a grid time series)
%     sta_thisday  vector of station-grid mean observations (same period)
%   OUTPUT
%     smapCDF2     column vector of corrected SMAP (NaN where failed)

% force column vectors so all internal algebra stays orientation-safe
smap_thisday = smap_thisday(:);
sta_thisday  = sta_thisday(:);

% plausible range for volumetric soil moisture
smap_thisday(smap_thisday > 0.6 | smap_thisday < 0.04) = nan;
smap_delNan = smap_thisday(~isnan(smap_thisday));

sta_thisday(sta_thisday > 0.6 | sta_thisday < 0.04) = nan;
sta_delNan = sta_thisday(~isnan(sta_thisday));

% empirical CDF ranks (Weibull plotting position)
pSmap = (1:numel(smap_delNan)).' / (numel(smap_delNan) + 1);
pSta  = (1:numel(sta_delNan)).'  / (numel(sta_delNan) + 1);

% map SMAP samples onto the station CDF ranks
try
    smapInt2sta = interp1(pSmap, sort(smap_delNan), pSta, 'linear', 'extrap');
catch
    smapCDF2 = nan(size(smap_thisday));
    return;
end

% difference between the two CDFs at the station ranks
DIFF = sort(sta_delNan) - smapInt2sta;

% 2nd-order polynomial fit of the difference as a function of the SMAP value
COEFF = polyfit(smapInt2sta, DIFF, 2);

% apply the correction to the full SMAP series and clip to [0,1]
smapCDF2 = polyval(COEFF, smap_thisday) + smap_thisday;
smapCDF2(smapCDF2 > 1 | smapCDF2 < 0) = nan;

end
