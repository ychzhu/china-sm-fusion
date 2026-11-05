function M = read_tgrid(path, varname)
%READ_TGRID  Read a (time,grid) or (time,station) NetCDF variable.
%   ncread returns the transposed (column-major) layout; this returns the
%   documented (nTime x nSite) double matrix with NaN for masked fill values.
M = double(ncread(path, varname));   % (nSite x nTime)
M = M.';                             % (nTime x nSite)
M(isnan(M)) = NaN;                   % normalise any fill
end
