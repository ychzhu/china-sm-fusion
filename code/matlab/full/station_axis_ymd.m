function axisYmd = station_axis_ymd()
%STATION_AXIS_YMD Continuous declared station calendar for V2.1.
%   V2.0 incorrectly inferred this axis from an incomplete legacy NetCDF,
%   silently dropping 2023-01-01 and 2024-01-01 from annual source CSVs.
days = datetime(2015,1,1):caldays(1):datetime(2026,10,1);
axisYmd = cellstr(string(days, 'yyyyMMdd'));
end
