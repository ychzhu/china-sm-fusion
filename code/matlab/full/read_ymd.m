function ymd = read_ymd(path)
%READ_YMD  Read the time_ymd variable as a cell array of date strings.
v = ncread(path, 'time_ymd');
ymd = strtrim(cellstr(num2str(v)));
end
