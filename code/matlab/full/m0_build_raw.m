function m0_build_raw()
%M0_BUILD_RAW  Build the 3596-station raw matrix (authoritative 4290-day axis).
%   MATLAB port of s0_build_raw.py. Reads the 12 yearly station CSVs and
%   gathers every station onto the released 3596 day axis.

P = paths_full();

% authoritative day axis from the released 3596 station product
axisYmd = read_ymd(P.station3596Nc);
nDay = numel(axisYmd);
dayPos = containers.Map(axisYmd, 1:nDay);
fprintf('authoritative day axis: %d %s -> %s\n', nDay, axisYmd{1}, axisYmd{end});

% 12 yearly CSVs
files = dir(fullfile(P.station, 'staD10_yearly', 'staD10_2*.csv'));
[~, ord] = sort({files.name});
files = files(ord);
nYear = numel(files);

meta = containers.Map('KeyType','double','ValueType','any');
for yi = 1:nYear
    f = fullfile(files(yi).folder, files(yi).name);
    d = readtable(f, 'TreatAsMissing', {'--',''});
    ids = int64(d.Station_ID_d);
    la = double(d.Lat); lo = double(d.Lon); pr = d.Province;
    for r = 1:numel(ids)
        if ~meta.isKey(double(ids(r)))
            meta(double(ids(r))) = {la(r), lo(r), pr{r}};
        end
    end
    fprintf('read %s rows=%d\n', files(yi).name, size(d,1));
end

allIds = sort(cellfun(@double, meta.keys));
nSta = numel(allIds);
staPos = containers.Map(allIds, 1:nSta);
lat = zeros(nSta,1); lon = zeros(nSta,1); prov = strings(nSta,1);
for k = 1:nSta
    v = meta(allIds(k));
    lat(k)=v{1}; lon(k)=v{2}; prov(k)=string(v{3});
end
fprintf('unique stations: %d\n', nSta);

raw = nan(nSta, nDay, 'single');
presentYear = false(nSta, nYear);
for yi = 1:nYear
    f = fullfile(files(yi).folder, files(yi).name);
    d = readtable(f, 'TreatAsMissing', {'--',''});
    vnames = d.Properties.VariableNames;
    dcols = vnames(startsWith(vnames, 'D'));
    ids = int64(d.Station_ID_d);
    rows = arrayfun(@(s) staPos(double(s)), ids);
    presentYear(rows, yi) = true;
    for c = 1:numel(dcols)
        dn = dcols{c}(2:end);
        if ~dayPos.isKey(dn); continue; end
        cj = dayPos(dn);
        vals = double(d.(dcols{c}));
        fin = isfinite(vals);
        for rr = find(fin(:)).'
            raw(rows(rr), cj) = single(vals(rr));
        end
    end
end

fprintf('raw finite values: %d\n', sum(isfinite(raw(:))));
fprintf('raw range: %.4f .. %.4f\n', min(raw(:)), max(raw(:)));
fprintf('raw values outside [0,0.6]: %d\n', sum(raw(:)<0 | raw(:)>0.6));

save(fullfile(P.work,'raw3596.mat'), 'raw','-v7.3');
save(fullfile(P.work,'ids3596.mat'), 'allIds');
save(fullfile(P.work,'lat3596.mat'), 'lat');
save(fullfile(P.work,'lon3596.mat'), 'lon');
save(fullfile(P.work,'present_year.mat'), 'presentYear');
metaTbl = table(allIds(:), lat, lon, prov, ...
    'VariableNames', {'cma_id','lat','lon','province'});
writetable(metaTbl, fullfile(P.work,'stations3596_meta.csv'));
fprintf('saved work_full raw artifacts\n');
end
