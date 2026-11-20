function m2_stage2_tca()
%M2_STAGE2_TCA  Triple Collocation error estimation on the 3596 domain.
%   Faithful MATLAB port of s2_stage2_tca.py: closed-form TCA error variances,
%   one global common axis, station_index/grid_index written 0-based.

P = paths_full();
w0=warning('off','all');

sYmd = station_axis_ymd(); nDay=numel(sYmd);

clean=double(load_mat(P,'obj_clean3596.mat','clean'));
lat=load_mat(P,'lat3596.mat','lat');
ids=load_mat(P,'ids3596.mat','allIds');
nSta=numel(ids);

glat=ncread(P.fusion3596Nc,'lat');
glon=ncread(P.fusion3596Nc,'lon');
nGrid=numel(glat);

smapFile=fullfile(P.smap,'smapOI_2day_20150331_20261017.nc');
mSm=read_tgrid(smapFile,'sm'); mYmd=read_ymd(smapFile);
eraFile=fullfile(P.reference,'era5_gcs_swvl1_36km_20150331_20261001.nc');
eSm=read_tgrid(eraFile,'sm'); eYmd=read_ymd(eraFile);

% nearest grid for each station (euclidean on lat/lon, matches cKDTree)
s2g=knnsearch([glat(:),glon(:)],[lat(:),load_mat(P,'lon3596.mat','lon')],'Distance','euclidean');

% global common axis: station intersect smap intersect era, sorted
sMap=containers.Map(sYmd,1:numel(sYmd));
mMap=containers.Map(mYmd,1:numel(mYmd));
eMap=containers.Map(eYmd,1:numel(eYmd));
common=sort(intersect(intersect(sYmd,mYmd),eYmd));
sIdx=cellfun(@(s)sMap(s),common);
mIdx=cellfun(@(s)mMap(s),common);
eIdx=cellfun(@(s)eMap(s),common);
fprintf('common days: %d %s -> %s\n',numel(common),common{1},common{end});

rows=[]; failN=0;
for k=1:nSta
    g=s2g(k);
    x1=clean(k,sIdx).';
    x2=mSm(mIdx,g);
    x3=eSm(eIdx,g);
    ok=isfinite(x1)&isfinite(x2)&isfinite(x3);
    if sum(ok)<100
        failN=failN+1; continue;
    end
    C=cov([x1(ok),x2(ok),x3(ok)]);
    q1=C(1,1)-C(1,2)*C(1,3)/C(2,3);
    q2=C(2,2)-C(1,2)*C(2,3)/C(1,3);
    q3=C(3,3)-C(1,3)*C(2,3)/C(1,2);
    if q1<=0 || q2<=0 || q3<=0
        failN=failN+1; continue;
    end
    rows=[rows; k-1, g-1, sqrt(q1),sqrt(q2),sqrt(q3), sum(ok)]; %#ok<AGROW>
end

t=array2table(rows,'VariableNames',{'station_index','grid_index','station_sigma','smap_sigma','era5_sigma','n_common'});
writetable(t,fullfile(P.work,'tca3596.csv'));

fprintf('stations with TCA: %d (failed/insufficient: %d)\n',height(t),failN);
pr=@(x,p)prctile(x,p);
for j=3:5
    name={'station','smap','era5'}; v=rows(:,j);
    fprintf('  %-7s median=%.4f mean=%.4f p10=%.4f p90=%.4f min=%.4f max=%.4f\n', ...
        name{j-2},median(v),mean(v),pr(v,10),pr(v,90),min(v),max(v));
end
fprintf('common days median=%.0f min=%.0f max=%.0f\n',median(rows(:,6)),min(rows(:,6)),max(rows(:,6)));
fprintf('saved tca3596.csv\n');
warning(w0);
end

function v=load_mat(P,file,var)
v=load(fullfile(P.work,file),var); v=v.(var);
end
