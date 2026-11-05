function m10_compare_refs()
%M10_COMPARE_REFS  Station-level comparison against independent references.
%   Faithful port of s10_compare_refs.py. Five systems (GLDAS, ERA5, SMAP raw,
%   SMAP CDF, Fusion) are sampled at each station's nearest grid and compared
%   with the objective-clean station observations on the common axis.

P = paths_full();
w0=warning('off','all');

Sg=load(fullfile(P.work,'grid_prep_full.mat'));
glat=Sg.glat; glon=Sg.glon; common=Sg.common;
smapBg=Sg.smapBg;
smapCdf=load(fullfile(P.work,'smap_cdf_full.mat'),'smapCdf'); smapCdf=smapCdf.smapCdf;
fusion=load(fullfile(P.work,'fusion_full.mat'),'fusion'); fusion=fusion.fusion;

slat=load(fullfile(P.work,'lat3596.mat'),'lat'); slat=slat.lat;
slon=load(fullfile(P.work,'lon3596.mat'),'lon'); slon=slon.lon;
cleaned=load(fullfile(P.work,'obj_clean3596.mat'),'clean'); cleaned=double(cleaned.clean);
staYmd=read_ymd(P.station3596Nc);
nSta=numel(slat); nCommon=numel(common);

% station observations aligned to common
staMap=containers.Map(staYmd,1:numel(staYmd));
sSel=cellfun(@(s)staMap(s),common);
staObs=cleaned(:,sSel);

% nearest grid for each station
ng=knnsearch([glat(:),glon(:)],[slat(:),slon(:)],'Distance','euclidean');

% independent references aligned to common
era5=align_ref(P,fullfile(P.reference,'era5_gcs_swvl1_36km_20150331_20261001.nc'),common);
gldas=align_ref(P,fullfile(P.reference,'gldas_noah_soilm_36km_20150331_20261001.nc'),common);

systems=struct( ...
    'GLDAS',gldas(:,ng).', 'ERA5',era5(:,ng).', ...
    'SMAP_raw',smapBg(:,ng).', 'SMAP_CDF',smapCdf(:,ng).', ...
    'Fusion',fusion(:,ng).');
names=fieldnames(systems);

fprintf('system vs station observations (pooled station-day pairs)\n');
for q=1:numel(names)
    a=systems.(names{q});
    valid=isfinite(a)&isfinite(staObs);
    x=a(valid); y=staObs(valid); d=x-y;
    fprintf('%-9s N=%9d MAE=%.4f RMSE=%.4f BIAS=%+.4f CORR=%.4f\n', ...
        names{q},numel(x),mean(abs(d)),sqrt(mean(d.^2)),mean(d),corr(x,y));
end
warning(w0);
end

function aligned=align_ref(P,file,common)
g=read_tgrid(file,'sm'); ymd=read_ymd(file);
rmap=containers.Map(ymd,1:numel(ymd));
sel=cellfun(@(s)rmap(s),common);
aligned=g(sel,:);
end
