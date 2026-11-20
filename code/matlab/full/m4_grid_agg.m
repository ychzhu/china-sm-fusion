function m4_grid_agg()
%M4_GRID_AGG  Station -> 36-km grid aggregation + cell error (rho=0).
%   Faithful port of s4_grid_agg.py. Bearing = every grid with >=1 station;
%   grid_obs is indexed by bearing position on the common (station intersect
%   SMAP) axis; smap_bg (full grid) is stored for the downstream stages.

P = paths_full();
w0=warning('off','all');

GRID_RES=P.gridRes; FALLBACK=P.staFallbackErr; SMAP_ERR=P.smapTcaErr;

% grid geometry (3596 fusion product)
glat=ncread(P.fusion3596Nc,'lat'); glon=ncread(P.fusion3596Nc,'lon');
nGrid=numel(glat);

% objective-cleaned station matrix + geometry
cleaned=load(fullfile(P.work,'obj_clean3596.mat'),'clean'); cleaned=double(cleaned.clean);
slat=load(fullfile(P.work,'lat3596.mat'),'lat'); slat=slat.lat;
slon=load(fullfile(P.work,'lon3596.mat'),'lon'); slon=slon.lon;
nSta=numel(slat);
staYmd=station_axis_ymd();

% SMAP 2-day background
smapFile=fullfile(P.smap,'smapOI_2day_20150331_20261017.nc');
smapBg=read_tgrid(smapFile,'sm'); smapYmd=read_ymd(smapFile);

% box membership: for each grid, stations within +/-0.09 in both lat and lon
members=cell(nGrid,1);
for gi=1:nGrid
    members{gi}=find(abs(slat-glat(gi))<=GRID_RES & abs(slon-glon(gi))<=GRID_RES);
end
nWith=cellfun(@numel,members);
bearing=find(nWith>0);
fprintf('bearing grids=%d (old 1714), assigned stations=%d\n',numel(bearing),sum(nWith>0));

% new TCA per-station sigma, keyed by station_index (0-based CSV)
tca=readtable(fullfile(P.work,'tca3596.csv'));
staSigma=nan(nSta,1);
staSigma(tca.station_index+1)=tca.station_sigma;

% static cell error rho=0
gridErr=nan(nGrid,1);
for b=1:numel(bearing)
    gi=bearing(b);
    s=staSigma(members{gi});
    s(~isfinite(s) | s<=0)=FALLBACK;
    gridErr(gi)=sqrt(sum(s.^2))/numel(members{gi});
end
fprintf('cell error median=%.4f mean=%.4f (old median 0.0386)\n',median(gridErr(bearing)),mean(gridErr(bearing)));

% common axis: station intersect smap, sorted
common=sort(intersect(staYmd,smapYmd));
sp=containers.Map(staYmd,1:numel(staYmd));
mp=containers.Map(smapYmd,1:numel(smapYmd));
selS=cellfun(@(s)sp(s),common);
selM=cellfun(@(s)mp(s),common);
staAligned=cleaned(:,selS);
smapBg=smapBg(selM,:);

nBearing=numel(bearing); nCommon=numel(common);
gridObs=nan(nBearing,nCommon);
for b=1:nBearing
    m=members{bearing(b)};
    gridObs(b,:)=mean(staAligned(m,:),1,'omitnan');
end
fprintf('common analysis days=%d grid_obs finite=%d\n',nCommon,sum(isfinite(gridObs(:))));
fprintf('stations per bearing cell: median=%.0f mean=%.1f max=%d\n', ...
    median(nWith(bearing)),mean(nWith(bearing)),max(nWith(bearing)));

smapErr=repmat(SMAP_ERR,nCommon,1);
save(fullfile(P.work,'grid_prep_full.mat'),'members','bearing','nWith','gridErr', ...
    'gridObs','glat','glon','common','smapBg','smapErr','-v7.3');
fprintf('saved grid_prep_full.mat\n');
warning(w0);
end
