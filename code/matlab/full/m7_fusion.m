function m7_fusion()
%M7_FUSION  Full-domain OI fusion (SMAP background error constant 0.046).
%   Faithful port of s7_fusion.py. Static OI weights use the correlation model
%   B=1/A^d with error conditioning; daily fusion adds the weighted residuals.

P = paths_full();
w0=warning('off','all');

A=10.0; THRESH_NUM=3; MIN_HIST=30;

Sg=load(fullfile(P.work,'grid_prep_full.mat'));
bearing=Sg.bearing; glat=Sg.glat; glon=Sg.glon;
gridErr=Sg.gridErr; gridObs=Sg.gridObs; SMAP_ERR=Sg.smapErr(1);
Sn=load(fullfile(P.work,'neighbors_full.mat'));
numAround=Sn.numAround; capIdx=Sn.capIdx; capDist=Sn.capDist;
smapCdf=load(fullfile(P.work,'smap_cdf_full.mat'),'smapCdf'); smapCdf=smapCdf.smapCdf;
nDays=size(smapCdf,1); nGrid=size(smapCdf,2);

errOBearing=gridErr(bearing);
blat=glat(bearing); blon=glon(bearing);

doOI=false(nGrid,1);
geoBp=cell(nGrid,1); geoW=cell(nGrid,1);
for gi=1:nGrid
    daysValid=sum(isfinite(smapCdf(:,gi)));
    if numAround(gi)<THRESH_NUM || daysValid<=MIN_HIST
        continue;
    end
    bp=capIdx(gi,:); bp=bp(bp>0);
    nA=numel(bp);
    la=blat(bp); lo=blon(bp);
    dmat=zeros(nA,nA);
    for r=1:nA
        dmat(r,:)=haversine(la(r),lo(r),la,lo);
    end
    B=1.0./(A.^dmat);
    corIk=1.0./(A.^capDist(gi,1:nA));
    errCon=errOBearing(bp)/SMAP_ERR;
    bad=~isfinite(errCon)|errCon<0;
    errCon(bad)=1.0;
    staticW=pinv(B+diag(errCon.^2))*corIk(:);
    doOI(gi)=true; geoBp{gi}=bp; geoW{gi}=staticW;
end
fprintf('doOI grids=%d\n',sum(doOI));

fusion=smapCdf;
fusionErr=repmat(SMAP_ERR,nDays,nGrid);
for di=1:nDays
    for gi=1:nGrid
        if ~doOI(gi); continue; end
        bp=geoBp{gi}; w=geoW{gi};
        resid=gridObs(bp,di)-smapCdf(di,bearing(bp)).';
        resid(~isfinite(resid))=0.0;
        fusion(di,gi)=smapCdf(di,gi)+sum(resid(:).*w);
        eDiff=errOBearing(bp)-SMAP_ERR;
        fusionErr(di,gi)=SMAP_ERR+sum(eDiff(:).*w);
    end
end

fprintf('fusion finite=%d (smap_cdf finite=%d)\n',sum(isfinite(fusion(:))),sum(isfinite(smapCdf(:))));
save(fullfile(P.work,'fusion_full.mat'),'fusion','-v7.3');
save(fullfile(P.work,'fusion_full_err.mat'),'fusionErr','-v7.3');
fprintf('saved fusion_full.mat and fusion_full_err.mat\n');
warning(w0);
end

function d=haversine(lat1,lon1,lat2,lon2)
p1=lat1*pi/180; p2=lat2*pi/180;
dphi=(lat2-lat1)*pi/180; dl=(lon2-lon1)*pi/180;
a=sin(dphi/2).^2+cos(p1).*cos(p2).*sin(dl/2).^2;
d=(2*asin(sqrt(a)))*180/pi;
end
