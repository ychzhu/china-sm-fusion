function m9_cv()
%M9_CV  Leak-free 10-fold CV (SMAP error constant 0.046).
%   Faithful port of s9_cv.py. The fold assignment is the fixed, archived
%   permutation in cv_fold_perm.csv (equals numpy RandomState(2026).permutation
%   (1932)), so both implementations use identical folds.

P = paths_full();
w0=warning('off','all');

A=10.0; THRESH=3; NFOLDS=10;

Sg=load(fullfile(P.work,'grid_prep_full.mat'));
bearing=Sg.bearing; glat=Sg.glat; glon=Sg.glon;
gridErr=Sg.gridErr; gridObs=Sg.gridObs; SMAP_ERR=Sg.smapErr(1);
Sn=load(fullfile(P.work,'neighbors_full.mat'));
capIdx=Sn.capIdx; capDist=Sn.capDist;
smapCdf=load(fullfile(P.work,'smap_cdf_full.mat'),'smapCdf'); smapCdf=smapCdf.smapCdf;
nEff=numel(bearing); nDays=size(gridObs,2);
blat=glat(bearing); blon=glon(bearing);
errO=gridErr(bearing);

% fixed permutation (0-based in CSV -> 1-based)
permT=readtable(fullfile(P.here,'cv_fold_perm.csv'));
perm1=permT.perm0+1;
edges=round(linspace(0,nEff,NFOLDS+1));
foldOf=zeros(nEff,1); foldV=cell(NFOLDS,1);
for k=1:NFOLDS
    m=perm1(edges(k)+1:edges(k+1));
    foldV{k}=m; foldOf(m)=k;
end
fprintf('fold sizes: %s\n',mat2str(cellfun(@numel,foldV)));

% per-bearing-cell geometry (full neighbour set, correlation matrix, corIk)
fullNbr=cell(nEff,1); fullB=cell(nEff,1); fullCorIk=cell(nEff,1);
for v=1:nEff
    gv=bearing(v);
    ci=capIdx(gv,:); ci=ci(ci>0);
    cd=capDist(gv,1:numel(ci));
    la=blat(ci); lo=blon(ci); nA=numel(ci);
    dmat=zeros(nA,nA);
    for r=1:nA
        dmat(r,:)=haversine(la(r),lo(r),la,lo);
    end
    fullNbr{v}=ci;
    fullB{v}=1.0./(A.^dmat);
    fullCorIk{v}=1.0./(A.^cd);
end

allB=[]; allF=[]; allO=[];
for k=1:NFOLDS
    V=foldV{k}; nv=numel(V);
    Bv=nan(nDays,nv); Fv=nan(nDays,nv); Ov=nan(nDays,nv);
    for ii=1:nv
        v=V(ii); gv=bearing(v);
        nbr=fullNbr{v};
        trainSel=foldOf(nbr)~=k;
        nbrT=nbr(trainSel);
        Bv(:,ii)=smapCdf(:,gv);
        Ov(:,ii)=gridObs(v,:).';
        if numel(nbrT)<THRESH
            Fv(:,ii)=Bv(:,ii); continue;
        end
        Bsub=fullB{v}(trainSel,trainSel);
        cIk=fullCorIk{v}(trainSel);
        errCon=errO(nbrT)/SMAP_ERR;
        bad=~isfinite(errCon)|errCon<0; errCon(bad)=1.0;
        w=pinv(Bsub+diag(errCon.^2))*cIk(:);
        for di=1:nDays
            resid=gridObs(nbrT,di)-smapCdf(di,bearing(nbrT)).';
            resid(~isfinite(resid))=0.0;
            Fv(di,ii)=smapCdf(di,gv)+sum(resid(:).*w);
        end
    end
    ok=isfinite(Ov)&isfinite(Bv)&isfinite(Fv);
    allB=[allB;Bv(ok)]; allF=[allF;Fv(ok)]; allO=[allO;Ov(ok)]; %#ok<AGROW>
    fprintf('fold %d: pairs=%d MAE bg %.4f -> OI %.4f\n',k,sum(ok(:)), ...
        mean(abs(Bv(ok)-Ov(ok))),mean(abs(Fv(ok)-Ov(ok))));
end

fprintf('===== 10-fold CV overall =====\n');
fprintf('total pairs=%d\n',numel(allO));
maeB=mean(abs(allB-allO)); maeF=mean(abs(allF-allO));
rmseB=sqrt(mean((allB-allO).^2)); rmseF=sqrt(mean((allF-allO).^2));
fprintf('MAE : bg %.4f -> OI %.4f (improve %.1f%%)\n',maeB,maeF,100*(1-maeF/maeB));
fprintf('RMSE: bg %.4f -> OI %.4f (improve %.1f%%)\n',rmseB,rmseF,100*(1-rmseF/rmseB));
fprintf('CORR: bg %.4f -> OI %.4f\n',corr(allB,allO),corr(allF,allO));
fprintf('BIAS: bg %+.4f -> OI %+.4f\n',mean(allB-allO),mean(allF-allO));
fprintf('pairs improved: %.1f%%\n',100*mean(abs(allF-allO)<abs(allB-allO)));
save(fullfile(P.work,'cv_results_full.mat'),'allB','allF','allO','-v7.3');
warning(w0);
end

function d=haversine(lat1,lon1,lat2,lon2)
p1=lat1*pi/180; p2=lat2*pi/180;
dphi=(lat2-lat1)*pi/180; dl=(lon2-lon1)*pi/180;
a=sin(dphi/2).^2+cos(p1).*cos(p2).*sin(dl/2).^2;
d=(2*asin(sqrt(a)))*180/pi;
end
