function m6_cdf()
%M6_CDF  Apply per-grid CDF bias correction on the full domain.
%   Faithful port of s6_cdf.py. A grid is corrected when it has >=3 capped
%   neighbours and >30 valid background values; the station reference is the
%   pooled grid_obs of its cap_idx neighbours. Reuses lib/cdf_adjusting2.

P = paths_full();
addpath(fullfile(P.codeDir,'lib'));
w0=warning('off','all');

THRESH_NUM=3; MIN_HIST=30;

Sg=load(fullfile(P.work,'grid_prep_full.mat'));
smapBg=Sg.smapBg; gridObs=Sg.gridObs;
nDays=size(smapBg,1); nGrid=size(smapBg,2);
Sn=load(fullfile(P.work,'neighbors_full.mat'));
numAround=Sn.numAround; capIdx=Sn.capIdx;

smapCdf=smapBg;
applied=0;
for gi=1:nGrid
    valid=isfinite(smapBg(:,gi));
    if numAround(gi)>=THRESH_NUM && sum(valid)>MIN_HIST
        ci=capIdx(gi,:); ci=ci(ci>0);
        staRef=reshape(gridObs(ci,:).',1,[]);
        corrected=cdf_adjusting2(smapBg(:,gi).',staRef);
        smapCdf(:,gi)=corrected(:);
        applied=applied+1;
    end
end

fprintf('CDF applied to %d / %d grids\n',applied,nGrid);
fprintf('smap background finite=%d CDF finite=%d\n',sum(isfinite(smapBg(:))),sum(isfinite(smapCdf(:))));
save(fullfile(P.work,'smap_cdf_full.mat'),'smapCdf','-v7.3');
fprintf('saved smap_cdf_full.mat\n');
warning(w0);
end
