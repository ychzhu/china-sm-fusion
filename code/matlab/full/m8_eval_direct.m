function m8_eval_direct()
%M8_EVAL_DIRECT  Full-domain direct (in-sample) validation.
%   Faithful port of s8_eval_direct.py. Compares the background and fusion on
%   the bearing cells with pooled and per-day metrics against grid_obs.

P = paths_full();
w0=warning('off','all');

Sg=load(fullfile(P.work,'grid_prep_full.mat'));
bearing=Sg.bearing; gridObs=Sg.gridObs;
smapCdf=load(fullfile(P.work,'smap_cdf_full.mat'),'smapCdf'); smapCdf=smapCdf.smapCdf;
fusion=load(fullfile(P.work,'fusion_full.mat'),'fusion'); fusion=fusion.fusion;
nDays=size(gridObs,2);

bgB=smapCdf(:,bearing).';   % nBearing x nDays
fuB=fusion(:,bearing).';
obs=gridObs;

names={'background','fusion'}; preds={bgB,fuB};
pooled=cell(1,2); daily=cell(1,2); counts=zeros(1,2);
for q=1:2
    pred=preds{q};
    valid=isfinite(obs)&isfinite(pred);
    counts(q)=sum(valid(:));
    pooled{q}=pooled_metrics(pred(valid),obs(valid));
    dmae=[]; dcorr=[];
    for di=1:nDays
        v=isfinite(obs(:,di))&isfinite(pred(:,di));
        if sum(v)<3; continue; end
        dmae(end+1)=mean(abs(pred(v,di)-obs(v,di))); %#ok<AGROW>
        if std(pred(v,di),1)>0 && std(obs(v,di),1)>0
            dcorr(end+1)=corr(pred(v,di),obs(v,di)); %#ok<AGROW>
        end
    end
    daily{q}=[mean(dmae),mean(dcorr)];
end

for q=1:2
    m=pooled{q};
    fprintf('=== %s === pairs=%d\n',names{q},counts(q));
    fprintf('  pooled MAE=%.4f RMSE=%.4f BIAS=%+.4f CORR=%.4f\n',m(1),m(2),m(3),m(4));
    fprintf('  daily  MAE=%.4f CORR=%.4f\n',daily{q}(1),daily{q}(2));
end
b=pooled{1}; f=pooled{2};
fprintf('in-sample improvement: MAE %.1f%%  RMSE %.1f%%\n',100*(1-f(1)/b(1)),100*(1-f(2)/b(2)));
warning(w0);
end

function m=pooled_metrics(pred,obs)
diff=pred-obs;
pv=pred-mean(pred); ov=obs-mean(obs);
corr=sum(pv.*ov)/sqrt(sum(pv.^2)*sum(ov.^2));
m=[mean(abs(diff)),sqrt(mean(diff.^2)),mean(diff),corr];
end
