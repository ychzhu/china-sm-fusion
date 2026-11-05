function m3_stage3()
%M3_STAGE3  Robust sigma threshold + error-driven QC (no circularity).
%   Faithful port of s3_stage3.py. Extreme sigma alone never deletes: it hard-
%   drops only with an objective Stage-1 fault; otherwise kept and down-weighted
%   in OI. The final flag equals the objective flag unchanged.

P = paths_full();
w0=warning('off','all');

tca=readtable(fullfile(P.work,'tca3596.csv'));
obj=readtable(fullfile(P.work,'obj_station_year.csv'),'TextType','string');
objFlag=load(fullfile(P.work,'obj_flag3596.mat'),'objFlag'); objFlag=objFlag.objFlag;
ids=load(fullfile(P.work,'ids3596.mat'),'allIds'); ids=ids.allIds;
nSta=numel(ids);

% robust sigma statistics from the station TCA distribution
sig=tca.station_sigma;
med=median(sig);
madv=median(abs(sig-med));
robustSD=1.4826*madv;
sigmaExtreme=med+3*robustSD;
fprintf('station sigma robust stats: n=%d median=%.4f MAD=%.4f robustSD=%.4f\n',numel(sig),med,madv,robustSD);
fprintf('p90=%.4f p95=%.4f p99=%.4f max=%.4f\n',prctile(sig,90),prctile(sig,95),prctile(sig,99),max(sig));
fprintf('NEW SIGMA_EXTREME=%.4f (old 0.1129)\n',sigmaExtreme);

% per-station sigma keyed by station_index (0-based in CSV)
sigByStation=nan(nSta,1);
sigByStation(tca.station_index+1)=tca.station_sigma;
idPos=containers.Map(double(ids),1:nSta);

nRow=height(obj);
fCma=zeros(nRow,1); fYear=zeros(nRow,1); fRep=zeros(nRow,1); fExp=zeros(nRow,1);
fComp=zeros(nRow,1); fMask=zeros(nRow,1); fSig=nan(nRow,1);
fDec=strings(nRow,1); fReas=strings(nRow,1);

for r=1:nRow
    k=idPos(double(obj.cma_id(r)));
    s=sigByStation(k);
    objectiveDrop=obj.decision(r)=="drop";
    extreme=isfinite(s) && s>sigmaExtreme;
    l3="";
    if extreme && objectiveDrop
        l3="L3_extreme_sigma";
    elseif extreme
        l3="HIGH_SIGMA_kept";
    end
    reasons=string(obj.reasons(r));
    if ismissing(reasons); reasons=""; end
    if strlength(l3)>0 && strlength(reasons)>0 && reasons~="nan"
        reasons=reasons+";"+l3;
    elseif strlength(l3)>0
        reasons=l3;
    end
    fCma(r)=double(obj.cma_id(r)); fYear(r)=double(obj.year(r));
    fRep(r)=double(obj.n_reported(r)); fExp(r)=double(obj.expected_days(r));
    fComp(r)=double(obj.active_completeness(r)); fMask(r)=double(obj.masked_days(r));
    if isfinite(s); fSig(r)=round(s,4); end
    fDec(r)=obj.decision(r); fReas(r)=reasons;
end

final=table(fCma,fYear,fRep,fExp,fComp,fMask,fSig,fDec,fReas, ...
    'VariableNames',{'cma_id','year','n_reported','expected_days','active_completeness', ...
    'masked_days','sigma','decision','reasons'});
writetable(final,fullfile(P.work,'final_station_year.csv'));

% final flag = objective flag (drop/mask unchanged), authoritative
save(fullfile(P.work,'final_flag3596.mat'),'objFlag');

nExtreme=sum(final.sigma>sigmaExtreme);
nKept=sum(contains(final.reasons,"HIGH_SIGMA_kept"));
nDropped=sum(contains(final.reasons,"L3_extreme_sigma"));
fprintf('\nstation-years with extreme sigma: %d (kept/down-weighted=%d, dropped w/ corroboration=%d)\n',nExtreme,nKept,nDropped);
fprintf('FINAL station-years: keep=%d drop=%d (%.1f%%)\n', ...
    sum(final.decision=="keep"),sum(final.decision=="drop"),100*mean(final.decision=="drop"));
fprintf('unique stations with >=1 keep year: %d\n',numel(unique(final.cma_id(final.decision=="keep"))));
warning(w0);
end
