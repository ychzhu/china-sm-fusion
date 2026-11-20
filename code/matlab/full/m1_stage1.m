function m1_stage1()
%M1_STAGE1  Objective hard-fault screening on the 3596 domain (TCA-independent).
%   MATLAB port of s1_stage1.py. Per-year neighbour filtering uses only the
%   stations present in that year.

P = paths_full();
w0 = warning('off','all');  % silence omitnan-on-empty warnings (Python also warns)

NEIGH_RADIUS=0.75; MIN_NEIGH=3; SEASON_MIN_DAYS=10;
RANGE_LO=0.0; RANGE_HI=0.6; RANGE_BAD_FRAC=0.05;
FLAT_STD=0.005; NEIGH_ACTIVE_STD=0.01;
JUMP_ABS=0.15; JUMP_NEIGH=0.05; JUMP_COUNT=3;
COMP_FLOOR=0.50; L2_MEAN_LO=-0.10; L2_MEAN_HI=0.15; R_THRESH=0.5;
KEEP=0; MASK=1; DROP=2;
NON_STATION_IDS=[1,2,3,11111,22222,33333,55555,66666,99999];

seasonNames={'DJF','MAM','JJA','SON'};
seasonMonths={[12,1,2],[3,4,5],[6,7,8],[9,10,11]};
growing=[false,true,true,true];

axisYmd = station_axis_ymd();
nDay = numel(axisYmd);
yrAll = cellfun(@(s) str2double(s(1:4)), axisYmd);
moAll = cellfun(@(s) str2double(s(5:6)), axisYmd);
years = (min(yrAll):max(yrAll)).';

Sraw=load(fullfile(P.work,'raw3596.mat')); raw=double(Sraw.raw);
Sla=load(fullfile(P.work,'lat3596.mat')); lat=Sla.lat;
Slo=load(fullfile(P.work,'lon3596.mat')); lon=Slo.lon;
Sid=load(fullfile(P.work,'ids3596.mat')); ids=Sid.allIds;
Spr=load(fullfile(P.work,'present_year.mat')); present=Spr.presentYear;
nSta=numel(ids);

globalNeigh = rangesearch([lat,lon],[lat,lon],NEIGH_RADIUS,'Distance','euclidean');

objFlag=zeros(nSta,nDay,'int8');
R=struct('cma_id',{},'year',{},'n_reported',{},'expected_days',{}, ...
    'individual_outage',{},'active_completeness',{},'masked_days',{}, ...
    'decision',{},'reasons',{});
ri=0;
for yi=1:numel(years)
    ymask = yrAll==years(yi);
    yidx = find(ymask);
    for i=1:nSta
        if ~present(i,yi); continue; end
        nonStation = ismember(double(ids(i)), NON_STATION_IDS);
        gn = globalNeigh{i}(:);
        nbr = gn(gn~=i & present(gn,yi));
        x = raw(i,ymask);
        reported = isfinite(x);
        nReported = sum(reported);
        hasNbr = numel(nbr) > MIN_NEIGH;
        dayFlag = zeros(1,numel(yidx),'int8');
        reasons = {};

        if hasNbr
            nbrRep = isfinite(raw(nbr,ymask));
            expected = reported | (mean(nbrRep,1) >= 0.5);
            nExpected = sum(expected);
            outage = sum(expected & ~reported);
            comp = nReported/nExpected;
            nmean = mean(raw(nbr,ymask),1,'omitnan');
        else
            nExpected = numel(yidx);
            outage = sum(~reported);
            comp = nReported/nExpected;
            nmean = nan(1,numel(yidx));
        end
        l0Fail = nExpected>=10 && comp < COMP_FLOOR;

        l1Fail=false; growingStuck=false;
        dx=[NaN,diff(x)]; dn=[NaN,diff(nmean)];
        for sn=1:4
            msel = ismember(moAll(yidx), seasonMonths{sn});
            if sum(msel)==0; continue; end
            xs=x(msel);
            valid=isfinite(xs);
            if sum(valid)<SEASON_MIN_DAYS; continue; end
            xv=xs(valid);
            sidxLocal=find(msel); sidxLocal=sidxLocal(valid);

            bad = xv<RANGE_LO | xv>RANGE_HI;
            if any(bad)
                dayFlag(sidxLocal(bad))=MASK;
                if sum(bad)/sum(valid) > RANGE_BAD_FRAC
                    l1Fail=true; reasons{end+1}='range_bad';
                end
            end

            if std(xv,1) < FLAT_STD
                if sn==1
                    dayFlag(sidxLocal)=MASK;
                    reasons{end+1}='winter_freeze';
                elseif growing(sn)
                    if hasNbr
                        nsub = nmean(msel); nsub = nsub(valid);
                        if std(nsub,1,'omitnan') > NEIGH_ACTIVE_STD
                            growingStuck=true;
                            reasons{end+1}=['growing_stuck_' seasonNames{sn}];
                        end
                    else
                        growingStuck=true;
                        reasons{end+1}=['growing_stuck_' seasonNames{sn}];
                    end
                end
            end

            if hasNbr
                jump = abs(dx)>JUMP_ABS & abs(dn)<JUMP_NEIGH;
                if sum(jump)>=JUMP_COUNT
                    l1Fail=true; reasons{end+1}='jump_bad';
                else
                    dayFlag(find(jump))=MASK;
                end
            end
        end

        l2Fail=false;
        if hasNbr
            mean3 = mean(x,'omitnan') - mean(nmean,'omitnan');
            if mean3<L2_MEAN_LO || mean3>L2_MEAN_HI
                bothm = reported & isfinite(nmean);
                if sum(bothm)>=2 && std(x(bothm),1)>0 && std(nmean(bothm),1)>0
                    if corr(x(bothm).',nmean(bothm).') < R_THRESH
                        l2Fail=true; reasons{end+1}='spatial_bad';
                    end
                elseif abs(mean3)>0.30
                    l2Fail=true; reasons{end+1}='spatial_bad';
                end
            end
        end

        if nonStation; reasons{end+1}='non_station_id'; end
        drop = nonStation || l0Fail || l1Fail || l2Fail || growingStuck;
        if l0Fail; reasons{end+1}='completeness_bad'; end
        if drop; dayFlag(:)=DROP; end
        objFlag(i,ymask)=dayFlag;

        reasons=uniqueStr(reasons);
        ri=ri+1;
        R(ri).cma_id=double(ids(i)); R(ri).year=years(yi);
        R(ri).n_reported=nReported; R(ri).expected_days=nExpected;
        R(ri).individual_outage=outage; R(ri).active_completeness=round(comp,4);
        R(ri).masked_days=sum(dayFlag==MASK);
        R(ri).decision=string(ternary(drop,'drop','keep'));
        R(ri).reasons=strjoin(reasons,';');
    end
end

res=struct2table(R);
writetable(res,fullfile(P.work,'obj_station_year.csv'));
save(fullfile(P.work,'obj_flag3596.mat'),'objFlag');
clean=single(nan(size(raw)));
keepMask = objFlag==KEEP;
clean(keepMask)=single(raw(keepMask));
save(fullfile(P.work,'obj_clean3596.mat'),'clean','-v7.3');

summ=[];
for yi=1:numel(years)
    sub=res(res.year==years(yi),:);
    rd=sum(sub.decision=="drop");
    summ=[summ; years(yi), size(sub,1), sum(sub.decision=="keep"), rd, ...
        round(100*rd/size(sub,1),1), sum(sub.masked_days)]; %#ok<AGROW>
end
summTbl=array2table(summ,'VariableNames',{'year','n_present','n_keep','n_drop','drop_pct','masked_days'});
writetable(summTbl,fullfile(P.work,'obj_yearly_summary.csv'));
disp(summTbl);
fprintf('TOTAL station-years=%d keep=%d drop=%d (%.1f%%)\n', ...
    height(res), sum(res.decision=="keep"), sum(res.decision=="drop"), ...
    100*mean(res.decision=="drop"));
fprintf('masked station-days=%d\n', sum(res.masked_days));
fprintf('clean-for-TCA finite=%d\n', sum(isfinite(clean(:))));
warning(w0);
end

function out=uniqueStr(c)
[~,idx]=unique(c,'stable');
out=c(sort(idx));
end

function v=ternary(cond,a,b)
if cond; v=a; else; v=b; end
end
