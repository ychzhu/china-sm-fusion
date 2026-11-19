function sens_table()
%SENS_TABLE  Collect the sweep into one CSV and print it.
P = paths_full();
tags = {'BASE','R070','R090','R126','R144','K03','K06','K15','K24', ...
        'A03','A05','A20','A30','SE040','SE052','SE060','H090','H180','P1','P3','M05'};
rows = {};
for k = 1:numel(tags)
    tag = tags{k};
    fo = fullfile(P.work,sprintf('sensoi_%s.mat',tag));
    fc = fullfile(P.work,sprintf('senscv_%s.mat',tag));
    if ~isfile(fo) && ~isfile(fc), continue; end
    r = struct('tag',tag,'radius',NaN,'maxNear',NaN,'a',NaN,'smapErr',NaN, ...
        'minHist',NaN,'polyOrder',NaN,'minNbr',NaN,'oiCells',NaN, ...
        'insitu_mae',NaN,'insitu_r',NaN,'insitu_bias',NaN, ...
        'cv_maeB',NaN,'cv_maeF',NaN,'cv_rmseF',NaN,'cv_rB',NaN,'cv_rF',NaN, ...
        'cv_impMae',NaN,'cv_closer',NaN);
    if isfile(fo)
        S = load(fo); V = S.V;
        r.radius=V.radius; r.maxNear=V.maxNear; r.a=V.a; r.smapErr=V.smapErr;
        r.minHist=V.minHist; r.polyOrder=V.polyOrder; r.minNbr=V.minNbr;
        r.oiCells = S.covStruct.cells;
        r.insitu_mae=S.dFix.mae; r.insitu_r=S.dFix.corr; r.insitu_bias=S.dFix.bias;
    end
    if isfile(fc)
        C = load(fc);
        r.cv_maeB=C.maeB; r.cv_maeF=C.maeF; r.cv_rmseF=C.rmseF;
        r.cv_rB=C.corrB; r.cv_rF=C.corrF; r.cv_impMae=C.impMae; r.cv_closer=100*C.fracCloser;
    end
    rows{end+1} = r; %#ok<AGROW>
end
T = struct2table([rows{:}]);
writetable(T, fullfile(P.work,'sens_summary.csv'));
disp(T(:,{'tag','oiCells','insitu_mae','insitu_r','cv_maeF','cv_rF','cv_impMae','cv_closer'}));
end
