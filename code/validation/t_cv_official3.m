function t_cv_official3(opts)
%T_CV_OFFICIAL3  Strict station-grouped CV, full-fold CDF reconstruction.
%
%   Per fold:
%     1. rebuild the CDF background field on EVERY grid cell from training-fold
%        neighbours only (production rule: +/-1.08 box, 9 nearest, >=3 nbrs,
%        >30 valid days), so background and analysis use one coherent field;
%     2. refit OI weights for every validation cell with same-fold neighbours
%        removed;
%     3. score against the held-out cell aggregate; retain posterior Aerr.
%
%   opts.brMode : 'grid'  (analysis background = production CDF field)
%                 'fold'  (analysis background = fold-specific CDF field)

if nargin < 1 || isempty(opts), opts = struct(); end
if ~isfield(opts,'brMode'), opts.brMode = 'fold'; end
if ~isfield(opts,'onlyFold'), opts.onlyFold = 0; end   % 0 = all folds
if ~isfield(opts,'verbose'), opts.verbose = 1; end

w0 = warning('off','all');
P  = paths_full();
V  = struct('radius',1.08,'near',9,'adecay',10,'smapErr',0.046, ...
            'minHist',30,'minNbr',3,'polyOrder',2);

Sg = load(fullfile(P.work,'grid_prep_full.mat'));
bearing = Sg.bearing(:); glat = double(Sg.glat(:)); glon = double(Sg.glon(:));
gridObs = double(Sg.gridObs); gridErr = double(Sg.gridErr(:));
smapBg  = double(Sg.smapBg);
smapCdfProd = double(load(fullfile(P.work,'smap_cdf_full.mat'),'smapCdf').smapCdf);
if strcmpi(opts.brMode,'grid'), bgField0 = smapCdfProd; else, bgField0 = smapBg; end

nGrid = numel(glat); nEff = size(gridObs,1); nDays = size(gridObs,2);
blat = glat(bearing); blon = glon(bearing); errO = gridErr;
bPosMap = -ones(nGrid,1); bPosMap(bearing) = (1:nEff)';
daysValid = sum(isfinite(smapBg),1);

va = readtable(fullfile(P.root,'datarelease-scienceDB','validation','cv_fold_assignments.csv'));
foldOfBearing = zeros(nEff,1);
for r = 1:height(va)
    foldOfBearing(va.bearing_position(r) + 1) = va.fold(r);
end
nCv = height(va);

S = load(fullfile(P.work,'sens_geo_R162K48.mat'));
capIdx = zeros(nGrid,V.near); capDist = nan(nGrid,V.near);
for gi = 1:nGrid
    ci = S.gIdx(gi,:); ci = ci(ci>0); cd = S.gDist(gi,1:numel(ci));
    keep = abs(blat(ci)-glat(gi)) <= V.radius & abs(blon(ci)-glon(gi)) <= V.radius;
    ci = ci(keep); cd = cd(keep);
    m = min(V.near,numel(ci));
    capIdx(gi,1:m) = ci(1:m); capDist(gi,1:m) = cd(1:m);
end

bgAll = nan(nCv,nDays); fuAll = nan(nCv,nDays); obAll = nan(nCv,nDays); aeAll = nan(nCv,nDays);
folds = unique(va.fold);
if opts.onlyFold > 0, folds = opts.onlyFold; end

for f = folds(:).'
    idx = find(va.fold == f);
    t0 = tic;
    % ---- fold-specific CDF field on all cells ----
    bgField = bgField0;
    nApplied = 0;
    for gi = 1:nGrid
        ci = capIdx(gi,:); ci = ci(ci>0);
        ci = ci(foldOfBearing(ci) ~= f);           % training neighbours only
        if numel(ci) < V.minNbr || daysValid(gi) <= V.minHist, continue; end
        ref = reshape(gridObs(ci,:).',1,[]);
        bgField(:,gi) = cdf_cell(smapBg(:,gi), ref, V.polyOrder);
        nApplied = nApplied + 1;
    end
    if opts.verbose, fprintf('fold %d: %d cells, CDF applied %d ... ', f, numel(idx), nApplied); end

    for ii = 1:numel(idx)
        r  = idx(ii); gv = va.grid_index(r) + 1; v = bPosMap(gv);
        if v <= 0, continue; end
        ci = capIdx(gv,:); ci = ci(ci>0); cd = capDist(gv,1:numel(ci));
        tr = foldOfBearing(ci) ~= f;
        nb = ci(tr); ndk = cd(tr); ndk = ndk(:);

        bgc = bgField(:,gv);
        ob  = reshape(gridObs(v,:), nDays, 1);
        if numel(nb) < V.minNbr
            ana = bgc; errOut = repmat(V.smapErr, nDays, 1);
        else
            nA = numel(nb); la = blat(nb); lo = blon(nb); dmat = zeros(nA,nA);
            for q = 1:nA, dmat(q,:) = hav(la(q),lo(q),la,lo); end
            Cmat = 1.0./(V.adecay.^dmat);
            cvec = 1.0./(V.adecay.^ndk);
            sN = errO(nb); bad = ~isfinite(sN) | sN <= 0; sN(bad) = 0.046;
            errCon = sN/V.smapErr;
            wv = pinv(Cmat + diag(errCon.^2)) * cvec(:);
            resid = gridObs(nb,:).' - bgField(:, bearing(nb));
            resid(~isfinite(resid)) = 0;
            ana = bgc + resid * wv;
            errOut = repmat(V.smapErr*sqrt(max(0, 1 - cvec(:).'*wv)), nDays, 1);
        end
        bgAll(r,:) = reshape(bgc,1,nDays);
        fuAll(r,:) = reshape(ana,1,nDays);
        obAll(r,:) = reshape(ob,1,nDays);
        aeAll(r,:) = reshape(errOut,1,nDays);
    end
    if opts.verbose, fprintf('%.1f s\n', toc(t0)); end
end

% ---- score only the folds that were run ----
runMask = ismember(va.fold, folds);
ok = isfinite(bgAll) & isfinite(fuAll) & isfinite(obAll) & repmat(runMask,1,nDays);
bb = bgAll(ok); ff = fuAll(ok); oo = obAll(ok);

CV = struct();
CV.brMode = opts.brMode; CV.folds = folds(:).';
CV.nPairs = numel(oo);
CV.maeB  = mean(abs(bb-oo));       CV.maeF  = mean(abs(ff-oo));
CV.rmseB = sqrt(mean((bb-oo).^2)); CV.rmseF = sqrt(mean((ff-oo).^2));
CV.corrB = pcorr(bb,oo);           CV.corrF = pcorr(ff,oo);
CV.biasB = mean(bb-oo);            CV.biasF = mean(ff-oo);
CV.fracCloser = mean(abs(ff-oo) < abs(bb-oo));
CV.impMae  = 100*(1-CV.maeF/CV.maeB);
CV.impRmse = 100*(1-CV.rmseF/CV.rmseB);

cellAerr = nan(nCv,1); cellAbsErr = nan(nCv,1); cellRMS = nan(nCv,1); cellErr = nan(nCv,1); cellN = zeros(nCv,1);
cellTrainNbr = zeros(nCv,1);
for r = 1:nCv
    if ~runMask(r), continue; end
    q = isfinite(obAll(r,:)) & isfinite(fuAll(r,:));
    if sum(q) < 30, continue; end
    cellN(r) = sum(q);
    cellAerr(r) = median(aeAll(r,q));
    cellAbsErr(r) = mean(abs(fuAll(r,q)-obAll(r,q)));
    cellRMS(r) = sqrt(mean((fuAll(r,q)-obAll(r,q)).^2));
    cellErr(r) = mean(fuAll(r,q)-obAll(r,q));
    gv = va.grid_index(r)+1; ci = capIdx(gv,:); ci = ci(ci>0);
    cellTrainNbr(r) = sum(foldOfBearing(ci) ~= va.fold(r));
end
CV.cellAerr = cellAerr; CV.cellAbsErr = cellAbsErr; CV.cellRMS = cellRMS;
CV.cellErr = cellErr; CV.cellN = cellN; CV.cellTrainNbr = cellTrainNbr;
qq = isfinite(cellAerr) & isfinite(cellAbsErr);
CV.rAerr = pcorr(cellAerr(qq), cellAbsErr(qq));
CV.rAerrRMS = pcorr(cellAerr(qq), cellRMS(qq));

tag = sprintf('cv3_%s', opts.brMode);
if opts.onlyFold > 0, tag = sprintf('%s_fold%d', tag, opts.onlyFold); end
save(fullfile(P.work,[tag '.mat']),'-struct','CV','-v7.3');

fprintf(['\n[%s] n=%d | MAE %.5f -> %.5f (%.2f%%) | RMSE %.5f -> %.5f (%.2f%%) | ', ...
         'r %.4f -> %.4f | BIAS %+.5f -> %+.5f | closer %.1f%%\n'], ...
    tag, numel(oo), CV.maeB, CV.maeF, CV.impMae, CV.rmseB, CV.rmseF, CV.impRmse, ...
    CV.corrB, CV.corrF, CV.biasB, CV.biasF, 100*CV.fracCloser);
fprintf('Aerr vs mean|error| r = %.3f ; Aerr vs RMS error r = %.3f (n=%d cells)\n', ...
    CV.rAerr, CV.rAerrRMS, sum(qq));
warning(w0);
end

function c = pcorr(x,y)
x = x(:) - mean(x); y = y(:) - mean(y);
c = sum(x.*y)/sqrt(sum(x.^2)*sum(y.^2));
end

function d = hav(lat1,lon1,lat2,lon2)
p1 = lat1*pi/180; p2 = lat2*pi/180;
dphi = (lat2-lat1)*pi/180; dl = (lon2-lon1)*pi/180;
aa = sin(dphi/2).^2 + cos(p1).*cos(p2).*sin(dl/2).^2;
d = (2*asin(sqrt(aa)))*180/pi;
end
