function sens_oi(tag, keepFusion)
%SENS_OI  Stage B: static OI weights, daily fusion, in-sample diagnostics.
%   Arrays:  smapCdf, fusion  (nDays x nGrid)
%            gridObs, smapCdfB  (nBearing x nDays)
%   The fusion is accumulated per neighbour slot, so every operand keeps an
%   explicit (nDays x nCell) shape and no implicit expansion is relied on.
if nargin < 2, keepFusion = false; end
w0 = warning('off','all');
P = paths_full();
Sc = load(fullfile(P.work,sprintf('senscdf_%s.mat',tag)));
V = Sc.V; smapCdf = Sc.smapCdf; capIdx = Sc.capIdx; capDist = Sc.capDist;
numAround = Sc.numAround; applied = Sc.applied;

Sg = load(fullfile(P.work,'grid_prep_full.mat'));
bearing = Sg.bearing; glat = double(Sg.glat); glon = double(Sg.glon);
gridObs = double(Sg.gridObs); gridErr = double(Sg.gridErr);
nGrid = numel(glat); nEff = numel(bearing);
nDays = size(smapCdf,1);
blat = glat(bearing); blon = glon(bearing); errO = gridErr(bearing);
bPosMap = -ones(nGrid,1); bPosMap(bearing) = 1:nEff;
daysValid = sum(isfinite(smapCdf),1);

doOI = false(nGrid,1);
wMat = zeros(nGrid,V.maxNear); posMat = zeros(nGrid,V.maxNear);
for gi = 1:nGrid
    if numAround(gi) < V.minNbr || daysValid(gi) <= V.minHist, continue; end
    bp = capIdx(gi,:); bp = bp(bp>0);
    nA = numel(bp); la = blat(bp); lo = blon(bp); dmat = zeros(nA,nA);
    for r = 1:nA, dmat(r,:) = hav(la(r),lo(r),la,lo); end
    B = 1.0./(V.a.^dmat);
    corIk = 1.0./(V.a.^capDist(gi,1:nA));
    errCon = errO(bp)/V.smapErr;
    bad = ~isfinite(errCon) | errCon < 0; errCon(bad) = 1.0;
    w = pinv(B + diag(errCon.^2)) * corIk(:);
    doOI(gi) = true;
    wMat(gi,1:nA) = w(:); posMat(gi,1:nA) = bp;   % bearing positions
end

baseFile = fullfile(P.work,'sens_base_domain.mat');
if ~isfile(baseFile), error('sens_base_domain.mat missing - run sens_prep first'); end
Bd = load(baseFile); baseMask = Bd.baseMask;

% ---- daily fusion, accumulated slot by slot ----
upd = find(doOI); nUpd = numel(upd);
inc = zeros(nDays, nUpd);                       % daily increment at updated cells
for j = 1:V.maxNear
    sel = find(posMat(upd,j) > 0);
    if isempty(sel), continue; end
    p = posMat(upd(sel), j);                    % bearing positions of the neighbours
    rg = gridObs(p, :).' - smapCdf(:, bearing(p));   % nDays x nSel
    rg(~isfinite(rg)) = 0;
    wj = reshape(wMat(upd(sel), j), 1, []);     % 1 x nSel
    inc(:, sel) = inc(:, sel) + rg .* wj;
end
fusion = smapCdf;
fusion(:, upd) = smapCdf(:, upd) + inc;
clear inc rg;

% ---- in-sample diagnostics ----
updB = doOI(bearing);
baseB = baseMask(bearing);              % restrict the fixed domain to bearing cells
bgB = smapCdf(:,bearing).';
fuB = fusion(:,bearing).';
dIn  = evalmetrics(bgB(updB,:),fuB(updB,:),gridObs(updB,:));
dFix = evalmetrics(bgB(baseB,:),fuB(baseB,:),gridObs(baseB,:));
covStruct = struct('cells',sum(doOI),'cellsBase',sum(baseMask), ...
    'fracDaysFiniteBase', mean(isfinite(fuB(baseB,:)),'all'));

cvIdx = find(baseMask); cvIdx = cvIdx(bPosMap(cvIdx) > 0);   % scored cells bear a station
nCv = numel(cvIdx);
cvBearing = bPosMap(cvIdx);
cvCapIdx = capIdx(cvIdx,:); cvCapDist = capDist(cvIdx,:);
cvErrO = errO;   % full bearing vector, indexed by bearing position

out = fullfile(P.work,sprintf('sensoi_%s.mat',tag));
if keepFusion
    save(out,'doOI','dIn','dFix','covStruct','V','applied','cvIdx','cvBearing', ...
         'cvCapIdx','cvCapDist','cvErrO','nCv','fusion','-v7.3');
else
    save(out,'doOI','dIn','dFix','covStruct','V','applied','cvIdx','cvBearing', ...
         'cvCapIdx','cvCapDist','cvErrO','nCv','-v7.3');
end
fprintf(['[%s] OI %d cells (base %d) | own-set MAE %.4f r %.4f | ', ...
         'fixed-set MAE %.4f r %.4f | CV bundle %d cells\n'], ...
    tag, sum(doOI), sum(baseMask), dIn.mae, dIn.corr, dFix.mae, dFix.corr, nCv);
warning(w0);
end

function D = evalmetrics(bgB, fuB, obs)
ok = isfinite(obs) & isfinite(bgB) & isfinite(fuB);
o = obs(ok); b = bgB(ok); f = fuB(ok);
D = struct();
D.nPairs = numel(o);
D.maeB  = mean(abs(b-o)); D.mae  = mean(abs(f-o));
D.rmseB = sqrt(mean((b-o).^2)); D.rmse = sqrt(mean((f-o).^2));
D.biasB = mean(b-o);  D.bias = mean(f-o);
D.corrB = pcorr(b,o); D.corr = pcorr(f,o);
D.fracBetter = mean(abs(f-o) < abs(b-o));
end

function c = pcorr(x,y)
x = x(:) - mean(x); y = y(:) - mean(y);
c = sum(x.*y)/sqrt(sum(x.^2)*sum(y.^2));
end

function d = hav(lat1,lon1,lat2,lon2)
p1 = lat1*pi/180; p2 = lat2*pi/180;
dphi = (lat2-lat1)*pi/180; dl = (lon2-lon1)*pi/180;
a = sin(dphi/2).^2 + cos(p1).*cos(p2).*sin(dl/2).^2;
d = (2*asin(sqrt(a)))*180/pi;
end
