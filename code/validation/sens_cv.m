function sens_cv(tag, dayStride)
%SENS_CV  Stage C: strict station-grouped 10-fold CV on the fixed domain.
%   Folds are fixed to the archived permutation (cv_fold_perm.csv) over the 1932
%   station-bearing cells, so every variant is scored on the same held-out
%   cells. CDF functions are inherited from the variant's stage A; only the OI
%   weights are refitted per fold from training neighbours.
if nargin < 2 || isempty(dayStride), dayStride = 3; end
w0 = warning('off','all');
P = paths_full();

Sc = load(fullfile(P.work,sprintf('senscdf_%s.mat',tag)));
V = Sc.V;
Sg = load(fullfile(P.work,'grid_prep_full.mat'));
bearing = Sg.bearing; gridObs = double(Sg.gridObs); glat = double(Sg.glat); glon = double(Sg.glon);
smapCdf = Sc.smapCdf;
Si = load(fullfile(P.work,sprintf('sensoi_%s.mat',tag)));
cvBearing = Si.cvBearing; cvCapIdx = Si.cvCapIdx; cvCapDist = Si.cvCapDist;
cvErrO = Si.cvErrO; nCv = Si.nCv;

blat = glat(bearing); blon = glon(bearing);
nEff = numel(bearing); nDays = size(gridObs,2);

permT = readtable(fullfile(P.here,'cv_fold_perm.csv'));
perm1 = permT.perm0 + 1;
NF = 10; edges = round(linspace(0,nEff,NF+1));
foldOf = zeros(nEff,1);
for k = 1:NF
    m = perm1(edges(k)+1:edges(k+1)); foldOf(m) = k;
end

dsel = 1:dayStride:nDays; nds = numel(dsel);
Bv = zeros(nCv*nds,1); Fv = Bv; Ov = Bv; n = 0;

for r = 1:nCv
    v = cvBearing(r);                     % bearing position of the validation cell
    k = foldOf(v);
    ci = cvCapIdx(r,:); ci = ci(ci>0);
    cd = cvCapDist(r,1:numel(ci));
    tr = foldOf(ci) ~= k;
    nb = ci(tr); nd = cd(tr);
    b = smapCdf(dsel, bearing(v));
    o = gridObs(v, dsel).';
    if numel(nb) < V.minNbr
        f = b;
    else
        nA = numel(nb); la = blat(nb); lo = blon(nb); dm = zeros(nA,nA);
        for q = 1:nA, dm(q,:) = hav(la(q),lo(q),la,lo); end
        B = 1.0./(V.a.^dm); corIk = 1.0./(V.a.^nd);
        errCon = cvErrO(nb)/V.smapErr;
        bad = ~isfinite(errCon) | errCon < 0; errCon(bad) = 1.0;
        w = pinv(B + diag(errCon.^2)) * corIk(:);
        rg = gridObs(nb, dsel).' - smapCdf(dsel, bearing(nb));
        rg(~isfinite(rg)) = 0;
        f = b + rg * w;
    end
    ok = isfinite(o) & isfinite(b) & isfinite(f);
    idx = n + (1:sum(ok));
    Bv(idx) = b(ok); Fv(idx) = f(ok); Ov(idx) = o(ok);
    n = n + sum(ok);
end
Bv = Bv(1:n); Fv = Fv(1:n); Ov = Ov(1:n);

CV = struct('tag',tag,'dayStride',dayStride,'nPairs',n,'nCells',nCv);
CV.maeB = mean(abs(Bv-Ov));  CV.maeF = mean(abs(Fv-Ov));
CV.rmseB = sqrt(mean((Bv-Ov).^2)); CV.rmseF = sqrt(mean((Fv-Ov).^2));
CV.corrB = pcorr(Bv,Ov); CV.corrF = pcorr(Fv,Ov);
CV.biasB = mean(Bv-Ov); CV.biasF = mean(Fv-Ov);
CV.fracCloser = mean(abs(Fv-Ov) < abs(Bv-Ov));
CV.impMae  = 100*(1-CV.maeF/CV.maeB);
CV.impCorr = CV.corrF - CV.corrB;

save(fullfile(P.work,sprintf('senscv_%s.mat',tag)),'-struct','CV');
fprintf(['[%s] CV n=%d cells=%d stride=%d | MAE %.5f -> %.5f (%.1f%%) | ', ...
         'RMSE %.5f -> %.5f (%.1f%%) | r %.4f -> %.4f | closer %.1f%%\n'], ...
    tag, n, nCv, dayStride, CV.maeB, CV.maeF, CV.impMae, ...
    CV.rmseB, CV.rmseF, 100*(1-CV.rmseF/CV.rmseB), CV.corrB, CV.corrF, 100*CV.fracCloser);
warning(w0);
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
