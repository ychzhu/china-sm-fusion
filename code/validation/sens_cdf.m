function sens_cdf(tag, genR, genK)
%SENS_CDF  Stage A: neighbour selection + CDF bias correction, one variant.
%   Selection, eligibility and the CDF reference pool follow the release exactly
%   (m5/m6): a +/-radius box with a cap of K nearest cells, >=3 neighbours and
%   >minHist valid background days; the reference is the pooled grid_obs of the
%   (present) neighbours inside the box, whether or not those cells bear a
%   station. Variants change only the parameters.
if nargin < 2 || isempty(genR), genR = 1.62; end
if nargin < 3 || isempty(genK), genK = 24; end
w0 = warning('off','all');
P = paths_full();
S  = load(fullfile(P.work,sprintf('sens_geo_R%dK%d.mat',round(genR*100),genK)));
Sg = load(fullfile(P.work,'grid_prep_full.mat'));
V  = variant(tag);

bearing = Sg.bearing; gridObs = double(Sg.gridObs); smapBg = double(Sg.smapBg);
glat = double(Sg.glat); glon = double(Sg.glon);
nGrid = numel(glat); nEff = numel(bearing);
blat = glat(bearing); blon = glon(bearing);
bPosMap = -ones(nGrid,1); bPosMap(bearing) = 1:nEff;

numAround = zeros(nGrid,1); capIdx = zeros(nGrid,V.maxNear); capDist = nan(nGrid,V.maxNear);
for gi = 1:nGrid
    ci = S.gIdx(gi,:); ci = ci(ci>0); cd = S.gDist(gi,1:numel(ci));
    if isempty(ci), continue; end
    keep = abs(blat(ci)-glat(gi)) <= V.radius & abs(blon(ci)-glon(gi)) <= V.radius;
    ci = ci(keep); cd = cd(keep);
    numAround(gi) = numel(ci);
    m = min(V.maxNear,numel(ci));
    capIdx(gi,1:m) = ci(1:m); capDist(gi,1:m) = cd(1:m);
end

daysValid = sum(isfinite(smapBg),1);
smapCdf = smapBg;                       % double; released as single later
applied = 0;
for gi = 1:nGrid
    if numAround(gi) < V.minNbr || daysValid(gi) <= V.minHist, continue; end
    p = capIdx(gi,:); p = p(p>0);   % already positions in the bearing list
    if isempty(p), continue; end
    ref = reshape(gridObs(p,:).',1,[]);
    smapCdf(:,gi) = cdf_cell(smapBg(:,gi), ref, V.polyOrder);
    applied = applied + 1;
end

out = fullfile(P.work,sprintf('senscdf_%s.mat',tag));
save(out,'smapCdf','capIdx','capDist','numAround','applied','tag','V','-v7.3');
fprintf('[%s] R%.2f K%d H%d P%d: CDF applied %d grids\n', ...
    tag, V.radius, V.maxNear, V.minHist, V.polyOrder, applied);
warning(w0);
end

function V = variant(tag)
V = struct('radius',1.08,'maxNear',9,'a',10,'smapErr',0.046, ...
           'minHist',30,'minNbr',3,'polyOrder',2);
switch tag
    case 'R070',  V.radius = 0.70;
    case 'R090',  V.radius = 0.90;
    case 'R126',  V.radius = 1.26;
    case 'R144',  V.radius = 1.44;
    case 'K03',   V.maxNear = 3;
    case 'K06',   V.maxNear = 6;
    case 'K15',   V.maxNear = 15;
    case 'K24',   V.maxNear = 24;
    case 'A03',   V.a = 3;
    case 'A05',   V.a = 5;
    case 'A20',   V.a = 20;
    case 'A30',   V.a = 30;
    case 'SE040', V.smapErr = 0.040;
    case 'SE052', V.smapErr = 0.052;
    case 'SE060', V.smapErr = 0.060;
    case 'H090',  V.minHist = 90;
    case 'H180',  V.minHist = 180;
    case 'P1',    V.polyOrder = 1;
    case 'P3',    V.polyOrder = 3;
    case 'M05',   V.minNbr = 5;
    case 'BASE',
    otherwise, error('unknown tag %s', tag);
end
end
