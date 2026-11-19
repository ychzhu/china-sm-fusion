function sens_prep(genR, genK)
%SENS_PREP  Verify the cached superset and define the fixed scoring domain.
if nargin < 1 || isempty(genR), genR = 1.62; end
if nargin < 2 || isempty(genK), genK = 24; end
P = paths_full();
geoFile = fullfile(P.work,sprintf('sens_geo_R%dK%d.mat',round(genR*100),genK));
if ~isfile(geoFile), sens_geo_build(genR,genK); end

Sg = load(fullfile(P.work,'grid_prep_full.mat'));
S  = load(geoFile);
bearing = Sg.bearing; glat = double(Sg.glat); glon = double(Sg.glon);
blat = glat(bearing); blon = glon(bearing);
smapBg = double(Sg.smapBg); nGrid = numel(glat);
daysValid = sum(isfinite(smapBg),1);

numBase = zeros(nGrid,1);
for gi = 1:nGrid
    ci = S.gIdx(gi,:); ci = ci(ci>0);
    numBase(gi) = min(sum(abs(blat(ci)-glat(gi))<=1.08 & abs(blon(ci)-glon(gi))<=1.08), 9);
end
baseMask = (numBase >= 3) & (daysValid(:) > 30);
save(fullfile(P.work,'sens_base_domain.mat'),'baseMask');
fprintf('fixed scoring domain (published settings): %d cells\n', sum(baseMask));
end
