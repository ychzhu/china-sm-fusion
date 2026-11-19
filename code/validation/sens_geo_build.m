function sens_geo_build(genR, genK)
%SENS_GEO_BUILD  Neighbour superset over the +/-genR box, sorted by distance.
%   The published chain selects neighbours by a +/-1.08 degree box (m5). Radius
%   variants therefore re-box the same superset; the haversine distance is used
%   only for ordering and for the correlation model.
if nargin < 1 || isempty(genR), genR = 1.62; end
if nargin < 2 || isempty(genK), genK = 24; end
P = paths_full();
Sg = load(fullfile(P.work,'grid_prep_full.mat'));
bearing = Sg.bearing; glat = double(Sg.glat); glon = double(Sg.glon);
nGrid = numel(glat);
blat = glat(bearing); blon = glon(bearing);

gNum = zeros(nGrid,1); gIdx = zeros(nGrid,genK); gDist = nan(nGrid,genK);
for gi = 1:nGrid
    cand = find(abs(blat-glat(gi)) <= genR & abs(blon-glon(gi)) <= genR);
    gNum(gi) = numel(cand);
    if isempty(cand), continue; end
    d = hav(glat(gi),glon(gi),blat(cand),blon(cand));
    [~,ord] = sortrows([d(:),cand(:)],[1,2]);
    kk = cand(ord); dd = d(ord);
    m = min(genK,numel(kk));
    gIdx(gi,1:m) = kk(1:m); gDist(gi,1:m) = dd(1:m);
end
S = struct('gNum',gNum,'gIdx',gIdx,'gDist',gDist,'genR',genR,'genK',genK);
save(fullfile(P.work,sprintf('sens_geo_R%dK%d.mat',round(genR*100),genK)),'-struct','S','-v7.3');
% sanity: published box
nb = zeros(nGrid,1);
for gi = 1:nGrid
    ci = gIdx(gi,:); ci = ci(ci>0);
    nb(gi) = min(sum(abs(blat(ci)-glat(gi))<=1.08 & abs(blon(ci)-glon(gi))<=1.08), 9);
end
fprintf('geo saved (genR=%.2f, K=%d). published box: >=3 nbrs = %d\n', genR, genK, sum(nb>=3));
end

function d = hav(lat1,lon1,lat2,lon2)
p1 = lat1*pi/180; p2 = lat2*pi/180;
dphi = (lat2-lat1)*pi/180; dl = (lon2-lon1)*pi/180;
a = sin(dphi/2).^2 + cos(p1).*cos(p2).*sin(dl/2).^2;
d = (2*asin(sqrt(a)))*180/pi;
end
