function m5_neighbors()
%M5_NEIGHBORS  Capped neighbour selection (R1=1.08, MAX_NEAR=9).
%   Faithful port of s5_neighbors.py. cap_idx holds 1-based bearing positions
%   (0 = none), cap_dist the haversine distances in degrees.

P = paths_full();
w0=warning('off','all');

R1=1.08; MAX_NEAR=9;

Sg=load(fullfile(P.work,'grid_prep_full.mat'));
bearing=Sg.bearing; glat=Sg.glat; glon=Sg.glon;
blat=glat(bearing); blon=glon(bearing);
nGrid=numel(glat);

numAround=zeros(nGrid,1);
capIdx=zeros(nGrid,MAX_NEAR);
capDist=nan(nGrid,MAX_NEAR);
for gi=1:nGrid
    cand=find(abs(blat-glat(gi))<=R1 & abs(blon-glon(gi))<=R1);
    numAround(gi)=numel(cand);
    if isempty(cand); continue; end
    dist=haversine_deg(glat(gi),glon(gi),blat(cand),blon(cand));
    % deterministic tie-break: distance first, then bearing position
    [~,ord]=sortrows([dist(:),cand(:)],[1,2]);
    k=ord(1:min(MAX_NEAR,numel(ord)));
    keep=cand(k);
    capIdx(gi,1:numel(keep))=keep;
    capDist(gi,1:numel(keep))=dist(k);
end

fprintf('grids with >=3 neighbors: %d\n',sum(numAround>=3));
fprintf('num_around min=%d median=%d max=%d\n',min(numAround),median(numAround),max(numAround));
save(fullfile(P.work,'neighbors_full.mat'),'numAround','capIdx','capDist','-v7.3');
fprintf('saved neighbors_full.mat\n');
warning(w0);
end

function d=haversine_deg(lat1,lon1,lat2,lon2)
p1=radians(lat1); p2=radians(lat2);
dphi=radians(lat2-lat1); dl=radians(lon2-lon1);
a=sin(dphi/2).^2+cos(p1).*cos(p2).*sin(dl/2).^2;
d=degrees(2*asin(sqrt(a)));
end

function r=radians(x); r=x*pi/180; end
function d=degrees(x); d=x*180/pi; end
