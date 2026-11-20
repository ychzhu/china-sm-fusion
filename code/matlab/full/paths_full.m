function P = paths_full()
%PATHS_FULL  Locations and constants for the MATLAB 3596 three-stage chain.
%
%   Mirrors code/python/full/paths_full.py. All paths are resolved relative to
%   this file. The released data are expected under <repo>/data; set the
%   environment variable SM_DATA_DIR to point at an extracted data archive.

here    = fileparts(mfilename('fullpath'));   % code/matlab/full
matDir  = fileparts(here);                    % code/matlab
codeDir = fileparts(matDir);                  % code
root    = fileparts(codeDir);                 % repository root

dataEnv = getenv('SM_DATA_DIR');
if ~isempty(dataEnv)
    dataRoot = dataEnv;
else
    dataRoot = fullfile(root, 'data');
end

P = struct();
P.here       = here;
P.matDir     = matDir;
P.codeDir    = codeDir;
P.root       = root;
P.data       = dataRoot;
P.station    = fullfile(dataRoot, 'station');
P.smap       = fullfile(dataRoot, 'smap_v09');
P.fusion     = fullfile(dataRoot, 'fusion');
P.reference  = fullfile(dataRoot, 'reference');
P.work       = fullfile(here, 'work_full');
if ~isfolder(P.work); mkdir(P.work); end

P.station3596Nc = fullfile(P.station, 'staD10_qc_3596_20150101_20261001.nc');
P.fusion3596Nc  = fullfile(P.fusion,  'fusion_sm_36km_3596_20150331_20261001.nc');

% V2.1: both fixed fallbacks use the three-decimal TCA median.
P.smapTcaErr      = 0.046;
P.staFallbackErr  = 0.046;
P.sigmaExtreme    = 0.1074;
P.gridRes         = 0.18;
end
