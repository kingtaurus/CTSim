function [sinoAtt, sinoPC] = loadLiverCTSinogram( fileName, geom, dataSize, fileType )
%  [sinoAtt, sinoPC] = loadLiverCTSinogram( fileName, geom, dataSize, fileType )
%
% Meng Wu at Stanford University
% 2013


if nargin < 3
    dataSize = [ geom.detSize(1) geom.detSize(2) geom.noViews ]; 
end

if nargin < 4
    fileType = 'single';
end


dataPath = 'C:\Users\febme_000\Documents\MATLAB\kvmv\data\liver-ct-data\';

fprintf('Load liver CT data from %s ...', [dataPath fileName]);

fid = fopen( [dataPath fileName, '.prep' ], 'r' ); 
a = fread(fid, inf, fileType); 
fclose( fid );
a = reshape( a, dataSize );

sinoAtt = rotateSinogram( a, 3 );

clear a;

fid = fopen( [dataPath fileName, '.scan' ], 'r' ); 
a = fread(fid, inf, fileType); 
fclose( fid );
a = reshape( a, dataSize );

sinoPC = rotateSinogram( a, 3 );

fprintf('done. \n\n');