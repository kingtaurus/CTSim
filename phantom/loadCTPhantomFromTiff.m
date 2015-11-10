%function [I, phan] = loadCTPhantomFromTiff( dataPath, dataName, geom )

dataPath = '\data\prostate-patient-ct-fiducial\';
dataName = 'prostateScans';

noSlices = 149;


info = dicominfo([ dataPath 'reference'] );


dx = info.PixelSpacing(1);
dy = info.PixelSpacing(2);
dz = info.SliceThickness;

nx = info.Rows;
ny = info.Columns;
nz = noSlices;


I = zeros( [nx ny nz], 'single' );

for i = 1:nz
    if mod(i,20) == 1
        fprintf('\tslice %d \n', i);
    end
    slice = imread([dataPath dataName '.tif'], 'Index', i);
    
    I(:,:,i) = slice;
end

if min( I(:) ) > -100
    I = I - 1000;
end


meta.NDims = 3;
meta.DimSize = size(I);
meta.ElementSpacing = [ dx, dy, dz];

writeMetaImage(I, [dataName '.mhd'], meta);





