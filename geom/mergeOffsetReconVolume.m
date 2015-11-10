function [ img, geom ] = mergeOffsetReconVolume( img1, img2, geom, offset, na )
% function img = mergeOffsetReconVolume( img1, img2, geom, offset, na )
%   Merge two offest reconstructed volume to a single volume
%   input:
%       img1    - reconstructed volume 1, at the top
%       img2    - reconstructed volume 2, at the bottom
%       geom    - geometry parameters
%       offset  - z offset of the two volumes
%       na      - number of the slices at centers to average
%   output:
%       img     - one big voluem
%       geom    - modified geometry parameters
%
%   Created by Meng Wu, 2014


if nargin < 3
    na = 2;
end

% new size of the image image volume
nz =  round(  geom.reconSize(3) / 2 + offset / 2 / geom.reconSpacing(3) );

img = zeros( [geom.reconSize(1) geom.reconSize(2) nz * 2], 'single'); 
img(:,:,1:nz) = img1( :,:, 1:nz);
img(:,:,end-nz+1:end) = img2( :,:,end-nz+1:end);


% average a few slices
img(:,:,nz-na+1:nz) = ( img(:,:,nz-na+1:nz) + img2( :,:,end-nz+2-na:end-nz+1) )/2;
img(:,:,nz+1:nz+na) = ( img(:,:,nz+1:nz+na) + img1( :,:,nz+1:nz+na) )/2;

geom.reconSize(3) = nz;

end