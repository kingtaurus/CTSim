function [img1, img2] = backProject2Mex( sino1, sino2, geom, M, i, projType )
% Backward projection using distance-driven method with ordered-subset in 2D
%       and 3d. Exact traspose to the forwardProjectDistanceDrivenMex
% input:
%       sino - projection result aka sinogram
%       geom - system geometry
%       M    - number of subsets
%       i    - subset number
% output:
%       img - image to project
%
% Meng Wu
% 2013.4

if nargin < 3
    M = 1;
    i = 0;
end

if nargin < 5
    projType = 'back,dd';
end

img1 = backProjectMex( sino1, geom, M, i, projType );
img2 = backProjectMex( sino2, geom, M, i, projType );
