function [sino1, sino2] = forwardProject2Mex( img1, img2, geom, M, i,  projType )
% Forward projection
% input:
%       img - image to project
%       geom - system geometry
%       M    - number of subsets
%       i    - subset number
% output:
%       sinp - projection result aka sinogram
%
% Meng Wu
% 2011.12

if nargin < 3
    M = 1;
    i = 0;
end

if nargin < 5
    projType = 'proj,dd';
end


sino1 = forwardProjectMex( img1, geom, M, i,  projType );
sino2 = forwardProjectMex( img2, geom, M, i,  projType );

end