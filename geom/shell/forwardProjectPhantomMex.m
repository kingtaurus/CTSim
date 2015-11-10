function sino = forwardProjectPhantomMex( img, phan, geom, projType )
% Forward project phantom
% input:
%       img - image to project
%       phan - image parameters including spacing and size
%       geom - system geometry
%       projType
% output:
%       sinp - projection result aka sinogram
%
% Meng Wu
% 2013.12
if nargin < 4
    projType = 'proj,tf';
end

geom.reconSpacing = phan.spacing;
geom.reconSize = phan.size;
geom.reconOffset = phan.offset;

sino = forwardProjectMex( img, geom, 1, 0, projType, true( phan.size(1:2) ) );

end