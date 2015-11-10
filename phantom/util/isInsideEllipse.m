function inside = isInsideEllipse(point, center, semiaxes, rotation)
% inside = isInsideEllipse(point, center, semiaxes, rotation)
%
% INPUT:
% point      N-dim. coordinates of point to check whether it is inside the
%            given ellipse
% center     N-dim. coordinates of ellipse center
% semiaxes   N-dim. vector of lengths of the ellipse's semi-axes
% rotation   N-dim. rotation matrix, specifying the rotation of the ellipse
%            w.r.t. the coordinate sytem
%
% OUTPUT:
% inside     logical return value, specifying whether the given point is
%            inside the given ellipse or not
%
% Copyright (c) 2012 by Andreas Keil, Stanford University.


%% Check Input Arguments

if nargin < 3
	error('Not enough arguments! Specify at least the point and center coordinates as well as the semi-axes radii.');
end

validateattributes(point, {'numeric'}, {'vector', 'integer'});
point = point(:);
noDims = length(point);

validateattributes(center, {'numeric'}, {'vector', 'real'});
if length(center) ~= noDims
	error('Dimension mismatch!');
end
center = center(:);

validateattributes(semiaxes, {'numeric'}, {'vector', 'real', 'positive'});
% if only one semi-axis is given, assume a circle
if isscalar(semiaxes) == 1
	semiaxes = semiaxes*ones(noDims, 1);
end
if length(semiaxes) ~= noDims
	error('Dimension mismatch!');
end
semiaxes = semiaxes(:);

if nargin < 4
	rotation = eye(noDims);
end
validateattributes(rotation, {'numeric'}, {'size', [noDims, noDims], 'real'});
if abs(det(rotation)-1) > 100*eps
	error('Given matrix is not a rotation matrix!');
end


%% Compute Return Value

inside = sum(((rotation'*(point-center))./semiaxes).^2) < 1;
