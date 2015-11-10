function sino = cosineRayWeighting( sino, geom, flag )
% function sino = cosineWeighting( sino, geom )
%   cosine weighting for cone beam geometry
%
% Meng Wu, 2014

if nargin < 3
    flag = 1;
end

nu = geom.detSize(1);
nv = geom.detSize(2);
du = geom.detSpacing(1);
dv = geom.detSpacing(2);
ou = geom.detOffset(1);
ov = geom.detOffset(2);

noViews = geom.noViews;
SDD = geom.SDD;

%sino = binningDetecotRows( sino, round(  geom.reconSpacing(3) / geom.detSpacing(2)  ) );

% array of detector positions in mm
u =  ( ( -(nu-1)/2:(nu-1)/2) + ou )  * du;
v = ( ( -(nv-1)/2:(nv-1)/2) + ov ) * dv ;

[uu, vv ] = meshgrid( u, v);

if geom.flatPanel
    weight = sqrt( SDD^2  + vv.^2 ) ./ sqrt(  uu.^2 + vv.^2 + SDD^2 );
else
    weight = SDD ./ sqrt( SDD.^2 + vv.^2 ) ;
end

if ~flag
    weight = 1./ weight;
end

for iview = 1:noViews
    % get a view and apply weighting
    sino(:, :,  iview) = weight .* sino(:, :, iview);
end