function sino = cosineWeighting( sino, geom )
% function sino = cosineWeighting( sino, geom )
%   cosine weighting for cone beam geometry
%
% Meng Wu, 2014

nu = geom.detSize(1);
du = geom.detSpacing(1);
ou = geom.detOffset(1);
noViews = geom.noViews;
SDD = geom.SDD;
SAD = geom.SAD;

if ndims( sino ) == 3
    
    nv = geom.detSize(2);
    dv = geom.detSpacing(2);
    ov = geom.detOffset(2);
    
    % array of detector positions in mm
    u =  ( ( -(nu-1)/2:(nu-1)/2) + ou )  * du;
    v = ( ( -(nv-1)/2:(nv-1)/2) + ov ) * dv ;
    
    [uu, vv ] = meshgrid( u, v);
    
    if geom.flatPanel
        weight = SAD * sqrt( 1.0 + ( uu.^2 / SDD^2 ) ) ./ sqrt(  uu.^2 + vv.^2 + SDD^2 );
    else
        weight = SAD / SDD * cos( uu ./ ( SDD * sqrt( 1 + (vv / SDD ).^2 )) );
    end
    
    for iview = 1:noViews
        % get a view and apply weighting
        sino(:, :,  iview) = weight .* sino(:, :, iview);
    end
    
else
    
    % array of detector positions in mm
    u = (( -(nu-1)/2:(nu-1)/2)+ geom.detOffset )' * du  ;
    
    % cosine weighting factors and filter
    if geom.flatPanel
        weight  = SAD ./ sqrt( u.^2 + SDD^2 );
    else
        weight  = SAD / SDD * cos( u / SDD );
    end
    
    
    for iview = 1:noViews
        % get a view and apply weighting
        sino(:, iview) = weight .* sino(:, iview);
    end
    
    
    
end

