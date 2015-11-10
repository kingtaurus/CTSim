function sino = filterHelicalSinogram( sino, H, weights, geom, tiltedFilter  )
% function sino = filterHelicalSinogram( sino, H, weights, geom, tiltedFilter  )
%
% Apply ramp filer and helical weighting to the sinogram
%
% Meng Wu
% 2014

u = ( ( -(geom.detSize(1)-1)/2:(geom.detSize(1)-1)/2) + geom.detOffset(1) ) * geom.detSpacing(1);
v = ( ( -(geom.detSize(2)-1)/2:(geom.detSize(2)-1)/2) + geom.detOffset(2) ) * geom.detSpacing(2);
[~, vv ] = meshgrid( u, v);

% compute the fan angle of each detector column
if geom.flatPanel
    gammas   = atan( - u / geom.SDD );
else
    gammas   = - u / geom.SDD ;
end

dbeta = abs(geom.betas(end) - geom.betas(1)) / (geom.noViews-1);
vt   = vv + repmat( gammas , geom.detSize(2), 1) / dbeta * geom.pitch * geom.detSize(2) * geom.detSpacing(2);

%in case not enough projections

if  size( weights, 3) == size( sino, 3 )
    sino = sino .* weights;
else
    ivoffset = floor( ( size( weights, 3) - size( sino, 3 ) ) / 2);
    sino = sino .* weights( :,:, ivoffset + 1: ivoffset + size(sino, 3) );
end


for iview = 1: size( sino, 3 )
    
    proj = sino(:,:,iview);
    
    if tiltedFilter
        proj = tilteProjction(  proj, vv , vt );
        proj =  filterFreq( proj, H, 2);
        sino(:,:,iview) = tilteProjction(  proj, vt , vv );
    else
        sino(:,:,iview) =  filterFreq( proj, H, 2);
    end
end


end