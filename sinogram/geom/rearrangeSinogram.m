function [ sinoOut ] = rearrangeSinogram( sinoIn, geom )

[nu, nv, noviews] = size( sinoIn );

if length( geom.detSize ) == 2
    sinoOut = zeros( [nv nu noviews], 'single');
    for i = 1:noviews
        sinoOut(:,:,i) = fliplr( sinoIn(:,:,i)' );
    end
else
    sinoOut = zeros( [nu noviews], 'single');
    for i = 1:noviews
        sinoOut(:,i) = flipud( squeeze( sinoIn( :, round((nv+1)/2),i) ) );
    end
end

end