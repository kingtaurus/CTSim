function sino = subtractSinoOffset( sino )

noViews = size( sino, 3 );


offsets = zeros( 1, noViews );

for iv = 1 : noViews
    
    leftSide = mean( sino(:, 1, iv ) );
    rightSide = mean( sino(:, end, iv ) );
    
    
    if abs( leftSide ) < abs( rightSide )
        
        offsets( iv ) =  leftSide;
    else
        offsets( iv ) =  rightSide;
    end
    
end

offsets = medfilt1( offsets , 20 );
for iv = 1 : noViews
    
    sino(:, :, iv ) =  sino(:, :, iv ) - offsets( iv );
    
end

end