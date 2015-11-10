function sino = phaserLargeDetectorBinning( sino )

if ndims(sino) ~= 3
    error('only works for 3D geomerty. \n');
end

nv = size( sino, 1 );

if mod( nv , 6 ) ~= 0
    error('incorrect detector height. \n');
end

k = nv / 6;

for i = 0 : k-1
   
    row = ( sino( 2*i+1,:,: ) + sino( 2*i+2,:,: ) ) / 2;
    
    sino( 2*i+1,:,: ) = row;
    sino( 2*i+2,:,: ) = row;
    
    
    row = ( sino( nv-(2*i),:,: ) + sino( nv-(2*i+1),:,: ) ) / 2;
    
    sino( nv-(2*i),:,: )    = row;
    sino( nv-(2*i+1),:,: )  = row;
    
    
end









end