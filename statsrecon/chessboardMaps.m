function maps = chessboardMaps( geom, numAS  )

n = sqrt( numAS );
if n ~= round( n )
    error('No. of subsets must be square of integer.\n');
end

maps = cell( numAS, 1);

for i = 1 : numAS
    
    map = false( geom.reconSize(1), geom.reconSize(2) );
    
    j = mod( i - 1, n ) + 1;
    k = floor( (i -1 )  / n ) + 1;
    
    map(j:n:end, k:n:end ) = true;
    
    map = map & geom.map;
    
    maps{i} = map;
    
    
end
    


end