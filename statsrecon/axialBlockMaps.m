function maps = axialBlockMaps( geom, numAS  )
% function maps = axialBlockMaps( geom, numAS  )
% generate the chess board map for axial block updates 
%
% Meng Wu at stanford University
% 2014

n = sqrt( numAS );
if n ~= round( n )
    error('No. of subsets must be square of integer.\n');
end

maps = cell( numAS, 1);

for i = 1 : numAS
    
    map = false( geom.reconSize(1), geom.reconSize(2) );
    
    j = mod( i - 1, n ) ;
    k = floor( (i -1 )  / n ) ;
    
    J = floor( geom.reconSize(1) / n  );
    K = floor( geom.reconSize(2) / n  );
    
    
    
    map( J*j+1:J*(j+1) ,  K*k+1:K*(k+1) ) = true;
    
    map = map & geom.map;
    
    maps{i} = map;
    
    
end
    


end