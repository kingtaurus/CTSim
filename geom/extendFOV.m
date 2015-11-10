function img = extendFOV( img, geom )

if ~geom.helicalScan && length( geom.detSize ) == 2 && length( geom.reconSize ) == 3
    
    k = round( ( geom.reconSize(3) - geom.detSize( 2 ) * geom.detSpacing( 2 ) / 2 / geom.reconSpacing(3) ) / 2);
    
    k = max(k , 0);
    
    img = extendVoi( img, k );
    
end


end