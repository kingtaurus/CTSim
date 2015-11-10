function img1 = updateImageAxialBlocks( img1, img2, map, type )


for iz = 1: size( img1 , 3)
    slice1 = img1(:, :, iz);
    slice2 = img2(:, :, iz);
    
    if type
        slice1( map) = slice2( map );
    else
        slice1( ~map ) = slice2( ~map );
    end
    
     img1(:, :, iz) = slice1;
    
    
end