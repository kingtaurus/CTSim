function J = upsampleVolumes( I, r )

J = zeros( (size(I)-1)*r , class(I) );

[x, y] =  meshgrid( (0:size(I,2)-1)*r, (0:size(I,1)-1)*r );
[vx, vy] =  meshgrid( (0:size(J,2)-1), (0:size(J,1)-1) );

h = fspecial('Gaussian', [3 3]*r , 0.2 * r  );

for i = 1:size( I, 3 )-1
    
    slice1 =  interp2(x, y, squeeze( I(:,:,i)), vx,vy, 'linear');
    if i == size(I,3)
        slice2 =  slice1;
    else
        slice2 =  interp2( x, y,  squeeze( I(:,:,i+1)),vx,vy, 'linear');
    end
    
    for j = 0:r-1
        J( :,:, (i-1)*r + j + 1 ) = medfilt2( ( (r-j) / r ) * slice1 + ( j / r ) * slice2, [1 1]) ;
    end
    
end

for j = 1 : size( J, 1 )
   J(j,:,:) = imfilter( squeeze(J(j,:,:)), h ); 
end

for j = 1 : size( J, 2 )
   J(:,j,:) = imfilter( squeeze(J(:,j,:)), h ); 
end

for j = 1 : size( J, 3 )
   J(:,:,j) = imfilter( squeeze(J(:,:,j)), h ); 
end

end