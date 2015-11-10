function sigma = computeStd( img, map )


imgTissue = img( map ); 

sigma = std( imgTissue(:) );

%fprintf('The std is %3i HU. \n', round(sigma) );

end