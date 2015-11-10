function error = computeAbsError( img1, img2, map )

d = abs( img1(map(:)) - img2(map(:)) );

error = mean( d );

fprintf('The abs error is %3i HU. \n', round(error) );


end