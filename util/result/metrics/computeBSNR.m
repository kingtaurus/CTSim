function bsnr = computeBSNR( img, innerBoundary, outerBoundary, backgroundMap )

sigma = computeTV( img, backgroundMap );

imgFeature= img( innerBoundary ); 
mu1 = mean(  imgFeature(:) );
imgFeature= img( outerBoundary ); 
mu2 = mean(  imgFeature(:) );

bsnr = abs( mu1 - mu2 ) ;

fprintf('The bsnr is %2.5f. \n', bsnr );


end