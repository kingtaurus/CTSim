function cnr = computeCNR( img, featureMap, backgroundMap )

imgBackground = img( backgroundMap ); 
sigma = std( imgBackground(:) );
mub = mean(  imgBackground(:) );

imgFeature= img( featureMap ); 
muf = mean(  imgFeature(:) );


cnr = abs( muf - mub ) / sigma;

fprintf('The cnr is %2.2f. \n', cnr );


end