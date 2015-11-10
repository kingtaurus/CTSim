function evaluateSoftTissue( img, spectrum, map )

img = convertMonoAttToHu( img, spectrum);

fprintf( [ inputname(1) ': \n\t'] );

computeSoftTissueTV( img, map );

fprintf('\t');
computeSoftTissueStd( img, map );

fprintf( '\n\n' );

end