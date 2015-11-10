function img = reconTwoPassConebeamArtifactCorrectionTissueBone( sinoAtt, geom, reconFunc, muSoft )

img = feval( reconFunc, sinoAtt );

imgBone = img;
imgBone( img > 0.25 * muSoft & img < 1.5 * muSoft ) = 0;
imgSoft = img - imgBone;

% compuate the bone only projection
sinoBone = forwardProjectMex( imgBone, geom );
sinoSoft = forwardProjectMex( imgSoft, geom );

% add reconstructed error
imgArtfBone = feval(  reconFunc, sinoBone ) - imgBone ;
imgArtfSoft = feval(  reconFunc, sinoSoft ) - imgSoft ;

img = img - imgArtfBone - 0.8 * imgArtfSoft;
%img( img < clusterThreshold(2) ) = imgSoft( img < clusterThreshold(2) );




end

