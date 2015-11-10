function img = reconTwoPassConebeamArtifactCorrectionBone( sinoAtt, geom, reconFunc, muSoft )

img = feval( reconFunc, sinoAtt );

imgBone = img - 1.5 * muSoft; % clusterThreshold(2);
imgBone( imgBone < 0 ) = 0;

% compuate the sinogram error
sinoSoft = sinoAtt - forwardProjectMex( imgBone, geom );

% add reconstructed error
imgSoft = feval( reconFunc, sinoSoft );
img = imgBone + imgSoft;

end

