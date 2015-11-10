load 'temp.mat';

[phan, map] = loadXCATPhantom(p);

[ geom ] = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e6);

[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);

%% Normal X-ray CT

% compute sinogram
sinoIntensity= simulatePhotonCountingData( phan, geom, spectrum, sinosDirKeV );
sinoAtt = computePhotonCountingSinogramAttunation( sinoIntensity, spectrum );

% corrections
sinoAtt = beamHardeningWarterCorrection(sinoAtt, spectrum);
sinoAtt = medianFilterSino( sinoAtt, 3 );

%% reconstruction

imgAttFBP = reconFBP( sinoAtt, geom, 'hamming' );
imgFBP = convertMonoAttToHu( imgAttFBP, spectrum) ;
figure; imshow( imgFBP, map.windowHu ); colormap gray;

%% PCXD

[ sinoPhotonCounts, energyBinSizes ] = simulatePCXDData( phan, geom, spectrum, sinosDirKeV, [0 50 80], true );

sinoAtt = computePCXDSinogramAttunation( sinoPhotonCounts, energyBinSizes );

%% reconstruction

imgAttFBP = reconFBP( sinoAtt{1}, geom, 'hamming' );
figure; imshow( imgAttFBP, map.windowAtt ); colormap gray;