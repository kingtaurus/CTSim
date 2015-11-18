load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadXCATPhantom(p);

[geom ] = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e6);

% Compute ground trut
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);

%% simulate projection data
sinoAtt = simulateAttenuationDataPhantom( phan, geom, spectrum, sinosDirKeV);

% compute sinogram
% sinoPhotonCounts = simulatePhotonCountingData( phan, geom, spectrum, sinosDirKeV, true );
% sinoAtt = computePhotonCountingSinogramAttunation( sinoPhotonCounts, spectrum );
% % corrections
% sinoAtt = beamHardeningWarterCorrection(sinoAtt, spectrum);
% sinoAtt = medianFilterSino( sinoAtt, 3 );

sinoAtt0 = sinoAtt;
geom0 = geom;

%% Trunctation correction

viewUpsamplingRate = 1;
truncationCorrection = true;
softTissueAtt = 0.15;
if truncationCorrection
    [ sinoAtt, geom ] = truncationCorrectionEllipicalPatch( sinoAtt0, geom0, softTissueAtt, 2, 1, 64 );
end


%%

imgAtract1 = reconATRACT1( sinoAtt0, geom0, 'hamming' );
figure; imagesc( imgAtract1(:,:,end/2)); colormap gray;

imgAtract2 = reconATRACT2( sinoAtt0, geom0, 'hamming' );
figure; imagesc( imgAtract2(:,:,end/2)); colormap gray;

imgAtractm = reconATRACTm( sinoAtt0, geom0, 'hamming' );
figure; imagesc( imgAtractm(:,:,end/2)); colormap gray;

imgFBP = reconFBP( sinoAtt0, geom0, 'hamming' );
figure; imagesc( imgFBP(:,:,end/2)); colormap gray;


%% full scan FBP reconstrunction

imgAttFBP = reconFBP( sinoAtt, geom, 'hamming', viewUpsamplingRate);
imgFBP = convertMonoAttToHu( imgAttFBP, spectrum) ;

figure; imshow( imgFBP(:,:,end/2), map.windowHu ); colormap gray;
% export_fig( 'FBP-Varian-60-center.jpg' );
% 
% figure; imshow( imgFBP(:,:,end), map.windowHu ); colormap gray;
% export_fig( 'FBP-Varian-60-top.jpg' );
% 
% figure; imshow( imgFBP(:,:,1), map.windowHu ); colormap gray;
% export_fig( 'FBP-Varian-60-bottom.jpg' );
% 
% figure; imdisp( fliplr( squeeze( imgFBP(end/2,1:2:end,:) ) ), map.windowHu ); colormap gray;
% export_fig( 'FBP-Varian-60-coronal.jpg' );

%% short scan

offsetAngle = 0;
[geomShort, sinoShort ]  = convertSinogram2ShortScan( geom, sinoAtt, offsetAngle );

%imgAttShort = reconFBP( sinoShort, geomShort, 'hamming', 1);
imgAttShort = reconShiftInvariantFBPShortScan( sinoShort, geomShort, 'hamming', 1);


imgFBPShort = convertMonoAttToHu( imgAttShort, spectrum) ;
figure; imshow( imgFBPShort(:,:,end/2), map.windowHu ); colormap gray;


% figure; imshow( imgFBPShort(:,:,end/2), map.windowHu ); colormap gray;
% export_fig( 'FBP-Varian-60-short-center.jpg' );
% 
% figure; imshow( imgFBPShort(:,:,end), map.windowHu ); colormap gray;
% export_fig( 'FBP-Varian-60-short-top.jpg' );
% 
% figure; imshow( imgFBPShort(:,:,1), map.windowHu ); colormap gray;
% export_fig( 'FBP-Varian-60-short-bottom.jpg' );
% 
% figure; imdisp( fliplr( squeeze( imgFBPShort(end/2,1:2:end,:) ) ), map.windowHu ); colormap gray;
% export_fig( 'FBP-Varian-60-short-coronal.jpg' );

%% 2D test case

[geom2d, sino2d ]  = convertSinogram2TowDimensional( geom0, sinoAtt0 );

if truncationCorrection
    [ sino2d, geom2d ] = truncationCorrectionEllipicalPatch( sino2d, geom2d, softTissueAtt, 3, 1, 64 );
end

imgAtt2d = reconFBP( sino2d, geom2d, 'hamming' );
imgFBP2d = convertMonoAttToHu( imgAtt2d, spectrum) ;
figure; imshow( imgFBP2d, map.windowHu ); colormap gray;

%% 2D test case

offsetAngle = 0;
[geomShort, sinoShort ]  = convertSinogram2ShortScan( geom0, sinoAtt0, offsetAngle );
[geomShort2d, sinoShort2d ]  = convertSinogram2TowDimensional( geomShort, sinoShort );

if truncationCorrection
    [ sinoShort2d, geomShort2d ] = truncationCorrectionEllipicalPatch( sinoShort2d, geomShort2d, softTissueAtt, 3, 1, 64 );
end

imgAttShort = reconFBP( sinoShort2d, geomShort2d, 'hamming' );
imgFBPShort = convertMonoAttToHu( imgAttShort, spectrum) ;
figure; imshow( imgFBPShort, map.windowHu ); colormap gray;
