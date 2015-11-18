load 'temp.mat';

systemOffset = 150; 
% Load simulation parameters and datas
[phan, map] = loadXCATPhantom(p);
[geom ] = loadProjectionGeometryCT( p );
spectrum = loadSpectraCT(p, geom, 1e6);

%% compute sinograms

phan.offset = [0 0 -systemOffset/2 / phan.spacing(3)];
sinoDir1 = [sinosDirKeV 'source' num2str(-systemOffset/2) 'mm/'];
sinoPhotonCounts = simulatePhotonCountingData( phan, geom, spectrum, sinoDir1, true );
sinoAtt1 = computePhotonCountingSinogramAttunation( sinoPhotonCounts, spectrum );

phan.offset = [0 0 +systemOffset/2 / phan.spacing(3)];
sinoDir2 = [sinosDirKeV 'source' num2str(systemOffset/2) 'mm/'];
sinoPhotonCounts = simulatePhotonCountingData( phan, geom, spectrum, sinoDir2, true );
sinoAtt2 = computePhotonCountingSinogramAttunation( sinoPhotonCounts, spectrum );

%% corrections
sinoAtt1 = beamHardeningWarterCorrection(sinoAtt1, spectrum);
sinoAtt2 = beamHardeningWarterCorrection(sinoAtt2, spectrum);
sinoAtt1 = medianFilterSino( sinoAtt1, 3 );
sinoAtt2 = medianFilterSino( sinoAtt2, 3 );

% reconstruction parameters
viewUpsamplingRate = 1;
softTissueAtt = 0.15;
window = 'hamming';

%Compute ground truth
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);


%% Full scan reconstruction

imgAttFBP1 = reconFBP( sinoAtt1, geom, window, viewUpsamplingRate);
figure; imshow( convertMonoAttToHu( imgAttFBP1(:,:,end/2), spectrum), map.windowHu ); 

imgAttFBP2 = reconFBP( sinoAtt2, geom, window, viewUpsamplingRate);
figure; imshow( convertMonoAttToHu( imgAttFBP2(:,:,end/2), spectrum), map.windowHu ); 

% combine 2 FDK
imgAttFBP = mergeOffsetReconVolume( imgAttFBP2, imgAttFBP1, geom, systemOffset, 3 );
imgHuFBP = convertMonoAttToHu( imgAttFBP, spectrum) ;
figure; imshow( flipud( squeeze( imgHuFBP(end/2,:,:) )' ), map.windowHu ); colormap gray;
figure; imshow( imgHuFBP(:,:,end/2), map.windowHu ); colormap gray;


% figure; imshow( imgHuFBP(:,:,end/2), map.windowHu ); colormap gray;
% export_fig( 'FBP-DS-center.jpg' );
% 
% figure; imshow( imgHuFBP(:,:,end/2+40), map.windowHu ); colormap gray;
% export_fig( 'FBP-DS-top.jpg' );
% 
% figure; imshow( imgHuFBP(:,:,end/2-40), map.windowHu ); colormap gray;
% export_fig( 'FBP-DS-bottom.jpg' );
% 
% figure; imshow(  flipud( squeeze( imgHuFBP(end/2,:,:) )' ) , map.windowHu ); colormap gray;
% export_fig( 'FBP-DS-coronal.jpg' );


%% Both doing short scans

[geomShort1, sinoShort1 ]  = convertSinogram2ShortScan( geom, sinoAtt1, 0 );
[geomShort2, sinoShort2 ]  = convertSinogram2ShortScan( geom, sinoAtt2, 200 );

imgAttFBPShort1 = reconFBP( sinoShort1, geomShort1, window, viewUpsamplingRate);
imgAttFBPShort2 = reconFBP( sinoShort2, geomShort2, window, viewUpsamplingRate);

imgFBPShort = mergeOffsetReconVolume( imgAttFBPShort2, imgAttFBPShort1, geomShort1, systemOffset, 3 );
imgHuFBPShort = convertMonoAttToHu( imgFBPShort, spectrum) ;
figure; imshow( flipud( squeeze( imgHuFBPShort(end/2,:,:) )' ), map.windowHu ); colormap gray;

% figure; imshow( imgHuFBPShort(:,:,end/2), map.windowHu ); colormap gray;
% export_fig( 'FBP-DS-short-center.jpg' );
% 
% figure; imshow( imgHuFBPShort(:,:,end/2+40), map.windowHu ); colormap gray;
% export_fig( 'FBP-DS-short-top.jpg' );
% 
% figure; imshow( imgHuFBPShort(:,:,end/2-40), map.windowHu ); colormap gray;
% export_fig( 'FBP-DS-short-bottom.jpg' );
% 
% figure; imshow( fliplr( squeeze( imgHuFBPShort(end/2,:,:) ) ), map.windowHu ); colormap gray;
% export_fig( 'FBP-DS-short-coronal.jpg' );
