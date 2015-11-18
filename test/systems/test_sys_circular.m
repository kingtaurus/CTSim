load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadXCATPhantom(p);

geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e6);
%sinoAtt = simulateAttenuationDataPhantom( phan, geom, spectrum, sinosDirKeV);

% % compute sinogram
[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDirKeV );

sinoRaw = phaserLargeDetectorBinning( sinoRaw );

sinoAtt =  processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);

if 1
    [ sinoAtt, geom ] = truncationCorrectionEllipicalPatch( sinoAtt, geom, 0.16, 4, 1, 64 );
end

%% full scan
 map.windowHu  = [-1000 1000];


imgAttFBP = reconFBP( sinoAtt, geom, 'hamming', 1, 1);
imgFBP = convertMonoAttToHu( imgAttFBP, spectrum) ;
figure; imdisp( imgFBP(:,:,end/2)', map.windowHu ); 
figure; imdisp( fliplr( squeeze(imgFBP(end/2,:,:))), map.windowHu ); 
figure; imdisp( imgFBP(:,:,:), map.windowHu ); 


%% short scan
offsetAngle = 0;

[geomShort, sinoShort ]  = convertSinogram2ShortScan( geom, sinoAtt, offsetAngle );
imgAttShort = reconFBP( sinoShort, geomShort, 'hamming' );
imgFBPShort = convertMonoAttToHu( imgAttShort, spectrum) ;
figure; imdisp( imgFBPShort(:,:,end/2)', map.windowHu ); 
figure; imdisp( fliplr( squeeze(imgFBPShort(end/2,:,:))), map.windowHu ); 
figure; imdisp( imgFBPShort(:,:,:), map.windowHu ); 

return