load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadMaterialsPhantom(p);

geom = loadProjectionGeometryShortScan( p, 200, 8.7 );

spectrum = loadSpectraCT(p, geom, 1e5);

% % compute sinogram
[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDir );

sinoAtt =  processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);


%% full scan
 map.windowHu  = [-1000 1000];
imgAttFBP = reconFBP( sinoAtt, geom, 'hamming', 1, 1);
imgFBP = convertMonoAttToHu( imgAttFBP, spectrum) ;
figure; imdisp( imgFBP(:,:,end/2)', map.windowHu ); 

return

%%

tubeCurrentProfile = tubeCurrentProfile ./ tubeCurrentProfile(1);

figure;
plot(tubeCurrentProfile)

writetable(table(tubeCurrentProfile),'tubeCurrentProfile_125kV.txt');