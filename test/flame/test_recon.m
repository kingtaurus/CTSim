% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e6);

%% load air scan data
dataPathAir = 'E:\Data\NasaFlame\Nov_6_2015_Study\combusion_10ppisic_60kV_50ma\air_scan\';

sinoAttAir = loadTableTopData( dataPathAir, geom );

%% load normal scan data

dataPath = 'E:\Data\NasaFlame\Nov_6_2015_Study\combusion_10ppisic_60kV_50ma\burn_02\';

%dataPath = 'E:\Data\NasaFlame\Nov_5_2015_Study\DiffusionFlameLaminar_1\';

sinoAtt = loadTableTopData( dataPath, geom );

%%

sinoAttAirBHC = beamHardeningMaterialCorrection(sinoAttAir, spectrum, 'Quartz', 10 );
geomAir = geom;
geomAir.reconOffset(3) = geom.reconOffset(3) + 0.1;

imgAir = reconFBP( sinoAttAirBHC, geomAir, 'hamming' );

figure(21); imdisp( imgAir, [0 0.5]   );

%%

sinoAttBHC = beamHardeningMaterialCorrection(sinoAtt, spectrum, 'Quartz', 10 );

imgKr = reconFBP( sinoAttBHC, geom, 'hamming' );

figure(22); imdisp( imgKr, [0 0.5] );

%%
imgAirReg = imgAir;

[optimizer, metric] = imregconfig('monomodal');
for i = 1 : size( imgAir, 3 )
    if mod(i, 50) == 0
       fprintf('(%i/%i)... ', i, size(imgAir, 3 ) ); 
    end
    slice = imregister(imgAir(:, :, i), imgKr(:,:,i), 'rigid', optimizer, metric);
    imgAirReg( :,:,i) = slice;
end

figure; imdisp( imgKr(:,end/2,:) - imgAirReg(:,end/2,:)  , [0 0.1] );