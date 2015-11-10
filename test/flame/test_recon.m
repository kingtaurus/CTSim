% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e6);

%% load air scan data
dataPathAir = 'E:\Data\NasaFlame\Nov_6_2015_Study\combusion_10ppisic_60kV_50ma\air_scan\';

sinoAttAir = loadTableTopData( dataPathAir, geom );

%% load normal scan data

dataPath = 'E:\Data\NasaFlame\Nov_6_2015_Study\combusion_10ppisic_60kV_50ma\burn_01\';

dataPath = 'E:\Data\NasaFlame\Nov_5_2015_Study\DiffusionFlameLaminar_1\';

sinoAtt = loadTableTopData( dataPath, geom );

%%

sinoAttAirBHC = beamHardeningMaterialCorrection(sinoAttAir, spectrum, 'Quartz', 10 );
geomAir = geom;
geomAir.reconOffset(3) = geom.reconOffset(3) + 0.5;

imgAir = reconFBP( sinoAttAirBHC, geomAir, 'hamming' );

figure(21); imdisp( imgAir  );

%%

sinoAttBHC = beamHardeningMaterialCorrection(sinoAtt, spectrum, 'Quartz', 10 );

imgKr = reconFBP( sinoAtt, geom, 'hamming' );

figure(22); imdisp( imgKr(:,end/2,:), [0 0.05] );

%%
imgAirReg = imgAir;

[optimizer, metric] = imregconfig('monomodal');
for i = 1 : size( imgAir, 3 )
    slice = imregister(imgAir(:, :, i), imgKr(:,:,i), 'rigid', optimizer, metric);
    imgAirReg( :,:,i) = slice;
end

figure; imdisp( abs( imgKr - imgAirReg ) , [0 0.1] );