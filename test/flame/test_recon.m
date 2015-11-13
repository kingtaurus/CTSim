% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e6);

%% load air scan data
dataPathAir = 'E:\Data\NasaFlame\Nov_10_2015_Study\3ppi_interface_60kV_50mA\air_cold_04\';

sinoAttAir = loadTableTopData( dataPathAir, geom );

sinoAttAir = beamHardeningMaterialCorrection(sinoAttAir, spectrum, 'Quartz', 10 );

%% load normal scan data

dataPath = 'E:\Data\NasaFlame\Nov_10_2015_Study\3ppi_interface_60kV_50mA\kr100_01\';

%dataPath = 'E:\Data\NasaFlame\Nov_5_2015_Study\DiffusionFlameLaminar_1\';

sinoAtt = loadTableTopData( dataPath, geom );

sinoAtt = beamHardeningMaterialCorrection(sinoAtt, spectrum, 'Quartz', 10 );

%%

geomAir = geom;
geomAir.reconOffset(3) = geom.reconOffset(3) ;

imgAir = reconFBP( sinoAttAir, geomAir, 'hamming' );

figure(21); 
if geom.reconSize(3) < 40
    imdisp( imgAir, [0 0.5]   );
else
    imdisp( imgAir(end / 2, :, : ), [0 0.5]   );
end

%%
imgKr = reconFBP( sinoAtt, geom, 'hamming' );

figure(22);
if geom.reconSize(3) < 40
    imdisp( imgKr, [0 0.5] );
else
    imdisp( imgKr(end / 2, :, : ), [0 0.5] );
end

%%
imgKrReg = imgKr;

[optimizer, metric] = imregconfig('monomodal');
for i = 1 : size( imgAir, 3 )
    if mod(i, 50) == 0
        fprintf('(%i/%i)... ', i, size(imgAir, 3 ) );
    end
    slice = imregister(imgKr(:,:,i), imgAir(:, :, i), 'rigid', optimizer, metric);
    imgKrReg( :,:,i) = slice;
end

figure(23); imdisp( imgKrReg(:,end/2,:) - imgAir(:,end/2,:)  , [-0.1 0.1] );