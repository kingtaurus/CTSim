load 'temp.mat';
close all;


%% Load simulation parameters and datas

%[phan, map, roi] = loadMaterialsPhantom(p);

geom = loadProjectionGeometryCT( p );
spectrum = loadSpectraCT(p, geom);

energyBin               = spectrum.energyBinLabels;
photonsPerEnergyBin     = spectrum.photonsPerEnergyBin * spectrum.DQE;
noEnergyBins            = length(energyBin);


material1 = 'Tissue_Soft_ICRU-44';
material2 = 'Bone_Cortical_ICRU-44';
material3 = 'Aluminum';

muSoftPoly = materialAttenuation( energyBin, material1, photonsPerEnergyBin );
muBonePoly = materialAttenuation( energyBin, material2, photonsPerEnergyBin );
muFiltPoly  = materialAttenuation( energyBin, material3, photonsPerEnergyBin );


figure;
semilogy(energyBin, muBonePoly ); hold on;
semilogy(energyBin, muSoftPoly, '-.' );
xlabel( 'Energy (KeV)', 'fontSize', 16);
ylabel( 'Linear Attenuation (cm^{-1})', 'fontSize', 16);
legend( 'Bone', 'Soft tissue', 'fontSize', 16)

figure;
plot(energyBin, photonsPerEnergyBin );
xlabel( 'Energy (KeV)', 'fontSize', 16);
ylabel( 'Fluence (Photons/pixel)', 'fontSize', 16);




%% compute effective energy and photon counts

muSoftLow = zeros( noEnergyBins );
muBoneLow = zeros( noEnergyBins );

muSoftMid = zeros( noEnergyBins );
muBoneMid = zeros( noEnergyBins );

muSoftHigh = zeros( noEnergyBins );
muBoneHigh = zeros( noEnergyBins );

for i = 1 : noEnergyBins
    for j = i : noEnergyBins
        
        
        energyLow =  sum ( photonsPerEnergyBin(1:i) .* energyBin(1:i) ) / sum(photonsPerEnergyBin(1:i));
        energyMid =  sum ( photonsPerEnergyBin(i:j) .* energyBin(i:j) ) / sum(photonsPerEnergyBin(i:j));
        energyHigh =  sum ( photonsPerEnergyBin(j:end) .* energyBin(j:end) ) / sum(photonsPerEnergyBin(j:end));
        
        muSoftLow(i, j) = materialAttenuation( energyLow, material1, 1 );
        muBoneLow(i, j) = materialAttenuation( energyLow, material2, 1 );
        
        muSoftMid(i, j) = materialAttenuation( energyMid, material1, 1 );
        muBoneMid(i, j) = materialAttenuation( energyMid, material2, 1 );
        
        muSoftHigh(i, j) = materialAttenuation( energyHigh, material1, 1 );
        muBoneHigh(i, j) = materialAttenuation( energyHigh, material2, 1 );
        
        muSoftLow(j, i) = muSoftLow(i, j);
        muBoneLow(j, i) = muBoneLow(i, j);
        
        muSoftMid(j, i) = muSoftMid(i, j);
        muBoneMid(j, i) = muBoneMid(i, j);
        
        muSoftHigh(j, i) = muSoftHigh(i, j);
        muBoneHigh(j, i) = muBoneHigh(i, j);
        
    end
    
end

%% compute optimal threshold
maxThickness1 = 30;
maxThickness2 = 10;
maxThickness3 = 0;

softThickness = linspace(0, maxThickness1, 8 );
boneThickness = linspace(0, maxThickness2, 5 );
totalPhotonLow = zeros( noEnergyBins );
totalPhotonMid = zeros( noEnergyBins );
totalPhotonHigh = zeros( noEnergyBins );

for i = 1 : noEnergyBins
    for j = i : noEnergyBins
        totalPhotonLow(i, j) = sum( photonsPerEnergyBin(1:i).*exp( - muFiltPoly(1:i) * maxThickness3) );
        totalPhotonMid(i, j) = sum( photonsPerEnergyBin(i:j).*exp( - muFiltPoly(i:j) * maxThickness3) );
        totalPhotonHigh(i, j) = sum( photonsPerEnergyBin(j:end).*exp( - muFiltPoly(j:end) * maxThickness3) );
        
        totalPhotonLow(j, i) = totalPhotonLow(i, j);
        totalPhotonMid(j, i) = totalPhotonMid(i, j);
        totalPhotonHigh(j, i) = totalPhotonHigh(i, j);
    end
    
end

%%

maxDiff = zeros( size(totalPhotonLow) );


for i = 1:length(softThickness)
    for j = 1: length(boneThickness)
        
        ts = softThickness(i);
        tb = boneThickness(j);
        
        sinoPhotonsTrue =  sum ( photonsPerEnergyBin .* ...
            exp( - muSoftPoly * ts - muBonePoly * tb - muFiltPoly * maxThickness3 ));
        
        sinoPhotonApprox = totalPhotonLow .* exp( - muSoftLow * ts - muBoneLow * tb ) ...
            + totalPhotonMid  .* exp( - muSoftMid * ts - muBoneMid * tb ) ...
            + totalPhotonHigh .* exp( - muSoftHigh * ts - muBoneHigh * tb );
        
        attenDiff =  abs( log( sinoPhotonsTrue ) - log( sinoPhotonApprox ) );
        
        maxDiff = max( maxDiff,  attenDiff );
        
    end
end

figure;
imagesc( maxDiff, [0 1] );

xlabel( 'Threshold Energy \tau (KeV)', 'fontSize', 16);
ylabel( 'Threshold Energy \tau (KeV)', 'fontSize', 16);
title( 'Difference', 'fontSize', 16);

%%

energies = meshgrid(energyBin );

optIndex = find( maxDiff == min( maxDiff(:)));

optDiff =  maxDiff( optIndex(1) );
optCutoff1 = energies( optIndex(1) );
optCutoff2 = energies( optIndex(2) );

fprintf( 'The optimal cutoff energy is %3.1f and %3.1f keV, with min max difference %2.2f. \n', optCutoff1, optCutoff2, optDiff);

photonCounts = [totalPhotonLow(optIndex) totalPhotonMid(optIndex) totalPhotonHigh(optIndex)];
muMaterail1  = [muSoftLow(optIndex) muSoftMid(optIndex) muSoftHigh(optIndex)];
muMaterail2  = [muBoneLow(optIndex) muBoneMid(optIndex) muBoneHigh(optIndex)];


