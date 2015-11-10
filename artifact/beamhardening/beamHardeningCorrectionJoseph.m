function sinoOut = beamHardeningCorrectionJoseph( imgInitial, geom, spectrum, material1, material2 )
% Compute Beam Hardening Correction imgInitial sinogram space using Joseph's method.
%   imgInitial:
%       imgInitial       - sinogram or image imgInitial attenuation
%       geom     - system geometry
%       spectrum - spectrum infomation
%   output:
%       sinoOut  - correction to add to original sinogram

if nargin < 4
    material1 = 'Tissue_Soft_ICRU-44';
    material2 = 'Bone_Cortical_ICRU-44';
end


if spectrum.energyIntegrating
    spectrum.photonsPerEnergyBin = spectrum.DQE * spectrum.photonsPerEnergyBin .* spectrum.energyBinLabels * spectrum.detectorGain ;
else
    spectrum.photonsPerEnergyBin = spectrum.DQE * spectrum.photonsPerEnergyBin * spectrum.detectorGain  ;
end

photonsTotal = spectrum.photonsTotal * spectrum.DQE;
energyBinLabels = spectrum.energyBinLabels;
photonsPerEnergyBin = spectrum.photonsPerEnergyBin * spectrum.DQE;

[muPolySoft, muMonoSoft] = materialAttenuation( energyBinLabels, material1, photonsPerEnergyBin);
[muPolyBone, muMonoBone] = materialAttenuation( energyBinLabels, material2, photonsPerEnergyBin);

%segment the material
imgInitial = extendVoi( imgInitial, 2 );
[phoSoft, phoBone] = materialSegmentationSharp( imgInitial, muMonoSoft, muMonoBone );

softProjLength = forwardProjectMex( phoSoft, geom );
boneProjLength = forwardProjectMex( phoBone, geom );

% create the look up table
[ sinoDiffTable, softLengthTable, boneLengthTable ] = attenuationLookupTable( muPolySoft, muMonoSoft, ...
    muPolyBone, muMonoBone, photonsPerEnergyBin, photonsTotal  );

softProjLength(softProjLength > 49.999 ) = 99.999;
boneProjLength(boneProjLength > 19.999 ) = 29.999;
sinoOut = interp2( softLengthTable, boneLengthTable, sinoDiffTable, softProjLength, boneProjLength);
sinoOut( isnan( sinoOut )) = 0;

end


function [ sinoDiffTable, softLengthTable, boneLengthTable ] = attenuationLookupTable( muPolySoft, muMonoSoft, ...
    muPolyBone, muMonoBone, photonsPerEnergyBin, photonsTotal  )

softLengths = linspace(-0.1,100,256);
boneLengths = linspace(-0.1,30,256);

[softLengthTable,  boneLengthTable] = meshgrid(softLengths, boneLengths );

sinoDiffTable = zeros( softLengthTable, boneLengthTable );

for i = 1:size(softLengthTable, 1)
    for j = 1:size(softLengthTable, 2)
        
        sinoDiffTable(i,j) = - log(  photonsTotal* exp( - softLengthTable(i,j) * muMonoSoft - boneLengthTable(i,j) * muMonoBone )) ...
            + log( sum( photonsPerEnergyBin.*exp(- softLengthTable(i,j) * muPolySoft - boneLengthTable(i,j) * muPolyBone )));
        
    end
end

end



