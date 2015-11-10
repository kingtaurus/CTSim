function sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile, waterCorrection )
% function sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile, waterCorrection )
%   Log and normalize the raw detector data with automatic exposure control
%   With an additional warter correction
%
% Meng Wu at Stanford University
% 2014

if nargin < 4
    waterCorrection = true;
end

fprintf('Process the Raw data: ... \n');

intensity = sinoRaw;
if ndims(sinoRaw) == 3
    
    for iv = 1 : size( sinoRaw, 3 )
        intensity(:, :, iv) = sinoRaw(:, :, iv) / tubeCurrentProfile( iv );
    end
    
else
    for iv = 1 : size( sinoRaw, 2 )
        intensity(:, iv) = sinoRaw(:, iv) / tubeCurrentProfile( iv );
    end
    
end

fprintf('\tconverting to attenuation integrals by -log(I/I_0)... \n');

if spectrum.compoundPoissonNoise
    sinoAtt = computeCompoundPoissonSinogramAttunation( intensity, spectrum );
else
    sinoAtt = computeSimplePoissonSinogramAttunation( intensity, spectrum );
end


if waterCorrection
    fprintf('\twater base beam hardening correction ... \n');
    sinoAtt = beamHardeningWarterCorrection(sinoAtt, spectrum, true );
end

fprintf('Done\n\n');