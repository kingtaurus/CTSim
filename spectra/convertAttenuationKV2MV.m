function imgAttMeV = convertAttenuationKV2MV( imgAttKeV, spectrumKeV, spectrumMeV, maxAtten )
% convert attenuation in KeV image to the corresponding values in MeV image
% input:
%       imgAttKeV   - attenuation in KeV image
%       spectrumKeV - X-ray spectrum infomation of KeV
%       spectrumMeV - X-ray spectrum infomation of MeV
% output:
%       imgAttMeV   - attenuation in MeV image
%
% Meng Wu at Stanford
% 2012.10

if spectrumKeV.energyBinLabels(end) == 119
    %muKeV = [0.00, 0.200, 0.225, 0.375, 0.480, 0.600, 1.200, 1.800];
    muKeV = [0.00, 0.200, 0.230, 0.407, 0.595, 1.000, 2.00];
else
    muKeV = [0 0.207 0.216 0.43 0.86];
end

if spectrumMeV.energyBinLabels(end) == 2484
    %muMeV = [0.00, 0.084, 0.090, 0.120, 0.140, 0.150, 0.20, 0.25];
    muMeV = [0.00, 0.081, 0.087, 0.114, 0.145, 0.20, 0.35];
else
    %muMeV = [0.00, 0.060, 0.066, 0.085, 0.105, 0.120, 0.170, 0.210];
    muMeV = [0.00, 0.053, 0.055, 0.071, 0.090, 0.150, 0.25];
end

muKeV = [0 0.206 0.215 0.433 0.866 1.3];
muMeV = [0 0.05 0.052 0.0697 0.125 0.19];


maxAtten = min( maxAtten,(muKeV(end) - 0.001));

% make sure no pixel is outside the lookup table
imgAttKeV( imgAttKeV <= 0 ) = 0;
imgAttKeV( imgAttKeV >= maxAtten ) =  maxAtten;


% attenuation conversion using lookup table
imgAttMeV = interp1( muKeV, muMeV, imgAttKeV);
