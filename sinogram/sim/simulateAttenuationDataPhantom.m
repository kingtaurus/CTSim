function sinoAtt = simulateAttenuationDataPhantom( phan, geom, spectrum, sinoDir, addCountingNoise )
% function sinoAtt = simulateAttenuationDataPhantom( phan, geom, spectrum, addCountingNoise )
% only simulate monochromatic data
% fixed tube mAs 
% input:
%       phan        - phantom parameters
%       geom        - geometry parameters
%       spectrum    - spectrum parameters
%       addCountingNoise    - default false
% output:
%       sinoAtt
%
% Meng Wu at Stanford University
% 2012 - 2013

if nargin < 5
    addCountingNoise = false;
end

spectrumEnergies        = spectrum.energyBinLabels;
spectrumPhotonsPerBin   = spectrum.photonsPerEnergyBin;
photonsTotal            = spectrum.photonsTotal;

if ~exist(sinoDir, 'dir'), mkdir(sinoDir); end
fprintf('Compute expected attenuation sinogram from numerical phantom ... \n');

sinoAttFilename = sprintf([sinosDir 'sinoAttData_%s_kVp%3.0f_FLUX%1.2g_NOISE%i.mhd'], ...
    phan.materialMappingName, max(spectrum.energyBinLabels), photonsTotal, addCountingNoise );

if exist(sinoAttFilename, 'file')
    fprintf(' use saved sinogram data. \n');
    sinoAtt = readMetaImage(sinoAttFilename);
    return;
end

if spectrum.useBowtie
    fprintf('Warning: bowtie filter is not supported.\n');
end

% get attunation phantom image
noMaterials = length( phan.materialIndexes );
imgAtt = zeros( size(phan.phantomDensities), 'single');
for m = 1:noMaterials
    
    [ ~, materialAtt] = materialAttenuation( spectrumEnergies, phan.materialNames{m}, spectrumPhotonsPerBin);
    % convert material phantom to attenuation image
    imgMaterial = single( phan.phantomMaterials == phan.materialIndexes{m} );
    
    if phan.useDensity
        imgAtt =  imgAtt + imgMaterial .* phan.phantomDensities * materialAtt;
    else
        imgAtt =  imgAtt + imgMaterial * materialAtt;
    end
end

if length( geom.reconSize ) == 2 && length( geom.detSize ) == 1
    phan.reconSpacing = phan.reconSpacing(1:2);
    phan.reconSize = phan.size(1:2);
    phan.reconOffset = phan.offset(1:2);
    imgAtt = imgAtt(:,:,round(end/2));
end

% forward projection
sinoAtt =  forwardProjectPhantomMex( imgAtt, phan, geom, 'proj,tf' );

% add Poisson noise if necessary
if addCountingNoise
    fprintf('\tsimulating counting noise... ');
    
    sinoPhotonCounts = poissrnd( photonsTotal * exp( - sinoAtt ) );
    sinoAtt = log(photonsTotal) - log( sinoPhotonCounts );
    
end

writeMetaImage(sinoAtt, sinoAttFilename);
fprintf('Done\n\n');

end