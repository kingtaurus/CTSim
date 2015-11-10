function sinoAtt = simulateAttenuationDataCT( phan, geom, spectrum, sinosDir,  addCountingNoise )
% function sinoAtt = simulateAttenuationDataCT( phan, geom, averageEnergy, photonsTotal, backgroundEvents )
% Compute projection image from CT scan data
% only simulate monochromatic data
% Meng Wu at Stanford University
% 2013

photonsTotal            = spectrum.photonsTotal;
averageEnergy           = spectrum.energyAverage;


fprintf('Compute expected attenuation sinogram from CT image ... \n');

sinoAttFilename = sprintf([sinosDir 'sinoAttData_CT_kVp%3.0f_FLUX%1.2g_NOISE%i.mhd'], ...
    max(spectrum.energyBinLabels), photonsTotal, addCountingNoise );


% get attunation phantom image
muWater = materialAttenuation( averageEnergy, 'Water_Liquid', 1);
if min( phan.phantomMaterials(:) ) <  -500
    imgAtt = ( phan.phantomMaterials + 1000 ) * muWater / 1000;
else
    imgAtt =  phan.phantomMaterials;
end

if length( geom.reconSize ) == 2 && length( geom.detSize ) == 1
    phan.reconSpacing = phan.reconSpacing(1:2);
    phan.reconSize = phan.size(1:2);
    phan.reconOffset = phan.offset(1:2);
    imgAtt = imgAtt(:,:,round(end/2));
end

if spectrum.useBowtie
    fprintf('Warning: bowtie filter is not supported.\n');
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