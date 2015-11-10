function printOutSpectrumInfo( spectrum )

fprintf('Spectra summary:\n');
fprintf('\tNumber of photons per pixel: %4.3g \n', spectrum.photonsTotal );
fprintf('\tNumber of photons per 1mm sqr at 1m pixel: %4.3g \n', spectrum.photonsPerMm2At1m );

if spectrum.useBowtie
    fprintf('\tUse bowtie filter with thickness from %2.2f to %2.2f mm. \n', min(spectrum.bowtieThickness(:)), max(spectrum.bowtieThickness(:)) );
    fprintf('\tBowtie flat field ratio from %1.3f to %1.3f. \n', max(spectrum.flatFieldRatio(:)), min(spectrum.flatFieldRatio(:)) );
end

if  isfield(spectrum,'automaticExposureControl')
    if spectrum.automaticExposureControl
        fprintf('\tUse automatic exposure control.\n');
    end
end

if isfield(spectrum,'energyIntegrating')
    if spectrum.energyIntegrating
        fprintf('\tUse energy integrating detector.\n');
    end
end

if isfield(spectrum,'compoundPoissonNoise')
    if spectrum.compoundPoissonNoise
        fprintf('\tUse compound Possion noise detector.\n');
    end
end

% print spectra summaries
photonsEffectiveThroughTissue    = computeResultingPhotons( spectrum.photonsPerEnergyBin ...
    * spectrum.DQE, spectrum.energyBinLabels, 'Tissue_Soft_ICRU-44', 200);

fprintf('\t Beam | # ph. (P) | # ph. (Eff) | # ph. (Tissue) | E_min (keV) | E_max (keV) | E_avg (keV) \n');
fprintf('\t------+-----------+-------------+----------------+-------------+-------------+-------------\n');
fprintf('\t keV  | %9.0g | %9.0g   | %9.0g      |   %7.2f   |   %7.2f   |   %7.2f   \n', ...
    spectrum.photonsTotal, spectrum.photonsTotal * spectrum.DQE, ...
    sum(photonsEffectiveThroughTissue), spectrum.energyBinLabels(1), spectrum.energyBinLabels(end), spectrum.energyAverage);


end