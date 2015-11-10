function projThicknessCounts = computeEffectiveThickness( materialMap1, materialMap2, rho, geom, spectrum, ...
    material1, material2, phis, thetas, thickness1, thickness2  )

downsampleRate = 4;
[ phi1, theta1 ] = attenuationCoefficientDecompose( material1, spectrum );
[ phi2, theta2 ] = attenuationCoefficientDecompose( material2, spectrum );

% downsample in recon
geom.detSize = round( geom.detSize / downsampleRate );
geom.detSpacing = geom.detSpacing * downsampleRate;
geom.betas = geom.betas(1:downsampleRate:end);
geom.couchZ = geom.couchZ(1:downsampleRate:end);
geom.noViews = length( geom.betas );

proj1 = forwardProjectMex( materialMap1.*rho, geom );
proj2 = forwardProjectMex( materialMap2.*rho, geom );

%figure; imdisp( proj1 );
%figure; imdisp( proj2 );


phisProj = proj1 * phi1 + proj2 * phi2;
thetaProj = proj1 * theta1 + proj2 * theta2;



%% Add bowtie into consideration
if spectrum.useBowtie
    
    bowtieThickness = spectrum.bowtieThickness(1:downsampleRate:end, 1:downsampleRate:end) / 10; % mm to cm
    [ phiFiltration, thetaFiltration ] = attenuationCoefficientDecompose( spectrum.bowtieMaterial, spectrum );
    
    if ndims( phisProj ) == 2
        for iv = 1 : size(phisProj, 2)
            phisProj(:, iv) = phisProj(:, iv) + bowtieThickness .* phiFiltration;
            thetaProj(:, iv) = thetaProj(:, iv) + bowtieThickness .* thetaFiltration;
        end
    else
        for iv = 1 : size(phisProj, 3)
            phisProj(:, :, iv) = phisProj(:, :, iv) + bowtieThickness .* phiFiltration;
            thetaProj(:, :, iv) = thetaProj(:, :, iv) + bowtieThickness .* thetaFiltration;
        end
        
    end
end

phisProj = phisProj(:);
thetaProj = thetaProj(:);

%% count number of pixels with given thickness
sampleLength = length( phis );
projThicknessCounts = zeros( sampleLength, 1 );

dphis = abs(phis(2) - phis(1)) /2 ;
dtheta = abs(thetas( sqrt(sampleLength) +1) - thetas(1)) / 2 ;

for i = 1 : sampleLength
    
    projThicknessCounts(i) =  sum( ( phisProj >= phis(i) - dphis ) & ( phisProj <= phis(i) + dphis ) ...
        & ( thetaProj >= thetas(i) - dtheta ) &  ( thetaProj <= thetas(i) + dtheta ) );
    
end

%%

if nargin < 10
    figure;
    imagesc(phis, thetas,  reshape(projThicknessCounts, [sqrt(sampleLength) sqrt(sampleLength)]) );
    colorbar; axis xy;
    xlabel 'Photon-Eletric Effects '; ylabel 'Compton Scattering';
    
else
    figure;
    imagesc(thickness1, thickness2,  reshape(projThicknessCounts, [sqrt(sampleLength) sqrt(sampleLength)]) );
    colorbar; axis xy;
    xlabel 'Material 1 (cm)'; ylabel 'Material 2 (cm)';
end