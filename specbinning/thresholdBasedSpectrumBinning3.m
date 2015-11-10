function [BinOpt, PhiOpt, ThetaOpt ] = thresholdBasedSpectrumBinning3( spectrum, weights, phis, thetas )
% Threshold based spectrum binning, three energy bins version
%   input:
%       spectrum    - spectrum infomation
%       weights     - weighted norm
%       phis        - list of photo-electric component line intergrals
%       thetas      - list of Compoton scattering component line intergrals
%   output:
%       BinOpt      - optimal bin sizes
%       PhiOpt      - optimal photo-electric component values
%       ThetaOpt    - optimal Compoton scattering values
%
%   Note:
%   1. This is the old version for spectrum binning, in which the number of
%   energy bins are fixed.
%   2. Since this approach is wrose than the generalized spectrum binning,
%   I didn't write a new version for it.
%
%
% Meng Wu at Stanford University
% 2013

normType = 1;

energyLabels            = spectrum.energyBinLabels;
photonsPerEnergyBin     = spectrum.photonsPerEnergyBinOriginal * spectrum.DQE;
energyAverage           = spectrum.energyAverage;
totalPhotons            = sum( photonsPerEnergyBin );
noEnergyBins            = length( photonsPerEnergyBin );

% fit the material to the photon electric and compton scattering
PhiAll = photoElectric(energyLabels)./photoElectric(energyAverage);
ThetaAll = kleinNishina(energyLabels)./kleinNishina(energyAverage);

% ground truth
y = log(  sum ( repmat( photonsPerEnergyBin , [1 length(thetas)] ) ...
    .*  exp( - PhiAll * phis' - ThetaAll * thetas' ) ) )';


scoreFinal = inf;
scoreL1 = zeros(noEnergyBins, noEnergyBins);
scoreL2 = zeros(noEnergyBins, noEnergyBins);
scoreLinf = zeros(noEnergyBins, noEnergyBins);

for i = 1:noEnergyBins
    for j = i:noEnergyBins
        
        b = [sum( photonsPerEnergyBin(1:i)) sum(photonsPerEnergyBin(i:j)) sum(photonsPerEnergyBin(j:end)) ];
        b = b * totalPhotons / sum(b);
        
        e1 =  sum ( photonsPerEnergyBin(1:i) .* energyLabels(1:i) ) / sum(photonsPerEnergyBin(1:i));
        e2 =  sum ( photonsPerEnergyBin(i:j) .* energyLabels(i:j) ) / sum(photonsPerEnergyBin(i:j));
        e3 =  sum ( photonsPerEnergyBin(j:end) .* energyLabels(j:end) ) / sum(photonsPerEnergyBin(j:end));
        
        Phi   = interp1geom(energyLabels, PhiAll, [e1 e2 e3]) ;
        Theta = interp1geom(energyLabels, ThetaAll, [e1 e2 e3]);
        
        x1 = b(1) * exp( - Phi(1) * phis -  Theta(1) * thetas );
        x2 = b(2) * exp( - Phi(2) * phis -  Theta(2) * thetas );
        x3 = b(3) * exp( - Phi(3) * phis -  Theta(3) * thetas );
        
        z = log(x1 + x2 + x3 ) - y;
        
        scoreLinf(i,j) = max( abs(z) .* double(weights > 0 ) );
        scoreL1(i,j) = mean( abs(z) .* weights);
        scoreL1(i,j) = sqrt( mean( z.^2.* weights )  );
        
        scoreLinf(j,i) = scoreLinf(i,j) ;
        scoreL1(j,i) = scoreL1(i,j) ;
        scoreL1(j,i) = scoreL1(i,j) ;
        
        if normType == 0
            score = scoreLinf(i,j);
        elseif normType == 1
            score = scoreL1(i,j);
        else
            score = scoreL2(i,j);
        end
        
        if score < scoreFinal
            BinOpt = b;
            PhiOpt = Phi;
            ThetaOpt = Theta;
            scoreFinal = score;
        end
        
    end
end

if true
    figure;
    imagesc( scoreL1, [0 0.5]); axis xy;
    xlabel( 'Threshold 1 (KeV)','FontSize',16);
    ylabel( 'Threshold 2 (KeV)','FontSize',16);
    colorbar;
    set(gca,'FontSize',16)
end

fprintf('Threshold-based spectrum binning with L%1.0f minimization with %2.5g. (Old version) \n', normType, scoreFinal );
if weights(1) ~= 1
    fprintf('\tNote: use weighted norm \n');
end
fprintf('\tBin \t|\tPhi \t|\tTheta\t \n');
fprintf('------------------------------------\n')
for i = 1 : 3
    fprintf('\t%1.3f\t|\t%1.3f\t|\t%1.3f\t \n',  BinOpt(i) / sum(BinOpt), PhiOpt(i), ThetaOpt(i));
end
fprintf( 'Done.\n\n' );



end


