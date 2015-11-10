function [BinOpt, PhiOpt, ThetaOpt ] = thresholdBasedSpectrumBinning2( spectrum, weights, phis, thetas )
% Threshold based spectrum binning, two energy bins version
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

scoreL1 = zeros(1, noEnergyBins);
scoreL2 = zeros(1, noEnergyBins);
scoreLinf = zeros(1, noEnergyBins);

for i = 1:noEnergyBins
    
    b = [sum( photonsPerEnergyBin(1:i))  sum(photonsPerEnergyBin(i:end)) ];
    b = b * totalPhotons / sum(b);
    
    e1 =  sum ( photonsPerEnergyBin(1:i) .* energyLabels(1:i) ) / sum(photonsPerEnergyBin(1:i));
    e2 =  sum ( photonsPerEnergyBin(i:end) .* energyLabels(i:end) ) / sum(photonsPerEnergyBin(i:end));
    
    Phi   = interp1geom(energyLabels, PhiAll, [e1 e2]) ;
    Theta = interp1geom(energyLabels, ThetaAll, [e1 e2]);
    
    x1 = b(1) * exp( - Phi(1) * phis -  Theta(1) * thetas );
    x2 = b(2) * exp( - Phi(2) * phis -  Theta(2) * thetas );
    
    z = log(x1 + x2) - y;
    
    scoreLinf(i) = max( abs(z) .* double(weights > 0 ) );
    scoreL1(i) = mean( abs(z) .* weights);
    scoreL2(i) = sqrt( mean( z.^2.* weights )  );
    
    if normType == 0
        score = scoreLinf(i);
    elseif normType == 1
        score = scoreL1(i);
    else
        score = scoreL2(i);
    end
    
    if score < scoreFinal
        BinOpt = b;
        PhiOpt = Phi;
        ThetaOpt = Theta;
        scoreFinal = score;
    end
    
end
if true
    figure;
    plot( scoreL1); hold on;
    plot( scoreL2, '--');
    plot( scoreLinf, '-.');
    xlabel( 'Threshold (KeV)', 'fontSize', 14);
    legend('L1', 'L2', 'Linf');
    set(gca,'FontSize',16)
end

fprintf('Threshold-based spectrum binning with L%1.0f minimization with %2.5g. (Old version) \n', normType, scoreFinal );
if weights(1) ~= 1
    fprintf('\tNote: use weighted norm \n');
end
fprintf('\tBin \t|\tPhi \t|\tTheta\t \n');
fprintf('------------------------------------\n')
for i = 1 : 2
    fprintf('\t%1.3f\t|\t%1.3f\t|\t%1.3f\t \n',  BinOpt(i) / sum(BinOpt), PhiOpt(i), ThetaOpt(i));
end
fprintf( 'Done.\n\n' );



end


