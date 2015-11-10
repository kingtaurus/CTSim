function [BinOpt, PhiOpt, ThetaOpt] = generalizedSpectrumBinning( spectrum, weights, numOfBins, phis, thetas )
% Generalized spectrum binning with arbitrary number of energy bins
%   input:
%       spectrum    - spectrum infomation
%       weights     - weighted norm
%       numOfBins   - number of energy bins
%       phis        - list of photo-electric component line intergrals
%       thetas      - list of Compoton scattering component line intergrals
%   output:
%       BinOpt      - optimal bin sizes
%       PhiOpt      - optimal photo-electric component values
%       ThetaOpt    - optimal Compoton scattering values
% Note:
%   1. New version with recursion
%   2. Very slow when number of energy bins is larger than three, so be
%   patient.
%
% Meng Wu at Stanford University
% 2014

energyBinLabels         = spectrum.energyBinLabels;
photonsPerEnergyBin     = spectrum.photonsPerEnergyBinOriginal * spectrum.DQE;
energyAverage           = spectrum.energyAverage;

% fit the material to the photon electric and compton scattering
PhiAll = photoElectric(energyBinLabels)./photoElectric(energyAverage);
ThetaAll = kleinNishina(energyBinLabels)./kleinNishina(energyAverage);


% truth attenuation line intergral
trueAtt =  log(  sum ( repmat( photonsPerEnergyBin , [1 length(thetas)] ) ...
    .*  exp( - PhiAll * phis' - ThetaAll * thetas' ) ) )';

params.trueAtt = trueAtt;
params.photonsPerEnergyBin = photonsPerEnergyBin;
params.energyBinLabels = energyBinLabels;
params.PhiAll = PhiAll;
params.ThetaAll = ThetaAll;
params.phis = phis;
params.thetas = thetas;
params.weights = weights;
params.numOfBins = numOfBins;

% recursively loop all binning combinations
[score, BinOpt, PhiOpt, ThetaOpt] = spectrumBinningSizes( params, 1, -100 );

% printout results
fprintf('Generalized spectrum binning with L1 minimization with %2.5g. (New version) \n', -score );
if weights(1) ~= 1
    fprintf('\tNote: use weighted norm \n');
end
fprintf('\tBin \t|\tPhi \t|\tTheta\t \n');
fprintf('------------------------------------\n')
for i = 1 : numOfBins
    fprintf('\t%1.3f\t|\t%1.3f\t|\t%1.3f\t \n',  BinOpt(i) / sum(BinOpt), PhiOpt(i), ThetaOpt(i));
end
fprintf( 'Done.\n\n' );

end

function [score, BinSizes, Phi, Theta ] = spectrumBinningSizes(  params, level, scoreOld, thresholds,  BinSizesOld, PhiOld, ThetaOld )
% doing spectrum binning for arbitrary number of bins, recusively


numOfBins = params.numOfBins;
if nargin < 4
    thresholds = zeros(1, numOfBins-1 );
    BinSizesOld = zeros(1, numOfBins );
    PhiOld = zeros(1, numOfBins );
    ThetaOld = zeros(1, numOfBins );
    
end

if level == numOfBins % tail recursion
    
    propotions = [thresholds 1] - [0 thresholds];
    
    % determine the engery bin sizes
    BinSizes = propotions * sum( params.photonsPerEnergyBin );
    
    % initial corresponding energy for each bins
    binAverageEnergies = zeros( 1, numOfBins );
    i = 1;    j = 1;
    for b = 1 : numOfBins
        while sum( params.photonsPerEnergyBin(1:j)  ) < sum( BinSizes(1:b) ) + 10 && j < length( params.photonsPerEnergyBin )
            j = j + 1;
        end
        binAverageEnergies(b) = sum ( params.photonsPerEnergyBin(i:j) .* params.energyBinLabels(i:j) ) / sum( params.photonsPerEnergyBin(i:j));
        i = j;
    end
    
    % initialized component parameters
    Phi   = interp1geom(params.energyBinLabels, params.PhiAll, binAverageEnergies);
    Theta = interp1geom(params.energyBinLabels, params.ThetaAll, binAverageEnergies);
    
    %gradient descent with backline tracking
    [errorOld, dPhi, dTheta] = spectrumBinningGradient( Phi, Theta, BinSizes, ...
        params.phis, params.thetas, params.trueAtt, params.weights );
    
    %gradient descent
    alpha = 1;
    while alpha > 1e-4
        
        [errorNew, dPhiNew, dThetaNew] = spectrumBinningGradient( Phi + alpha * dPhi, Theta + alpha * dTheta, BinSizes, ...
            params.phis, params.thetas, params.trueAtt, params.weights );
        if   errorOld  < errorNew + 1e-6
            alpha = alpha / 2; % reduce step size
        else
            Phi =  Phi + alpha * dPhi; % gradient descent
            Theta = Theta + alpha * dTheta;
            dPhi = dPhiNew;
            dTheta = dThetaNew;
            errorOld = errorNew;
        end
    end
    
    % record the local optimal
    score = - errorOld ;
    
else
    
    incrementSize = 0.01 * numOfBins;
    if level == 1
        thresholdStart = 0.1;
    else
        thresholdStart = thresholds( level - 1 ) + 0.1;
    end
    
    thresholdStop = 1 - 0.1 * (numOfBins - level) ;
    
    % recursion loop to set thresholds
    for t = thresholdStart : incrementSize : thresholdStop
        
        thresholds( level ) = t;
        [score, BinSizes, Phi, Theta] = spectrumBinningSizes( params, level + 1,  scoreOld, thresholds, BinSizesOld, PhiOld, ThetaOld  );
        
        % save better binning values
        if score > scoreOld
            scoreOld = score;
            BinSizesOld = BinSizes;
            PhiOld = Phi;
            ThetaOld = Theta;
        end
        
    end
    
    score = scoreOld;
    BinSizes = BinSizesOld;
    Phi     = PhiOld;
    Theta   = ThetaOld;
end

end


function [error, dPhi, dTheta] = spectrumBinningGradient( Phi, Theta, BinSizes, phis, thetas, trueAtt, weights )
% gradient of weighted norm against phis and thetas

normType = 1;
numOfBins = length( BinSizes );
numOfCombination = length( phis );

x = zeros(numOfCombination, numOfBins);
for ib = 1: numOfBins
    x( :, ib ) =  BinSizes( ib ) * exp( - Phi( ib ) * phis -  Theta( ib ) * thetas );
end

% attenuation line itergral
y = sum( x, 2);
z = log( y ) - trueAtt;

% computed weighted norm
if normType == 2
    error = sqrt( mean( z.^2.*weights));
else
    error = mean( abs(z).*weights);
    z = sign( z );
end
z = z.*weights;

% compute gradient
dPhi = zeros(1, numOfBins);
dTheta = zeros(1, numOfBins);
for ib = 1: numOfBins
    dPhi( ib ) = mean( z .* x(:,ib) .* phis ./ y ) ;
    dTheta( ib ) = mean( z .* x(:,ib) .* thetas ./ y ) ;
end

end

