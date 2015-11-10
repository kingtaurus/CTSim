function [BinOpt, PhiOpt, ThetaOpt ] = generalizedSpectrumBinning3( spectrum, weights, phis, thetas )
% Generalized spectrum binning, three energy bins version
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
%   2. New version see function generalizedSpectrumBinning()
%
%
% Meng Wu at Stanford University
% 2013


normType = 1;

energyLabels            = spectrum.energyBinLabels;
photonsPerEnergyBin     = spectrum.photonsPerEnergyBin * spectrum.DQE;
energyAverage           = spectrum.energyAverage;
totalPhotons            = sum( photonsPerEnergyBin );

% fit the material to the photon electric and compton scattering
PhiAll = photoElectric(energyLabels)./photoElectric(energyAverage);
ThetaAll = kleinNishina(energyLabels)./kleinNishina(energyAverage);

% ground truth
y = log(  sum ( repmat( photonsPerEnergyBin , [1 length(thetas)] ) ...
    .*  exp( - PhiAll * phis' - ThetaAll * thetas' ) ) )';


scoreFinal = 1;
for ratio1 = 0.2:0.05:0.8
    for ratio2 = ratio1+0.02:0.05:0.9
        
        b = [ratio1*totalPhotons, (ratio2-ratio1)*totalPhotons, (1-ratio2)*totalPhotons];
        
        i = 1;
        j = 1;
        while sum( photonsPerEnergyBin(1:i) ) < b(1) - 10
            i = i+1;
        end
        
        while sum( photonsPerEnergyBin(1:j) ) < b(2) + b(1) - 10
            j = j+1;
        end

        e1 =  sum ( photonsPerEnergyBin(1:i) .* energyLabels(1:i) ) / sum(photonsPerEnergyBin(1:i));
        e2 =  sum ( photonsPerEnergyBin(i:j) .* energyLabels(i:j) ) / sum(photonsPerEnergyBin(i:j));
        e3 =  sum ( photonsPerEnergyBin(j:end) .* energyLabels(j:end) ) / sum(photonsPerEnergyBin(j:end));
        
        %initial with effective energy and some randomness
        Phi   = interp1geom(energyLabels, PhiAll, [e1 e2 e3]) ;
        Theta = interp1geom(energyLabels, ThetaAll, [e1 e2 e3]);
        
        [scoreOld, dPhi, dTheta] = computeGradient( Phi, Theta, phis, thetas, b, y, weights, normType );
        %gradient descent
        alpha = 1;
        while alpha > 1e-4
            
            [scoreNew, dPhiNew, dThetaNew] = computeGradient( Phi + alpha * dPhi, Theta + alpha * dTheta, phis, thetas, b, y, weights, normType );     
            if   scoreOld  < scoreNew + 1e-6
                alpha = alpha / 2;
            else
                Phi =  Phi + alpha * dPhi;
                Theta = Theta + alpha * dTheta;
                dPhi = dPhiNew;
                dTheta = dThetaNew;
                scoreOld = scoreNew;
            end
        end
        
        if scoreOld < scoreFinal
            BinOpt = b;
            PhiOpt = Phi;
            ThetaOpt = Theta;
            scoreFinal = scoreOld;
        end
  
        
    end
end

fprintf('Generalized spectrum binning with L%1.0f minimization with %2.5g. (Old version) \n', normType, scoreFinal );
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


function [score, dPhi, dTheta] = computeGradient( Phi, Theta, phis, thetas, b, y, weights, normType )

x1 = b(1) * exp( - Phi(1) * phis -  Theta(1) * thetas );
x2 = b(2) * exp( - Phi(2) * phis -  Theta(2) * thetas );
x3 = b(3) * exp( - Phi(3) * phis -  Theta(3) * thetas );

z = log(x1 + x2 + x3) - y;
if normType == 2
    score = sqrt( mean( z.^2.*weights));
else
    score = mean( abs(z).*weights);
    z = sign( z );
end

z = z.*weights;

dPhi =  [  mean( z .* x1 .* phis ./ (x1+x2+x3)) ,  mean( z .*  x2 .* phis ./ (x1+x2+x3)),  mean( z .*  x3 .* phis ./ (x1+x2+x3)) ];
dTheta =  [ mean( z .*  x1 .* thetas ./ (x1+x2+x3)) ,  mean( z .* x2 .* thetas ./ (x1+x2+x3)) ,  mean( z .* x3 .* thetas ./ (x1+x2+x3)) ];

end

