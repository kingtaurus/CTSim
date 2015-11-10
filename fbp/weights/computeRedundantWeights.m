function redundantWeights = computeRedundantWeights( gammas, betas, iview, nv, shortScanWeights )
% function redundantWeights = computeRedundantWeights( gammas, betas, iview, nv, segmentLength )
%
%
% Meng Wu, 2014.3

noViews     = length( betas );
Delta       = abs(betas(end) - betas(1)) / (noViews - 1) * noViews;
betaMin     = min(betas);

if shortScanWeights
    weights = shortScanSliverWeight( gammas, betas(iview) - betaMin,  Delta - pi );
else
    weights = helicalSimpleRedundancyWeights( gammas, betas(iview) - betaMin, Delta );
end
% weights = helicalSimpleRedundancyWeights( gammas, betas(iview) - betaMin, Delta );

redundantWeights = repmat( weights, nv, 1);

end