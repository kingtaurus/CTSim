function weights = helicalConeAngleWeights( gammas, beta, Delta, tanAlpha, s, uu, vv, segmentLength  )

a = Delta / 2;
weights = zeros( size( tanAlpha) );

for i = 1:length(gammas)
    
    gamma = gammas(i);    
    if abs( beta + pi - 2 * gamma )  < a
        angleCongj = + pi - 2 * gamma;
    elseif abs( beta - pi - 2 * gamma )  < a
        angleCongj = - pi - 2 * gamma;
    end
    
    zConj = angleCongj * s;
    tanAlphaConj = interp2(  uu, vv, tanAlpha, -uu(:,i),  vv(:,i) - zConj );
    tanAlphaConj( isnan( tanAlphaConj ) ) =  1;
    tanAlphaConj( tanAlphaConj == 0 ) = 1;
    
    if strcmpi(segmentLength, 'short')
        w1 = shortScanSliverWeight( gamma, beta,  Delta - pi );
        w2 = shortScanSliverWeight( - gamma, beta + angleCongj,  Delta - pi );
        
        weights(:,i) =  ( w1 * tanAlphaConj ) ./ (  w1 * tanAlphaConj + w2 * tanAlpha(:,i) ) ;
        
    else
        weights(:,i) =  tanAlphaConj ./ ( tanAlphaConj + tanAlpha(:,i) ) ;
    end
    
end

end