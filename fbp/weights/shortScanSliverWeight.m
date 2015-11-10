function weight = shortScanSliverWeight( gammas, beta,  Tau)
% function weight = shortScanSliverWeight( gammas, beta,  Tau)
%
% Meng Wu at Stanford University
% 2013


weight = zeros( size(gammas) );
for i = 1 : length( gammas )
    
    gamma = gammas(i);
    if 0 <= beta && beta < Tau + 2 * gamma
        weight(i) = ( sin( pi /4 * beta / ( Tau/2 + gamma)) )^2;
    elseif Tau + 2 * gamma <= beta && beta < pi + 2 * gamma
        weight(i) = 1;
    elseif pi + 2 * gamma <= beta && beta < pi + Tau
        weight(i) = ( sin( pi /4 * ( pi + Tau - beta) / ( Tau/2 - gamma)) )^2;
    else
        weight(i) = 0;
    end
end

end


