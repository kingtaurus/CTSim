function sinoOut = binningDetecotRows( sinoIn, k )
% combine multiples rows of the detector for FDK method
%
% Meng Wu at Stanford Univeristy
% 2012 - 2013

if nargin < 2
    k = 2;
end

if k <= 1
    sinoOut = sinoIn;
    return;
end

if mod( k , 2) ==  0
    h = [ 1/2; ones(k - 1, 1); 1/2 ] / k ;
else
    h = ones(k, 1) / k;
end

if length( size( sinoIn )) == 3
    
    sinoOut = sinoIn;
    
    for iv = 1: size( sinoIn, 3)
       
        % entire sinogram
        sinoOut(:,:,iv) = filter2( h, sinoIn(:,:,iv), 'same');
        
    end
    
end    
    
    
end