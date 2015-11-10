function I = immedian3( I, n )
% Filter the sinogram with 2d median filter
% Input:
%       I -  sinogram to process
%       n    -  half size of median filter
% Output:
%       I -  filter result
%
% Meng Wu
% 2013

if length(n) == 1
    n = [n n];
end

if ismatrix( I )
    
    I = medfilt2(I, n);
    
elseif ndims( I ) == 3
    
    for iv = 1:size( I, 3)
        
        sinoView        = squeeze(I(:,:,iv));
        I(:,:,iv)    = medfilt2(sinoView, n, 'symmetric');
        
    end
    
end



end