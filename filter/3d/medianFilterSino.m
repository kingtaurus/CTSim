function sino = medianFilterSino( sino, n )
% Filter the sinogram with 2d median filter
% Input:
%       sino -  sinogram to process
%       n    -  half size of median filter
% Output:
%       sino -  filter result
%
% Meng Wu
% 2013

if length(n) == 1
    n = [n n];
end

if ismatrix( sino )
    
    sino = medfilt2(sino, n);
    
elseif ndims( sino ) == 3
    
    for iv = 1:size( sino, 3)
        
        sinoView        = squeeze(sino(:,:,iv));
        sino(:,:,iv)    = medfilt2(sinoView, n, 'symmetric');
        
    end
    
end



end