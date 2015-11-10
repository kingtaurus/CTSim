function sinoOut = reverseSingoram( sinoIn )
% function sinoOut = reverseSingoram( sinoIn )
%
% Meng Wu at Stanford University
% 2013


sinoOut = zeros( size( sinoIn ), 'single' );

% rotate the sinogram
if ndims( sinoIn ) == 3
    
    noviews = size(sinoIn, 3);
    for iv = 1 : noviews
        sinoOut(:,:,iv) =  sinoIn(:,:,noviews-iv+1);
    end
    
else % 1D detector
    
    noviews = size(sinoIn, 2);
    for iv = 1 : noviews
        sinoOut(:, iv) =  sinoIn(:, noviews-iv+1);
    end
    
end

end