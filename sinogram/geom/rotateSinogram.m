function sinoOut = rotateSinogram( sinoIn, rotate, transpose, flipLat, flipLng  )
% function sinoOut = rotateSinogram( sinoIn, rotate, transpose, flipLng, flipLat  )
%
% Meng Wu at Stanford University
% 2013

if nargin < 2
    rotate = 0;
end

if nargin < 3
    transpose = 0;
end

if nargin < 4
    flipLat = 0;
end

if nargin < 5
    flipLng = 0;
end


% rotate the sinogram
if ndims( sinoIn ) == 3
    
    [n, m, noviews] = size( sinoIn );
    
    if mod( rotate + transpose, 2 ) ==1
        sinoOut = zeros( [m n noviews], 'single');
    else
        sinoOut = zeros( [n m noviews], 'single');
    end
    
    for iv = 1 : noviews
        
        proj = sinoIn(:,:,iv);
        
        proj = rot90( proj, rotate );
        
        if transpose
            proj = proj' ;
        end
        
        if flipLat
            proj = fliplr( proj );
        end
        
        if flipLng
            proj = flipud( proj );
        end
        
        sinoOut(:,:,iv) = proj;
    end
    
else % 1D detector
    
    [n, noviews] = size( sinoIn );
    
    sinoOut = zeros( [n noviews], 'single');
    
    for iv = 1 : noviews
        
        proj = sinoIn(:, iv);
        
        if  mod( rotate, 4 ) == 2
            proj = rot90( proj, 2 );
        end
        
        sinoOut(:,iv) = proj;
    end
    
end

end



