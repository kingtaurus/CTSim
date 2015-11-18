function [sinoRaw, sinoPrime, sinoScatter] = estimateScatterFromCollimatedPixels( sinoRaw, ...
    sinoPrime, sinoScatter, veiwDownSampleRate, scatterCorrection  )
% function [sinoRaw, sinoPrime, sinoScatter] = combineScatterSimulation( sinoRaw, ...
%    sinoPrime, sinoScatter, veiwDownSampleRate  )
%
%
% Meng Wu at Stanford University
% 2014

if nargin < 5
    scatterCorrection = false;
end

fprintf('Adding Monte Carlo scatter signals:')


% preprocess the scatter results
sinoScatter( sinoScatter < 0 ) = 0;
for iv = 1:size( sinoScatter, 3)
   
    
    sinoScatter(:,:,iv) = correctScatterImageArtifacts( squeeze(sinoScatter(:,:,iv)) );
    
    sinoPrime(:,:,iv) = correctScatterImageArtifacts( squeeze(sinoPrime(:,:,iv)) );
    
end

% re-arrange
sinoPrime = rotateSinogram( sinoPrime, 0, 1 );
sinoPrime = reverseSingoram( sinoPrime );

sinoScatter = rotateSinogram( sinoScatter, 0, 1 );
sinoScatter = reverseSingoram( sinoScatter );

tSubStep = tic;

% view by view scaling
for iv = 1 : size( sinoScatter, 3) - 1
    
    % print progress figure
    if mod(iv, 50) == 0 || iv == 1
        progressStr = sprintf('at (%i/%i)... ',  iv, size( sinoScatter, 3) );
        fprintf(progressStr);
    end
    
    for isubv = 1 : veiwDownSampleRate
        
        a = (isubv - 1) / veiwDownSampleRate;
        
        % get the corresponding frames
        viewPrime   = (1 - a) * sinoPrime(:,:,iv) + a * sinoPrime(:,:,iv + 1) ;
        viewScatter = (1 - a) * sinoScatter(:,:,iv) + a * sinoScatter(:,:,iv + 1) ;
        viewRaw = sinoRaw( :, :, (iv-1)* veiwDownSampleRate + isubv );
        
        % find the scale
        if size( viewRaw, 1 ) == size( viewPrime, 1 )
            scale = mean( viewRaw ) / mean( viewPrime );
        end
        
        if isubv == 1
            sinoPrime(:,:,iv)   = viewPrime * scale;
            sinoScatter(:,:,iv) = viewScatter * scale;
        end
        
        if scatterCorrection
            sinoRaw( :, :, (iv-1)* veiwDownSampleRate + isubv ) = viewRaw + poissrnd( viewScatter * scale ) - viewScatter * scale;
        else
            sinoRaw( :, :, (iv-1)* veiwDownSampleRate + isubv ) = viewRaw + poissrnd( viewScatter * scale );
        end
        
    end
    
end

fprintf('(%is)\n', round(toc(tSubStep)));

end


function view = correctScatterImageArtifacts( view )

h = fspecial('Gaussian', 8, 2);
vpad = 4;
upad = 16;
cpad = 4;

% reduce extreme values
view = medfilt2(view, [5 5], 'symmetric');

% correct edge
view(:, 1:vpad ) = repmat( mean( view(:, vpad+1:2*vpad) , 2), [1 vpad]);
view(:, end-vpad+1:end ) = repmat( mean( view(:, end-2*vpad+1:end-vpad), 2 ), [1, vpad] );

view( 1:upad, : ) = repmat( mean( view(upad+1:2*upad, :)), [upad 1]);
view( end-upad+1:end, : ) = repmat( mean( view(end-2*upad+1:end-upad, :)), [upad 1]);

view( end/2-cpad+1:end/2+cpad, : ) = repmat( mean( view(end/2-2*cpad:end/2-cpad, :) + view(end/2+cpad:end/2+2*cpad, :) )/2, [cpad*2 1]);

view  = imfilter(view, h, 'replicate');

view( view < 0 ) = 0;
end

