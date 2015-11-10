%% NLMeansFilter.m

%% NLMeansFilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uhat = NLMeansFilter(v,ws,wn,h,s)
% Apply s nonlocal means filter with search window half-width ws,
% neighborhood window half-width wn, averaging weight standard
% deviation h and similarity norm standard deviation s to the image v
% Inputs:
%    v    : the input image
%    ws   : the search window half-width
%    wn   : the neighborhood window half-width
%    s    : the averaging weight standard deviation
%    h    : the similarity norm standard deviation
% Outputs:
%    uhat : the output image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uhat] = NLMeansFilter(v,ws,wn,s,h,map)

    if nargin < 6
        map = ones( size(v), 'logical');
    end

    % get the size of the image
    n = length(v);
    
    % create s padded version of the noisy image
    V = padarray(v , wn * [1 1]);
    
    % compute the Gaussian kernel Ga
    [X,Y] = meshgrid(-wn:+wn , -wn:+wn);
    Ga    = exp(-(X.^2 + Y.^2) / (2*s^2));

    Ga = Ga / sum( Ga(:) );
    % iterate over the pixels
    uhat = zeros(n);
    for x = 1:size(v,1)
        for y = 1:size(v,2)
            
            if map(x,y) == 0
                uhat(x, y) = v(x , y);
                continue;
            end
            
            % compute the search window
            iMin = max(x-ws,1);
            iMax = min(x+ws,size(v,1));
            jMin = max(y-ws,1);
            jMax = min(y+ws,size(v,2));
            Ns   = v(iMin:iMax , jMin:jMax);
            Nxy = V(((x-wn):(x+wn))+wn , ((y-wn):(y+wn))+wn);
            % iterate over the pixels in the search window
            W = zeros(iMax-iMin+1 , jMax-jMin+1);
            
            for i = iMin:iMax
                for j = jMin:jMax
                    % compute the similiarity between N(x,y) and N(i,j)
                    
                    Nij = V(((i-wn):(i+wn))+wn , ((j-wn):(j+wn))+wn);
                    S   = sum(Ga(:) .* (Nxy(:) - Nij(:)).^2);
                   
                    % compute the (unnormalized) weight of pixel (i,j)
                    W(i-iMin+1,j-jMin+1) = exp(-S / (2*h^2));
                end
            end
            % compute the weighted average over the search neighborhood
            uhat(x,y) = sum(W(:) .* Ns(:)) / sum(W(:));
        end
    end
end