%% BilateralFilter.m

%% BilateralFilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uhat = BilateralFilter(v,K,N)
% Apply a bilateral filter with window half-width w, distance standard
% deviation sd and similarity standard deviation ss to the image v
% Inputs:
%    v    : the input image
%    w    : the window half-width
%    sd   : the distance standard deviation
%    ss   : the similarity standard deviation
% Outputs:
%    uhat : the output image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uhat] = BilateralFilter(v,w,sd,ss)
    % get the size of the image
    dim = size(v);
    
    % compute the distance component of the weights
    [X,Y] = meshgrid(-w:+w , -w:+w);
    G = exp(-(X.^2 + Y.^2) / (2*sd^2));
    
    % iterate over the pixels
    uhat = zeros(size(v), class(v));
    for x = 1:dim(1)
        for y = 1:dim(2)
            % get the neighborhood of (x,y)
            iMin = max(x-w,1);
            iMax = min(x+w,dim(1));
            jMin = max(y-w,1);
            jMax = min(y+w,dim(2));
            Nxy  = v(iMin:iMax , jMin:jMax);

            % compute the similarity component of the weights
            H = exp(-(v(x,y) - Nxy).^2 / (2*ss^2));

            % compute the (unnormalized) weights
            B = H .* G((iMin:iMax)-x+w+1 , (jMin:jMax)-y+w+1);
            
            % compute the weighted average over the neighborhood
            uhat(x,y) = sum(B(:) .* Nxy(:)) / sum(B(:));
        end
    end
end