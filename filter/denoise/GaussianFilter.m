%% GaussianFilter.m

%% GaussianFilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uhat = GaussianFilter(v,w,s)
% Apply a Gaussian filter with window half-width w and standard
% deviation s to the image v
% Inputs:
%    v    : the input image
%    w    : the half-width of the filter window
%    s    : the standard deviation of the Gaussian kernel
% Outputs:
%    uhat : the output image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uhat] = GaussianFilter(v,w,s)
    % construct the Gaussian filter
    G = fspecial('gaussian' , 2*w+1 , s);
    
    % apply the Gaussian filter
    uhat = imfilter(v , G);
end