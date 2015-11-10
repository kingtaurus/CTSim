%% EGGFilter.m

%% EGGFilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uhat = EGGFilter(v,ws,wn,s,a)
% Apply an edge-preserving generalized Gaussian filter with search window 
% half-width ws, neighborhood window half-width wn, averaging weight standard
% deviation s and similarity norm standard deviation a to the image v
% Inputs:
%    v    : the input image
%    ws   : the search window half-width
%    wn   : the neighborhood window half-width
%    s    : the averaging weight standard deviation
%    a    : the similarity norm standard deviation
%    b    : the spatial standard deviation
% Outputs:
%    uhat : the output image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uhat, ex, ey] = EGGFilter(v, ws, wn, s, a, b)
% get the size of the image
n = length(v);

[Ex, Ey] = getEdges( v );
X = -wn:wn;

% create a padded version of the noisy image
V = padarray(v , wn * [1 1]);

% iterate over the pixels
uhat = zeros(n);

for x = 1:n
    for y = 1:n
        % compute the search window
        iMin = max(x-ws,1);
        iMax = min(x+ws,n);
        jMin = max(y-ws,1);
        jMax = min(y+ws,n);
        ex = Ex(x,y);
        ey = Ey(x,y);
        
        Ns   = v(iMin:iMax , jMin:jMax);
        
        gx = exp( - ex * X.^2  / (2 * a^2) );
        gy = exp( - ey * X'.^2 / (2 * a^2) );
        G = gy * gx;
        G = G / sum( G(:).^2 );


        % iterate over the pixels in the search window
        W = zeros(iMax-iMin+1 , jMax-jMin+1);
        for i = iMin:iMax
           
            for j = jMin:jMax
                % compute the similiarity between N(x,y) and N(i,j)
                Nxy = V(((x-wn):(x+wn))+wn , ((y-wn):(y+wn))+wn);
                Nij = V(((i-wn):(i+wn))+wn , ((j-wn):(j+wn))+wn);
                S   = sqrt(sum(G(:) .* (Nxy(:) - Nij(:)).^2)) ;
                
                % compute the (unnormalized) weight of pixel (i,j)
                W(i-iMin+1,j-jMin+1) = exp( - S / (2*s^2) ) ...
                    * exp( -( ey*(i-x)^2 + ey*(j-y)^2) / (2*b^2) );
            end
        end
        % compute the weighted average over the search neighborhood
        uhat(x,y) = sum(W(:) .* Ns(:)) / sum(W(:));        
    end
end

end

%%
function [ex, ey] = getEdges( u )

H = fspecial('prewitt');

ex = imfilter(u , H');
ex = abs( ex ) / 60 ;
ex(ex < 1 ) = 1 ;

ey = imfilter(u, H);
ey = abs(ey) / 60 ;
ey(ey < 1 ) = 1;

end
