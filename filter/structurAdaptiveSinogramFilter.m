function sinoOut = structurAdaptiveSinogramFilter( sinoAtt, geom, w, s  )

fprintf('Structur adaptive sinogram filter ... \n');

geom.reconSize(3) = geom.reconSize(3) + 4;

ImgPrior = reconFBP( sinoAtt, geom, 'hamming');

ImgPrior = extendFOV( ImgPrior, geom );

sinoPrior =  forwardProjectMex( ImgPrior, geom );

sinoOut = sinoAtt;

for iv = 1 : size( sinoPrior, 3 )
    sinoOut( :, :, iv) = sasFilter( squeeze( sinoAtt( :, :, iv) ), squeeze( sinoPrior( :, :, iv) ), w,s );
end

fprintf('\t done.\n');


end


function uhat = sasFilter( u, v, w, s )

% get the size of the image
dim = size(u);

% compute the Gaussian kernel Ga
[X,Y] = meshgrid(-w:+w , -w:+w);
Ga    = exp(-(X.^2 + Y.^2 ) / (2*s^2));

% create s padded version of the noisy image
V = padarray(v , w * [1 1]);
U = padarray(u , w * [1 1]);

Ga = Ga / sum( Ga(:) );
% iterate over the pixels
uhat = zeros(size(u), class(u));
for x = w+1:dim(1)-w
    for y = w+1:dim(2)-w
        % get the neighborhood of (x,y)
        iMin = max(x-w,1);
        iMax = min(x+w,dim(1));
        jMin = max(y-w,1);
        jMax = min(y+w,dim(2));
        vv  = V(((x-w):(x+w))+w , ((y-w):(y+w))+w);
        uu = U(((x-w):(x+w))+w , ((y-w):(y+w))+w);
        % iterate over the pixels in the search window
        S = zeros(iMax-iMin+1 , jMax-jMin+1);
        
        for i = iMin:iMax
            for j = jMin:jMax
                
                % compute the similiarity between N(x,y) and N(i,j)
                VV = V(((i-w):(i+w))+w , ((j-w):(j+w))+w);
                S(i-iMin+1,j-jMin+1)    = sum( abs(vv(:) - VV(:)));
                
            end
        end
        
        
        % compute the similarity component of the weights
        S = 1 - S / sum( S(:) );
        
        % compute the (unnormalized) weights
        f = Ga(:) .* S(:);
        
        % compute the weighted average over the neighborhood
        uhat(x,y) = sum(f(:) .* uu(:)) / sum(f(:));
    end
end

end