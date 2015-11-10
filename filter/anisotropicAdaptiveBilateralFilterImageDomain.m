function imgOut = anisotropicAdaptiveBilateralFilterImageDomain( imgIn, imgVarH, imgVarV, w, sd, ss,  beta )

fprintf('Anisotropic adaptive bilateral filter in image domain with sd = %2.2f, ss = %2.2f, beta = %2.2f ... \n', sd, ss, beta);

imgOut = zeros( size(imgIn), class(imgIn) );

for iz = 1 : size( imgIn, 3 )
    imgOut(:,:,iz) = anisotropicAdaptiveBilateralFilter( squeeze( imgIn(:,:,iz) ), ...
        squeeze( imgVarH(:,:,iz) ), squeeze( imgVarV(:,:,iz) ),  w, sd, ss, beta );
end

fprintf('\t done.\n');

end

function [J] = anisotropicAdaptiveBilateralFilter(I, Vh, Vv, w,sd, ss, beta)
    % get the size of the image
    dim = size(I);
    
    % compute the distance component of the weights
    [X,Y] = meshgrid(-w:+w , -w:+w);
    G = exp(-(X.^2 + Y.^2) / (2*sd^2));
    [Ah, Av] = estimateCorr( w );

    % iterate over the pixels
    J = zeros(size(I), class(I));
    for x = 1+w:dim(1)-w
        for y = 1+w:dim(2)-w
            % get the neighborhood of (x,y)
            iMin = max(x-w,1);
            iMax = min(x+w,dim(1));
            jMin = max(y-w,1);
            jMax = min(y+w,dim(2));
            Nxy  = I(iMin:iMax , jMin:jMax);

            % compute the similarity component of the weights
            H = exp(- ( (I(x,y) - Nxy).^2  / ( 2 * ss^2) ) );

            H = H .* exp(- ( Ah * Vh(x,y) + Av * Vv(x, y)) / ( 2 * beta * (Vh(x,y) + Vv(x,y))) );
            % compute the (unnormalized) weights
            B = H .* G;
            
            % compute the weighted average over the neighborhood
            J(x,y) = sum(B(:) .* Nxy(:)) / sum(B(:));
        end
    end
end


function [Ax, Ay] = estimateCorr( w )


[X,Y] = meshgrid(-w:+w , -w:+w);

R2 = sqrt( X.^2 + Y.^2 );
R2( w + 1, w + 1 ) = 1;

Ax = abs(  X ./ R2 );
Ay = abs(  Y ./ R2 );



end

