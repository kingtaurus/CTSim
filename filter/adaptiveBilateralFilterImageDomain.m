function imgOut = adaptiveBilateralFilterImageDomain( imgIn, imgVar, w, sd, beta )

fprintf('Bilateral filter in image domain with sd = %2.2f, beta = %2.2f ... \n', sd, beta);

imgOut = zeros( size(imgIn), class(imgIn) );

for iz = 1 : size( imgIn, 3 )
    
    imgOut(:,:,iz) = adaptiveBilateralFilter( squeeze( imgIn(:,:,iz) ), squeeze( imgVar(:,:,iz) ),  w, sd, beta );
    
end

fprintf('\t done.\n');

end

function [J] = adaptiveBilateralFilter(I, V, w,sd,beta)
    % get the size of the image
    dim = size(I);
    
    % compute the distance component of the weights
    [X,Y] = meshgrid(-w:+w , -w:+w);
    G = exp(-(X.^2 + Y.^2) / (2*sd^2));
    
    % iterate over the pixels
    J = zeros(size(I), class(I));
    for x = 1:dim(1)
        for y = 1:dim(2)
            % get the neighborhood of (x,y)
            iMin = max(x-w,1);
            iMax = min(x+w,dim(1));
            jMin = max(y-w,1);
            jMax = min(y+w,dim(2));
            Nxy  = I(iMin:iMax , jMin:jMax);

            % compute the similarity component of the weights
            H = exp(-(I(x,y) - Nxy).^2 / ( 2 * beta * V(x,y) ));

            % compute the (unnormalized) weights
            B = H .* G((iMin:iMax)-x+w+1 , (jMin:jMax)-y+w+1);
            
            % compute the weighted average over the neighborhood
            J(x,y) = sum(B(:) .* Nxy(:)) / sum(B(:));
        end
    end
end