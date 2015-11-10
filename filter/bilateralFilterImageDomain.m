function imgOut = bilateralFilterImageDomain( imgIn, w, sd, ss )

fprintf('Bilateral filter in image domain with sd = %2.2f, ss = %2.2f ... \n', sd, ss);

imgOut = zeros( size(imgIn), class(imgIn) );

for iz = 1 : size( imgIn, 3 )
    
    view = squeeze( imgIn(:,:,iz) );
    imgOut(:,:,iz) = BilateralFilter( view, w, sd, ss );
    
end

fprintf('\t done.\n');

end

