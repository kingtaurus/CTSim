function imgOut = nonLocalMeanFilterImageDomain( imgIn, ws, wn, s, h )

fprintf('Non local mean filter in image domain with s = %2.2f, h = %2.2f ... \n', s, h);


imgOut = zeros( size(imgIn), class(imgIn) );

for iz = 1 : size( imgIn, 3 )
    
    view = squeeze( imgIn(:,:,iz) );
    
    map = view > -500;
    
    imgOut(:,:,iz) = NLMeansFilter( view, ws, wn, s, h, map );
    
end

fprintf('\t done.\n');

end