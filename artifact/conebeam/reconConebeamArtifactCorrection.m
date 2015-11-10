function img = reconConebeamArtifactCorrection( sinoAtt, sinoSd, geom, reconFunc, nitn, shrinkSize )



img = feval( reconFunc, sinoAtt );


for itn = 1 : nitn
    
    % compuate the sinogram error
    sinoDiff = sinoAtt - forwardProjectMex( img, geom );
    
    % shrinkage the noise
    sinoDiff = shrinkage( sinoDiff,  shrinkSize(1) );

    
    %sinoDiff( sinoDiff >  shrinkSize(2) ) =  shrinkSize(2); 
    %sinoDiff( sinoDiff < -shrinkSize(2) ) = -shrinkSize(2); 
    
    
    % add reconstructed error
    img = img - feval( reconFunc, sinoDiff ) / 2;
    
    
end

end


function a = shrinkage(a, b)

a = sign(a) .* max( abs(a) - b, 0 );

end