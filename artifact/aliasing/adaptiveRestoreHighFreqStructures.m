function img = adaptiveRestoreHighFreqStructures( img, sinoSparse, geomSparse, T, rate, nitn, threshold, geomUp   )

detDSR = [2 4];
t = T;
srchWidth = 8;

for i = 1: nitn
    
    sinoResd = sinoSparse - forwardProjectMex( img, geomSparse  );
    
    sinoUp = upsamplingViewFeaturePreserving3( sinoResd, geomSparse, rate, srchWidth, detDSR, threshold );
    
    imgResidue = reconFBP( sinoUp, geomUp, 'hamming' );
    %
    
    
    if i < nitn
        imgResidue = immedian3( imgResidue, [3 3]);
        imgResidue = softThresholding(imgResidue, 0, t );
    end
    
    img = combineHighFreqStructures( img, imgResidue, geomSparse, true );
    
    img( img < 0 ) = 0;
    
    %figure(101); imdisp( squeeze( sinoUp(end/2,:,:) )', [-0.5 0.5]);
    figure(102); imdisp( imgResidue(:,:, 4:8:end ), [-0.01 0.01])
    figure(103); imdisp( img(:,:, 4:8:end ), [0.1 0.3])
    
    t = 0.5 * t;
end