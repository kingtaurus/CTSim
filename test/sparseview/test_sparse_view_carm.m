
% dataHeader = 'PigPBV1.m';
dataHeader = 'NeuroDataPatient4.m';

[ sinoRaw, sinoAtt, geom, spectrum ] = loadProcessedCarmCTData( dataHeader, 1 );

sinoAtt = rotateSinogram( sinoAtt, 0, 0, 1, 0 );

%sinoAtt =  reverseSingoram( sinoAtt );

sinoAtt = subtractSinoOffset( sinoAtt );

[ sinoAtt, geom ] = truncationCorrectionWaterCylinderFitting( sinoAtt, geom, 0.2, 4, 16 );
figure; imdisp( sinoAtt(end/2,:,:) );

%%
% Stanford FDK reconstruction

img_fbp = reconFBP( sinoAtt, geom );

img_fbp_r = rotateSinogram( img_fbp, 3 );
figure; imdisp( img_fbp_r(:,:,8:16:end), [0. 0.4] );

%% Downsample projection data to mimic sparse view case

[geomSparse, sinoSparse ]  = convertSinogram2SparseViewScan( geom, sinoAtt, 2  );

img_sparse_fbp = reconFBP( sinoSparse, geomSparse );

img_sparse_fbp_r = rotateSinogram( img_sparse_fbp, 3 );
figure; imdisp( img_sparse_fbp_r(:,:,8:16:end), [0.11 0.31] );


%% linear interpolation method
sinoInterp = upsamplingViews( sinoSparse, geomSparse, 2 );

img_sparse_fbp_interp = reconFBP( sinoInterp, geom, 'hamming', 1, 1);

img_sparse_fbp_interp_r = rotateSinogram( img_sparse_fbp_interp, 3 );
figure; imdisp( img_sparse_fbp_interp_r(:,:,8:16:end), [0.11 0.31] );
title 'FBP linear interpolation reconstruction (sparse views by4)';

%% 1.  The fist multi-resoltion piecewise linear frequency response FBP

img_fbp_ms = reconFBPMultiResolution( sinoSparse, geomSparse, 'hamming', 8, 0.6, 0);

img_fbp_ms_r = rotateSinogram( img_fbp_ms, 3 );
figure; imdisp( img_fbp_ms_r(:,:,8:16:end),  [0.11 0.31] );

%% 2. Extract high frequency structures in the sinogram

sino_low_freq = forwardProjectMex( img_fbp_ms, geomSparse  );

figure; imdisp( squeeze( sino_low_freq(end/2,:,:))' );
title 'Projections of low frequency structures';

% Subtract the non-aliasing component out of the true sinogram

sinoResd = sinoSparse - sino_low_freq;

figure; imdisp( squeeze( sinoResd(end/2,:,:))');
title 'Projections of structures that may cause view aliasing artifacts';

%% 3. View interpolation that also perserve the shapes

% The trick of the structure preserving interpolation is directly shifting
% and copying to the center the structures that are similar frow two
% neighboring views.
threshold = 0.8;
upSamplingRate = 2;
srchWidth = 6;
detDSR = [2 4];

sinoUp = upsamplingViewFeaturePreserving3( sinoResd, geomSparse, upSamplingRate, srchWidth, detDSR, threshold );

% The movement of the high frequency structurs in the singram now looks smooth.
figure; imdisp( squeeze(sinoUp(end/2,:,:))' );
title 'Interpolated high freqency structures';

% 
img_residue = reconFBP(  sinoUp, geom, 'hamming', 1, 1);

img_residue_r = rotateSinogram( img_residue, 3 );
figure; imdisp( img_residue_r(:,:,8:16:end), [-0.1 0.1] );
title 'Restored high frequency component in the image space';

% combine to get final result
img_fbp_new = combineHighFreqStructures( img_fbp_ms, img_residue, geom, true );

img_fbp_new_r = rotateSinogram( img_fbp_new, 3 );
figure; imdisp( img_fbp_new_r(:,:,8:16:end), [0.11 0.31] );
title 'Proposed sparse view FBP reconstruction';

%%

img_new_itr = adaptiveRestoreHighFreqStructures( img_fbp_ms, sinoSparse, geomSparse, 0.005, 2, 3, 0.8, geom );

img_new_itr_r = rotateSinogram( img_new_itr, 3 );
figure; imdisp( img_new_itr_r(:,:,8:16:end), [0.11 0.31] );
title 'Proposed sparse view FBP reconstruction';

figure; imdisp( img_new_itr_r(:,:,8:16:end) - img_fbp_r(:,:,8:16:end), [-0.1 0.1] );
title 'Restored high frequency component in the image space';

return;
%%
i = 84;

imgs = zeros( 360, 360, 4 );

imgs(:,:, 1) = img_fbp_r(:,:,i);
imgs(:,:, 2) = img_sparse_fbp_r(:,:,i);
imgs(:,:, 3) = img_sparse_fbp_interp_r(:,:,i);
imgs(:,:, 4) = img_fbp_new_r(:,:,i);

figure; imdisp( imgs(151:270, 121:240 ), [0.1 0.3] );

a = img_fbp_r( 151:270, 121:240 ,i);
std(  a(:) ) * 2000

a = img_sparse_fbp_r( 151:270, 121:240 ,i);
std(  a(:) ) * 2000

a = img_sparse_fbp_interp_r( 151:270, 121:240 ,i);
std(  a(:) ) * 2000

a = img_fbp_new_r( 151:270, 121:240 ,i);
std(  a(:) ) * 2000

%%

i = 64;

figure; imdisp( img_fbp_r(171:270, 121:220 ,i) );

w = fspecial('gaussian', 11, 1.5);
L = 0.6;

ssim_index(img_fbp_r(171:270, 121:220 ,i) , img_sparse_fbp_r(171:270, 121:220 ,i), [0.01 0.03], w, L)

ssim_index(img_fbp_r(171:270, 121:220,i), img_sparse_fbp_interp_r(171:270, 121:220 ,i), [0.01 0.03], w, L)

ssim_index(img_fbp_r(171:270, 121:220,i), img_fbp_new_r(171:270, 121:220 ,i), [0.01 0.03], w, L)

%%

imgs = zeros( 360, 360, 4 );
    
figure
filename = 'neuro4.gif';
for n = 1 : 128
   

imgs(:,:, 1) = img_fbp_r(:,:,n);
imgs(:,:, 2) = img_sparse_fbp_r(:,:,n);
imgs(:,:, 3) = img_sparse_fbp_interp_r(:,:,n);
imgs(:,:, 4) = img_fbp_new_r(:,:,n);

imdisp( imgs, [0.1 0.3] );
    
drawnow
frame = getframe(1);
im = frame2im(frame);
[A,map] = rgb2ind(im,256); 
	if n == 1;
		imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.2);
	else
		imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.2);
	end
end



