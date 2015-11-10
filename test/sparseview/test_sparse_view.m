%% test_sparse_view.m
% Compare several reconstruction methods for sparse view CT
% Meng Wu at University of Erlangen-Nuremburg
% 2014.10

%% Simulate normal CT data
clear;
close all;
load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadXCATPhantom(p);

geom = loadProjectionGeometryShortScan( p, 210, 0 );

% geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e6);

% compute sinogram
[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDir, 1, 1, 2 );

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);

weights = computeWeightsPwls( sinoRaw, 0, spectrum.electronicNoise );

map.windowAtt = [0.05 0.35];


%% Reconstruction without view downsampling

% Filter backprojection with hamming window

img_fbp = reconFBP( sinoAtt, geom, 'ram-lak', 1, 1);

figure; imdisp( img_fbp(:,:,ceil(end/2)), map.windowAtt );
title 'FBP reconstruction hamming window (full views)';
evaluateSoftTissue( img_fbp, spectrum, map );


%% Downsample projection data to mimic sparse view case
% We first consider downsample by 4 which is a moderate sparse view case,
% hard enough for interpolation and FBP reconstruction. The TV iterative
% reconstruction can deal it well.

rate = 6;
[geomSparse, sinoSparse, weightsSparse ]  = convertSinogram2SparseViewScan( geom, sinoAtt, rate, weights );

img_sparse_fbp = reconFBP( sinoSparse, geomSparse, 'hamming', 1, 1);

figure; imdisp( img_sparse_fbp(:,:,ceil(end/2)), map.windowAtt );
title 'FBP reconstruction (sparse views by4)';
evaluateSoftTissue( img_sparse_fbp, spectrum, map );


% linear interpolation method
[sinoInterp, geomInterp] = upsamplingViews( sinoSparse, geomSparse, rate );

img_sparse_fbp_interp = reconFBP( sinoInterp, geomInterp, 'hamming', 1, 1);

figure; imdisp( img_sparse_fbp_interp(:,:,ceil(end/2)), map.windowAtt );
title 'FBP linear interpolation reconstruction (sparse views by4)';
evaluateSoftTissue( img_sparse_fbp_interp, spectrum, map );



%% 2-Pass FBP based view undersampling artifact correction algorithms
% Proposed by Meng Wu
%
% The scans with not enough projections can caused by limited detector
% readout rate, or performed intentionally for dose reduction.
% Undersampling in projections with total variation constrained iterative
% reconstructions has shown ability of provide good image quality, but
% usually have long computation time [1]. Analytical methods based on
% interpolation in projection space have been proposed too [2-3]. The cause
% of the view aliasing artifacts in the CT reconstruction is some high
% frequency structures are not supported by the angular sampling. Although
% the interpolation method are often used to increase angular sampling
% rate, the high frequency structures are not well preserved in the
% interpolation. There are method call directional sinogram interpolation
% method to generate more projections by optimized (iterative)
% double-orientation estimation in sinogram space and directional
% interpolation. The challenge of this approach is the high frequency
% structure may hide inside the line integral (the total attenuation).
% Therefore, we want to take out the part that not causing aliasing out of
% the sinogram before the interpolation.
%
% H. Zhang and J.-J. Sonke, “Directional sinogram interpolation for
% sparse angular acquisition in cone-beam computed tomography.,?J. Xray.
% Sci. Technol., vol. 21, no. 4, pp. 481?6, Jan. 2013.

%% 1.  The fist multi-resoltion piecewise linear frequency response FBP
%
% The fist task it try to reconstruction as image that dose not has
% asliasing. In other word, the image will have locally highest resoltion
% that is support by the scan view sampling. The image will be
% multi-resoltuion. Closer to the center of grantry rotation have hight
% resoltuion. The reconstruction as such image, we spilt the ramp filter
% into different components in the infrequency segment.

img_fbp_ms = reconFBPMultiResolution( sinoSparse, geomSparse, 'hamming', 8, 0.5, 0 );

figure; imdisp( img_fbp_ms(:,:,ceil(end/2)), map.windowAtt );
title 'Multi-resoltionFBP reconstruction (sparse views by4)';
evaluateSoftTissue( img_fbp_ms, spectrum, map );

%%
% img_blur = azimuthalBluring( img_sparse_fbp_interp, geom, 8 );
%
% imdisp( img_blur, map.windowAtt ) ;
%
% imdisp( img_blur - img_sparse_fbp_interp ) ;
%
%
% sino_low_freq = forwardProjectMex( img_fbp_ms, geomInterp );
%
% sinoUp_linear = squeeze(sinoInterp(end/2,:,:))' - squeeze(sino_low_freq(end/2,:,:))';
%
% figure; imdisp( sinoUp_linear, [-0.5 0.5] );
% title 'Linearly interpolated high freqency structures';
%
%
% img_residue_linear = reconFBP( sinoInterp - sino_low_freq, geomInterp, 'hamming', 1, 1);
% figure; imshow( img_residue_linear(:,:,ceil(end/2)), [-0.1 0.1] );
% title 'Linear interpolation Restored high frequency component in the image space';

%% 2. Extract high frequency structures in the sinogram
% Forwar project the non-aliasing FBP reconstruction

sino_low_freq = forwardProjectMex( img_fbp_ms, geomSparse  );

% Subtract the non-aliasing component out of the true sinogram

sinoResd = sinoSparse - sino_low_freq;

figure; imdisp( squeeze( sinoResd(end/2,:,:))', [-0.5 0.5] );
title 'Projections of structures that may cause view aliasing artifacts';


%% 3. View interpolation that also perserve the shapes

% The trick of the structure preserving interpolation is directly shifting
% and copying to the center the structures that are similar frow two
% neighboring views.
threshold = 1.0;
upSamplingRate = rate;
srchWidth = rate ;
detDSR = [2 4];

[sinoUp, geomUp] = upsamplingViewFeaturePreserving3( sinoResd, geomSparse, upSamplingRate, srchWidth, detDSR, threshold );

% The movement of the high frequency structurs in the singram now looks smooth.
figure; imdisp( squeeze(sinoUp(end/2,:,:))', [-0.5 0.5] );
title 'Interpolated high freqency structures';

img_residue = reconFBP(  sinoUp, geomUp, 'hamming', 1, 1);

figure; imdisp( img_residue(:,:,ceil(end/2)), [-0.05 0.05] );
title 'Restored high frequency component in the image space';

% combine to get final result
% img_fbp_new = img_residue + img_fbp_ms;
img_fbp_new = combineHighFreqStructures( img_fbp_ms, img_residue, geom, false );

figure; imdisp( img_fbp_new(:,:,ceil(end/2)), map.windowAtt );
title 'Proposed sparse view FBP reconstruction';
evaluateSoftTissue( img_fbp_new, spectrum, map );


%%

img_new_itr = adaptiveRestoreHighFreqStructures( img_fbp_ms, sinoSparse, geomSparse, 0.02, rate, 3,  1, geomUp );

figure; imdisp( img_new_itr(:,:,ceil(end/2)), map.windowAtt );
title 'Proposed sparse view FBP reconstruction with 3 iterations';
evaluateSoftTissue( img_new_itr, spectrum, map );

figure; imdisp( img_new_itr(:,:,ceil(end/2)) - img_fbp_ms(:,:,ceil(end/2)), [-0.05 0.05] );
title 'Restored high frequency component in the image space';


return;

%%

close all;

impop( img_fbp, map.windowAtt);
impop( img_fbp_ms, map.windowAtt);
impop( img_sparse_fbp_interp, map.windowAtt);
impop( img_fbp_new, map.windowAtt);
impop( img_new_itr, map.windowAtt);

%%

img_error_interp = img_sparse_fbp_interp - img_fbp;
img_error_interp = img_error_interp(:,:,end/2);
mean( abs( img_error_interp(:) ))
impop( img_error_interp, [-0.05 0.05] );

img_error_new = img_fbp_new - img_fbp;
img_error_new = img_error_new(:,:,end/2);
impop( img_error_new, [-0.05 0.05]);
mean( abs( img_error_new(:) ))

img_err_itr = img_new_itr - img_fbp;
img_err_itr = img_err_itr(:,:,end/2);
impop( img_err_itr, [-0.05 0.05] ) ;
mean( abs( img_err_itr(:) ))

return

%%

a = squeeze( sinoUp(end/2,:,:) );
b = a(:, 31:33);
b(:, 1) = b(:,1) + 0.3 ;
b(:, 3) = b(:,3) - 0.3;
figure;  plot( b );
legend( 'Projection 1', 'Interpolated Projection', 'Projection 2' );


c = b;
c(:, 2) = ( c(:, 1) + c(:, 3) ) /2 ;
figure;  plot( c );
legend( 'Projection 1', 'Interpolated Projection', 'Projection 2' );
%% wavelet based combine approach
addpath( genpath('..\Tools'));
W = Wavelet;

%img_fbp_ms1 = reconFBP2dMultiResolution( sinoSparse, geomSparse, 'hamming', 8, 0.4 );
%img_fbp_ms2 = reconFBP2dMultiResolution( sinoSparse, geomSparse, 'hamming', 8, 0.5 );

% compute wavlet of the non-aliasing image and residue image
im_W = W * ( img_residue_linear );
im_R = W * img_residue;

% find the maximum 20 percents wavelet coefficients
m = sort(  abs( im_W(:) )  .* abs( im_R(:) ) ,'descend');
ndx = floor(length(m) * 0.05);
thresh1 = m(ndx);

% use selected wavelet cofficients to restor residue image
im_W_th = im_R .* (abs(im_W) > thresh1) ;
figure; imshow( abs( im_W_th ), [0 0.1] );
img_fbp_residue_wavelet = W' * im_W_th;

% combine togethon
img_fbp_new_wavelet = img_fbp_ms + img_fbp_residue_wavelet;

figure; imshow( img_fbp_residue_wavelet(:,:,ceil(end/2)), [-0.1 0.1] );
title 'Wavelwt restored residue image';

figure; imshow( img_fbp_new_wavelet(:,:,ceil(end/2)), map.windowAtt );
title 'Proposed sparse view FBP reconstruction 2 (sparse views by4)';
evaluateSoftTissue( img_fbp_new_wavelet, spectrum, map );

