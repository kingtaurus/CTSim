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

geom = loadProjectionGeometryShortScan( p, 223, 90 );

% geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e6);

% compute sinogram
[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDir, 1, 1, 2 );

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);

weights = computeWeightsPwls( sinoRaw, 0, spectrum.electronicNoise );


%% Reconstruction without view downsampling 

% Filter backprojection with hamming window

img_fbp = reconFBP( sinoAtt, geom, 'hamming', 1, 1);

figure; imshow( img_fbp(:,:,ceil(end/2)), map.windowAtt ); 
title 'FBP reconstruction hamming window (full views)';
evaluateSoftTissue( img_fbp, spectrum, map );

%% Downsample projection data to mimic sparse view case 
% We first consider downsample by 4 which is a moderate sparse view case,
% hard enough for interpolation and FBP reconstruction. The TV iterative
% reconstruction can deal it well.

rate = 10;

[geomSparse, sinoSparse, weightsSparse ]  = convertSinogram2SparseViewScan( geom, sinoAtt, rate, weights );

img_sparse_fbp = reconFBP( sinoSparse, geomSparse, 'hamming', 1, 1 );

figure; imshow( img_sparse_fbp(:,:,ceil(end/2)), map.windowAtt ); 
title 'FBP reconstruction (sparse views)';
evaluateSoftTissue( img_sparse_fbp, spectrum, map );


%% Linear interpolation method

[sinoInterp, geomInterp] = upsamplingViews( sinoSparse, geomSparse, rate );

img_sparse_fbp_interp = reconFBP( sinoInterp, geomInterp, 'hamming', 1, 1);

figure; imshow( img_sparse_fbp_interp(:,:,ceil(end/2)), map.windowAtt );
title 'FBP linear interpolation reconstruction';
evaluateSoftTissue( img_sparse_fbp_interp, spectrum, map );


%% Iterative reconstruction 
% TV regularized iterative reconstruction is very common method for sparse
% view CT reconstruction. This method usually require a stronger L1 based
% regularization to obtain a smooth image. The proper values of the turning
% parameter are often hard to determine. 
%
% [1] J. Tang, B. E. Nett, and G.-H. Chen, “Performance comparison between
% total variation (TV)-based compressed sensing and statistical iterative
% reconstruction algorithms.,” Phys. Med. Biol., vol. 54, no. 19, pp.
% 5781–5804, 2009.   
nitn = 100;
numos = 4;
img0 = img_sparse_fbp_interp;

img_sparse_tv_1 = reconPwlsLALMOs14( sinoSparse, weightsSparse, geomSparse, 100, 'isotv', nitn, 1, numos, img0  );

figure; imshow( img_sparse_tv_1(:,:,ceil(end/2)), map.windowAtt );
title 'PWLS with TV reconstruction beta = 100';
evaluateSoftTissue( img_sparse_tv_1, spectrum, map );

img_sparse_tv_2 = reconPwlsLALMOs14( sinoSparse, weightsSparse, geomSparse, 1000, 'isotv', nitn, 1, numos, img0  );

figure; imshow( img_sparse_tv_2(:,:,ceil(end/2)), map.windowAtt );
title 'PWLS with TV reconstruction beta = 200';
evaluateSoftTissue( img_sparse_tv_2, spectrum, map );

img_sparse_tv_3 = reconPwlsLALMOs14( sinoSparse, weightsSparse, geomSparse, 500, 'isotv', nitn, 1, numos, img0  );

figure; imshow( img_sparse_tv_3(:,:,ceil(end/2)), map.windowAtt );
title 'PWLS with TV reconstruction beta = 500';
evaluateSoftTissue( img_sparse_tv_3, spectrum, map );

img_sparse_tv_4 = reconPwlsLALMOs14( sinoSparse, weightsSparse, geomSparse, 1000, 'isotv', nitn, 1, numos, img0  );

figure; imshow( img_sparse_tv_4(:,:,ceil(end/2)), map.windowAtt );
title 'PWLS with TV reconstruction beta = 1000';
evaluateSoftTissue( img_sparse_tv_4, spectrum, map );


img_sparse_tv_5 = reconPwlsLALMOs14( sinoSparse, weightsSparse, geomSparse, 2000, 'isotv', nitn, 1, numos, img0  );

figure; imshow( img_sparse_tv_5(:,:,ceil(end/2)), map.windowAtt );
title 'PWLS with TV reconstruction beta = 2000';
evaluateSoftTissue( img_sparse_tv_5, spectrum, map );
return;


%% Failed attempted 2: Iterative reconstruction with interpolated view
% Use interploated projections to constrain the sparse view reconstruction,
% then one can reduce the strength of the TV regularization on the image
% domain.

[sinoInterp2, geomInterp4] = upsamplingViews( sinoSparse, geomSparse, 2 );
weightsInterp2 = upsamplingViews( weightsSparse, geomSparse, 2 );

img_interp_tv = reconPwlsLALMOs( sinoInterp4, weightsInterp2, geomInterp2, 200, 'isotv', nitn, 1, 4, img0  );

figure; imshow( img_interp_tv(:,:,ceil(end/2)), map.windowAtt );
title 'PWLS with TV reconstruction beta = 200 (interplated views)';
evaluateSoftTissue( img_interp_tv, spectrum, map );

% However, the TV based reconstruction often produces estimated projections
% at the missing angles that are more accurate than the interploation
% method. Thus the added constrains will introduce more errors in the
% reconstruction.     
% 
% This result has be shown in the 2014 CT meeing by Prof. Pan's group.


% Let us get forward projections at the missed angles for IR-TV
% reconstructions

sino_sparse_tv = forwardProjectMex( img_sparse_tv_2, geom );
sinoInterp = upsamplingViews( sinoSparse, geomSparse, rate );

% Show the differences from true sinogram
figure; imshow( squeeze( sinoAtt(end/2,:,:) - sinoInterp(end/2,:,:))', [-0.2 0.2] );
title 'Difference between true sinogram and linear interpolated sinogram';

figure; imshow( squeeze( sinoAtt(end/2,:,:) - sino_sparse_tv(end/2,:,:))', [-0.2 0.2] );
title 'Difference between true sinogram and TV estimated sinogram';


%% Failed attempted 2
% use random subsampling to mimic true compressed sensing, then average 
% the reconstructions with different subsampling to get the final image.
% The approach does not work well. There are several following
% observations:
%   1. The subsampled reconstruction with 0.9 rate is very similar to the
%   TV constrained reconstruction
%   2. As the sample rate goes down the reconstruction becomes very
%   unstable.
%   3. The approach get even more unstable to the number of the veiws are
%   too small. 
%   Therefore, in CT reconstruction, if you got some data, use all of them

beta = 200;
notest = 25;
nitns = 25;

delta = 0.9; %sub random sample rate
img_sparse_tv_90 = reconPwlsSeNesterovSqsResample( sinoSparse, weightsSparse, geomSparse, beta, 'isotv', nitns, delta, 2, img0  );
for i = 1:notest-1 
img_sparse_tv_90 = img_sparse_tv_90 + reconPwlsSeNesterovSqsResample( sinoSparse, weightsSparse, geomSparse, beta, 'isotv', nitns, delta, 2, img0  );
end
% average
img_sparse_tv_90 = img_sparse_tv_90 / notest;

figure; imshow( img_sparse_tv_90(:,:,ceil(end/2)), map.windowAtt );
title 'PWLS with TV reconstruction 90% (sparse views by 4)';
evaluateSoftTissue( img_sparse_tv_90, spectrum, map );

delta = 0.8; %sub random sample rate
img_sparse_tv_80 = reconPwlsSeNesterovSqsResample( sinoSparse, weightsSparse, geomSparse, beta, 'isotv', nitns, delta, 2, img0  );
for i = 1:notest-1 
img_sparse_tv_80 = img_sparse_tv_80 + reconPwlsSeNesterovSqsResample( sinoSparse, weightsSparse, geomSparse, beta, 'isotv', nitns, delta, 2, img0  );
end
% average
img_sparse_tv_80 = img_sparse_tv_80 / notest;

figure; imshow( img_sparse_tv_80(:,:,ceil(end/2)), map.windowAtt );
title 'PWLS with TV reconstruction 80% (sparse views by 4)';
evaluateSoftTissue( img_sparse_tv_80, spectrum, map );

delta = 0.7; %sub random sample rate
img_sparse_tv_70 = reconPwlsSeNesterovSqsResample( sinoSparse, weightsSparse, geomSparse, beta, 'isotv', nitns, delta, 2, img0  );
for i = 1:notest-1 
img_sparse_tv_70 = img_sparse_tv_70 + reconPwlsSeNesterovSqsResample( sinoSparse, weightsSparse, geomSparse, beta, 'isotv', nitns, delta, 2, img0  );
end
% average
img_sparse_tv_70 = img_sparse_tv_70 / notest;

figure; imshow( img_sparse_tv_70(:,:,ceil(end/2)), map.windowAtt );
title 'PWLS with TV reconstruction 80% (sparse views by 4)';
evaluateSoftTissue( img_sparse_tv_70, spectrum, map );

return;

%% Iterative reconstruction (penalized weighted least square)
nitn = 50;
numos = 16;

% Huber penalty function
img_huber = reconPwlsSeNesterovSqs( sinoAtt, weights, geom, 100, 'isotv', nitn, 1e-3, numos );

figure; imshow( img_huber(:,:,ceil(end/2)), map.windowAtt ); 
title 'PWLS with Huber reconstruction (full views)';
evaluateSoftTissue( img_huber, spectrum, map );

return;

