function [sino, redundantPixel1, redundantPixel2] = sinogramGainLinearCorrectionRedundant( sino, geom, secondSinoIndex, degree )
% function [sino, redundantPixel1, redundantPixel2] = ...
%   sinogramGainLinearCorrectionRedundant( sino, geom, secondSinoIndex, degree )
% note:
%   degree - degree of polynomial fit
%
%
% Meng Wu @ Stanford Univercity
% 2014.2


fprintf( 'Sinogram gain linaer correction using redundant data : \n');

minbetas = min(geom.betas);

nu = geom.detSize(1);
du = geom.detSpacing(1);
ou = geom.detOffset(1);

noViews = geom.noViews;
SDD = geom.SDD;
SAD = geom.SAD;
betas = geom.betas - min(geom.betas);


%sino = binningDetecotRows( sino, round(  geom.reconSpacing(3) / geom.detSpacing(2)  ) );

% array of detector positions in mm
u = ( ( -(nu-1)/2:(nu-1)/2) + ou ) * du ;
gamma = atan( - u / SDD );
Tau  = abs(geom.betas(end) - geom.betas(1)) - pi;


weightMap = zeros( geom.detSize(1), noViews );
beta = zeros( geom.detSize(1), noViews );
theta = zeros( geom.detSize(1), noViews );

% label pixels in the fan beam geomerty
for iv = 1:noViews
    weightMap(:,iv) =  shortScanSliverWeight( gamma, betas(iv),  Tau);
    beta(:,iv) =  betas(iv);
    theta(:,iv) = + gamma ;
end


% find coresponding redundant pixels
redundantSino1 = weightMap < 0.9 & weightMap > 0.1;
redundantSino1(:,secondSinoIndex:end) = 0;
%redundantSino2 = weightMap < 0.8 & weightMap > 0.1;
%redundantSino2(:,1:secondSinoIndex-1) = 0;


beta2 = beta(redundantSino1) + pi - 2 * theta(redundantSino1 );
gamma2 = - theta(redundantSino1 );

[bb, gg] = meshgrid( betas(secondSinoIndex:end), gamma );


k = 5;
if ndims( sino ) == 3

    numberOfRedundentPixels = sum(redundantSino1(:));
    redundantPixel1 = zeros( numberOfRedundentPixels * k, 1);
    redundantPixel2 = zeros( numberOfRedundentPixels * k, 1);
    
    for i = 0:k-1
        centralRow = squeeze( sino(end/2 + i - (k-1)/2,:,:));
        redundantPixel1(i*numberOfRedundentPixels+1 : (i+1)*numberOfRedundentPixels) = centralRow( redundantSino1(:) );
        redundantPixel2(i*numberOfRedundentPixels+1 : (i+1)*numberOfRedundentPixels) = interp2( bb, gg, centralRow(:,secondSinoIndex:end), beta2, gamma2 );
        
    end
else
    redundantPixel1 = sino( redundantSino1(:) );
    redundantPixel2 = interp2( bb, gg, sino(:,secondSinoIndex:end), beta2, gamma2 );
    
end

validPixels = redundantPixel1 > 0  & redundantPixel2 > 0 & ~isnan( redundantPixel2 );
redundantPixel1 = redundantPixel1( validPixels );
redundantPixel2 = redundantPixel2( validPixels );


fprintf( '\tNumber of redundant pixels: %i \n', length(redundantPixel2));
fprintf( '\tPercentage of total pixels: %2.3f \n', length(redundantPixel2)/numel(sino) );

% figure;
% imagesc( weightMap ); xlabel 'beta', ylabel 'gamma'
%

n = length(redundantPixel1);
if degree == 1
    X = [redundantPixel2 ones(n, 1)  ];
    p = polyfit(redundantPixel2, redundantPixel1, 1);
    fprintf('\tModel 1: sino1 = %2.2f * sino2 + %2.2f \n', p(1), p(2) );
    
elseif degree == 2
    X = [redundantPixel2.^2 redundantPixel2 ones(n, 1)  ];
    p = polyfit(redundantPixel2, redundantPixel1, 2);
    fprintf('\tModel 2: sino1 = %2.2f * sino2 ^ 2 +  %2.2f * sino2 +  %2.2f \n', p(1), p(2), p(3)  );
elseif degree == 3
    X = [redundantPixel2.^3 redundantPixel2.^2 redundantPixel2 ones(n, 1)  ];
        p = polyfit(redundantPixel2, redundantPixel1, 3);
    fprintf('\tModel 3: sino1 = %2.2f * sino2 ^ 3 + %2.2f * sino2 ^ 2 +  %2.2f * sino2 +  %2.2f \n', p(1), p(2), p(3), p(3) );
else
    Error('Degree is too high')
end

Z = (X' * X )^(-1);
residues = redundantPixel1 - polyval(p, redundantPixel2);
sigma =  sqrt( sum( residues.^2 ) / (n-1));
for i = 1:degree+1
    z = p(i) / ( sigma * sqrt( Z(i,i) ) ) ;
    fprintf('\tp(%i) = %2.3f, \tt-test = %2.3f, \tp-val = %2.3e. \n', i,  p(i), z, 1-tcdf(abs(z), n-degree-1) );
end
fprintf('\tStandard Error = %2.3f, R-squared = %2.4f\n', sigma, 1 - sigma^2 / std(redundantPixel1) );

if ndims( sino ) == 3
    sino(:,:,secondSinoIndex:end) = polyval(p, sino(:,:,secondSinoIndex:end));
else
    sino(:,secondSinoIndex:end) = polyval(p, sino(:,secondSinoIndex:end));
end

% angular spacing gain correction
if length( size( squeeze( sino ) ) ) == 3
    sino(:,:,1 : geom.noViews1)     = sino(:,:,1 : geom.noViews1) * geom.dbeta1 / geom.dbetaAverage;
    sino(:,:,geom.noViews1+1:end )  = sino(:,:,geom.noViews1+1:end ) * geom.dbeta2 / geom.dbetaAverage;
else
    sino(:,1 : geom.noViews1)     = sino(:,1 : geom.noViews1) * geom.dbeta1 / geom.dbetaAverage;
    sino(:,geom.noViews1+1:end )  = sino(:,geom.noViews1+1:end ) * geom.dbeta2 / geom.dbetaAverage;
end


figure;
x = 0:0.1:max(redundantPixel2);
y = polyval(p, x );
plot(redundantPixel2, redundantPixel1, '.' ); hold on;
%plot(x, y, 'r', 'lineWidth', 2 );
xlabel 'MV', ylabel 'KV';



fprintf('Done!\n\n\n');

%figure; imdisp( squeeze( sino(end/2,:,:)));

end