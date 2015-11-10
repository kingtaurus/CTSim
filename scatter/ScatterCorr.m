function [scattercorr] = ScatterCorr( allscatter, ncorr)
% Simple scatter estimation:
% This is for the scatter estimation in numerical simulations, where 
% allscatter is the scatter data from the detector. Assume only certain 
% rows are in the exposed region and others are in the collimator shadow. 
% Use the scatter data in shadow to estimate the scatter data in the 
% exposed region.
%
% Input: 
%   allscatter  -  2D scaled scatter data
%   ncorr       -  number of rows (in z direction) used for scatter estimation
%
% Output:
%   scattercorr -  2D scatter estimation
%
% For example, if allscatter is 96*512, only (5:92) rows reprsent the
% primary + scatter data, and others represent scatter data.
%
% He Yang
% 1.29.2015

if (nargin < 2)
    ncorr = 4;
end

ufilt = 9; % The width of kernel for lateral smoothing
ufilt2 = floor(ufilt/2);


% Case II: use only two rows, i.e. the (ncorr)th and (n2-ncorr+1)th row
[n1,n2] = size(allscatter);
scattercorr = zeros(n1,n2);
scattercorr(1:ncorr,:) = allscatter(1:ncorr,:);  
scattercorr(n1-ncorr+1:n1,:) = allscatter(n1-ncorr+1:n1,:);

% 1D median filter for the edges
%x1 = medfilt1(allscatter(ncorr,:),10);
%x2 = medfilt1(allscatter(n1-ncorr+1,:),10);
x1 = allscatter(ncorr,:);
x2 = allscatter(n1-ncorr+1,:);


x = [ncorr, n1-ncorr+1]';    % known rows
xin = [(ncorr+1):(n1-ncorr)]'; % unknown rows
%x = [ 1:ncorr, (n2-ncorr+1):n2];
%xin = (ncorr+1):(n2-ncorr);
% scattercorr(:,1:ncorr) = medfilt1(allscatter(:,1:ncorr,1),10);
% scattercorr(:,n2-ncorr+1:n2) = medfilt1(allscatter(:,n2-ncorr+1:n2),10);

%tSubStep = tic;
% Preprocess the scatter results
%allscatter( allscatter < 0 ) = 0.;
%view = zeros(scattercorr);


    %allscatter(:,:,iv) = correctScatterImageArtifacts( squeeze(allscatter(:,:,iv)) );
    
% % Columnwise scatter estimation (sinograms have been rotated)
% for j = 1:n2
%     y = [x1(j) x2(j)]';
%     %y = scattercorr(x,j);
%     %pp = polyfit(x,y,1); % Linear fitting
%     %pp = polyfit(x,y,2); % Quadratic fitting
%     %scattercorr(xin,j) = polyval(pp,xin);
%     scattercorr(xin,j) = interp1(x,y,xin);
% end

% Columnwise scatter estimation2: use the smallest estimation
for j = 1:n2
    y = [x1(j) x2(j)]';
    scattercorr(xin,j) = interp1(x,y,xin);
    scattercorr(xin,j) = repmat(min(x1(j),x2(j)),[n1-2*ncorr 1]);
end

    
%view = scattercorr;

% 2D smoothing
%scattercorr(ncorr+1:n1-ncorr,:) = medfilt2(scattercorr(ncorr+1:n1-ncorr,:), [7 7]);

% % Lateral smoothing
% for j = 1:n2
%     if (j<=1+ufilt2)
%         scattercorr(xin,j) = mean( view(xin, 1:j+ufilt2),2 );
%     elseif (j>1+ufilt2)&&(j<n2-2*ufilt2)
%         scattercorr(xin,j) = mean( view(xin, j-ufilt2:j+ufilt2),2 );
%     else
%         scattercorr(xin,j) = mean( view(xin, j-ufilt2:n2),2 );
%     end
% end
%   


scattercorr( scattercorr < 0 ) = 0.;  





%fprintf('(%is)\n', round(toc(tSubStep)));

end


% function view = correctScatterImageArtifacts( view )
% 
% h = fspecial('Gaussian', 8, 2);
% vpad = 4;
% upad = 16;
% cpad = 4;
% 
% % reduce extreme values
% view = medfilt2(view, [5 5], 'symmetric');
% 
% % correct edge
% view(:, 1:vpad ) = repmat( mean( view(:, vpad+1:2*vpad) , 2), [1 vpad]);
% view(:, end-vpad+1:end ) = repmat( mean( view(:, end-2*vpad+1:end-vpad), 2 ), [1, vpad] );
% 
% view( 1:upad, : ) = repmat( mean( view(upad+1:2*upad, :)), [upad 1]);
% view( end-upad+1:end, : ) = repmat( mean( view(end-2*upad+1:end-upad, :)), [upad 1]);
% 
% view( end/2-cpad+1:end/2+cpad, : ) = repmat( mean( view(end/2-2*cpad:end/2-cpad, :) + view(end/2+cpad:end/2+2*cpad, :) )/2, [cpad*2 1]);
% 
% view  = imfilter(view, h, 'replicate');
% 
% view( view < 0 ) = 0;
% end

