% Test simple scatter estimation:
% This is for the scatter correction in numerical simulations, where 
% allscatter is the scatter data from the detector. Assume only certain 
% rows are in the exposed region and others are in the collimator shadow. 
% Use the scatter data in shadow to estimate the scatter data in the 
% exposed region.
%
% Input: 
%   allscatter  -  scaled scatter data obtained from Monte Carlo simulations
%   ncorr       -  number of rows (in z direction) used for scatter estimation
%
% For example, if allscatter is 512*96*8, only 512*(35:60)*8 reprsents the
% primary + scatter data, and others represent scatter data.
%
% He Yang
% 2.3.2015
% 
% if (nargin < 2)
%     ncorr = 5;
% end

clear all; close all;
%% Load scatter data all_scatter(1:512, 1:96, 1:72)
load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Pelvis5_Cone_corrected.mat');

ix = 1:512; iy=1:96;
[ixx, iyy] = meshgrid(iy,ix);
figure; mesh(ixx,iyy,all_scatter(:,:,1)); title('allscatter(:,:,1)');

figure; plot(all_scatter(:,1,1)); title('allscatter(:,1,1)');
figure; plot(all_scatter(:,end,1)); title('allscatter(:,end,1)');

x1 = medfilt1(all_scatter(:,1,1),40);
figure; plot(x1);


%% Correct scatter image artifacts
% Case I: start from the top and bottom rows
[n1,n2,n3] = size(all_scatter);
scattercorr = zeros(n1,n2);

ncorr = 3;

% 1D median filter for the edges
scattercorr(:,1:ncorr) = medfilt1(all_scatter(:,1:ncorr,1),40);
scattercorr(:,n2-ncorr+1:n2) = medfilt1(all_scatter(:,n2-ncorr+1:n2,1),40);
%scattercorr(:,1:ncorr,1) = medfilt2(all_scatter(:,1:ncorr,1),[5 5], 'symmetric');

x = [ 1:ncorr, (n2-ncorr+1):n2]; % known columns
xin = (ncorr+1):(n2-ncorr);      % unknown columns

% Rowwise scatter estimation
for i = 1:n1
    y = scattercorr(i,x);
    %pp = polyfit(x,y,2); % Quadratic fitting
    pp = polyfit(x,y,1); % linear fitting
    scattercorr(i, xin) = polyval(pp,xin);        
end

% Lateral smoothing
ufilt = 9; % The width of kernel for lateral smoothing
ufilt2 = floor(ufilt/2);

view = scattercorr;
for i = 1:n1
    if (i<=1+ufilt2)
        scattercorr(i,xin) = mean( view(1:i+ufilt2, xin) );
    elseif (i>1+ufilt2)&&(i<n1-2*ufilt2)
        scattercorr(i,xin) = mean( view(i-ufilt2:i+ufilt2, xin) );
    else
        scattercorr(i,xin) = mean( view(i-ufilt2:n1, xin) );
    end
end

scattercorr( scattercorr < 0 ) = 0.;

figure; plot(scattercorr(round(end/2),:)); title('Estimated Scatter');

figure; mesh(ixx,iyy,scattercorr); title('Estimated Scatter Data 2D'); 
xlabel('Z direction'); ylabel('X direction'); xlim([1 96]); ylim([1 512]); zlim([0 80]);

figure; plot(all_scatter(round(end/2),:,1)); title('Scatter');

figure; mesh(ixx,iyy,all_scatter(:,:,1)); title('Scatter Data 2D');
xlabel('Z direction'); ylabel('X direction'); xlim([1 96]); ylim([1 512]); zlim([0 80]);


%% Correct scatter image artifacts
% Case II: use only two rows, i.e. the (ncorr)th and (n2-ncorr+1)th row
[n1,n2,n3] = size(all_scatter);
scattercorr = zeros(n1,n2);

ncorr = 4;

% 1D median filter for the edges
scattercorr(:,1:ncorr) = medfilt1(all_scatter(:,1:ncorr,1),10);
scattercorr(:,n2-ncorr+1:n2) = medfilt1(all_scatter(:,n2-ncorr+1:n2),10);

x = [ncorr, n2-ncorr+1];    % known columns
xin = (ncorr+1):(n2-ncorr); % unknown columns

% Rowwise scatter estimation
for i = 1:n1
    y = scattercorr(i,x);
    pp = polyfit(x,y,1); % linear fitting
    scattercorr(i, xin) = polyval(pp,xin);        
end

% Lateral smoothing
ufilt = 9; % The width of kernel for lateral smoothing
ufilt2 = floor(ufilt/2);

view = scattercorr;
for i = 1:n1
    if (i<=1+ufilt2)
        scattercorr(i,xin) = mean( view(1:i+ufilt2, xin) );
    elseif (i>1+ufilt2)&&(i<n1-2*ufilt2)
        scattercorr(i,xin) = mean( view(i-ufilt2:i+ufilt2, xin) );
    else
        scattercorr(i,xin) = mean( view(i-ufilt2:n1, xin) );
    end
end

scattercorr( scattercorr < 0 ) = 0.;

figure; plot(scattercorr(round(end/2),:)); title('Estimated Scatter');

figure; mesh(ixx,iyy,scattercorr); title('Estimated Scatter Data 2D'); 
xlabel('Z direction'); ylabel('X direction'); xlim([1 96]); ylim([1 512]); zlim([0 80]);

figure; plot(all_scatter(round(end/2),:,1)); title('Scatter');

figure; mesh(ixx,iyy,all_scatter(:,:,1)); title('Scatter Data 2D');
xlabel('Z direction'); ylabel('X direction'); xlim([1 96]); ylim([1 512]); zlim([0 80]);

%% Test ScatterCorr.m

load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Pelvis5_Cone_corrected.mat');
noView = 45;
scatter1 = all_scatter(:,:,noView)';  % Input: 2D scatter data
ix = 1:96; iy=1:512;
[ixx, iyy] = meshgrid(iy,ix);
figure; mesh(ixx,iyy,scatter1); title('Scatter Data 45');
%xlabel('Z direction'); ylabel('X direction');
xlabel('X direction'); ylabel('Z direction');

scattercorr2 = ScatterCorr(scatter1);
figure; mesh(ixx,iyy,scattercorr2); title('Scatter Estimation');
%xlabel('Z direction'); ylabel('X direction');
xlabel('X direction'); ylabel('Z direction');

err1 = scatter1 - scattercorr2;
figure; mesh(ixx,iyy,err1); title('Scatter Estimation Error');
%xlabel('Z direction'); ylabel('X direction');
xlabel('X direction'); ylabel('Z direction');

view2 = 45;
figure(20); plot(scatter1(:,view2));
figure(21); plot(scattercorr2(:,view2));
figure(22); plot(err1(:,view2));


%%
% %% add path
% fprintf('\tAdding needed paths... ');
% addpath(genpath('../CTSim'));
% fprintf('\tdone.\n');



ufilt = 9; % The width of kernel for lateral smoothing
ufilt2 = floor(ufilt/2);

scattercorr = all_scatter;

[n1,n2,n3] = size(all_scatter);
x = [ 1:ncorr, (n2-ncorr+1):n2];
xin = (ncorr+1):(n2-ncorr);

tSubStep = tic;
% Preprocess the scatter results
all_scatter( all_scatter < 0 ) = 0.;
view = zeros(scattercorr);

for iv = 1:n3
    all_scatter(:,:,iv) = correctScatterImageArtifacts( squeeze(all_scatter(:,:,iv)) );
    
    % Columnwise scatter estimation
    for i = 1:n1
        y = all_scatter(i,x,iv);
        pp = polyfit(x,y,2); % Quadratic fitting
        scattercorr(i, xin, iv) = polyval(pp,xin);        
    end
    
    % Lateral smoothing
    for i = 1:n1
        if (i<=1+ufilt2)
            view(i,x,iv) = mean( squeeze(scattercorr(1:i+ufilt2, x, iv)) );
        elseif (i>1+ufilt2)&&(i<n1-2*ufilt2)
            view(i,x,iv) = mean( squeeze(scattercorr(i-ufilt2:i+ufilt2, x, iv)) );
        else
            view(i,x,iv) = mean( squeeze(scattercorr(i-ufilt2:n1, x, iv)) );
        end
    end
    
end




fprintf('(%is)\n', round(toc(tSubStep)));

%end


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

