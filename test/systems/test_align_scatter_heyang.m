load 'temp.mat'

[phan, map] = loadXCATPhantom(p);
%[phan, map] = loadXCATPhantomDistorted(p);

%% first of all generate prime simulation data, here is the helical case

% turns = 3;
% pitch = 1/2;

turns = 2;
%pitch = 0.9;
pitch = 0.3;


[ geom ] = loadProjectionGeometryHelicalCT( p, turns, pitch );

spectrum = loadSpectraCT(p, geom, 1e6);

sinosHelicalDir = [sinosDirKeV sprintf( '%i-turns-%i-pitch', turns, round( pitch * 100 ) ), '/'];

sinoRaw = simulateCTRawData( phan, geom, spectrum, sinosHelicalDir  );

%% circular scan case

phan.offset(1) = 0;
phan.offset(2) = 0;

[ geom ] = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e6);

sinoRaw = simulateCTRawData( phan, geom, spectrum, sinosDirKeV  );



%% example of reading erica scatter simulation file

% create space to save image
%dims = [96 512 25];
dims = [512 96 8];
sinoPrime = zeros( dims, 'single' );
sinoAll = zeros( dims, 'single' );
sinoGrid = zeros( dims, 'single' );


for i = 1 : dims(3)
    
    filename = sprintf( 'Scout_helical.dat_00%02i', i - 1 );
    disp( filename );

    [curved_primary, curved_all, curved_grid] = readScatterSimulationResults( filename, dims );
    
    k = i;
    % the direcation may be reversed
    % k = dims(3) - i + 1;
    
    
    sinoPrime(:,:,k) = single( (   curved_primary' ) );
    sinoAll(:,:,k) = single( (   curved_all' ) );
    sinoGrid(:,:,k) = single( ( curved_grid' ) );
    
end

clear curved_primary curved_all curved_grid;

%% visual comparison

%close all;

sinoPrime_center = rotateSinogram( sinoPrime( 33:64, :, : ), 0, 0, 0, 0);

sinoPrime_center = reverseSingoram( sinoPrime_center );

figure( 'Name', 'sinoPrime_center' ); imdisp( log( sinoPrime_center + 0.001 ));

figure( 'Name', 'sinoRaw' ); imdisp( log( sinoRaw ) );

% figure( 'Name', 'scatters' ); imdisp( sinoGrid - sinoPrime, [0 1] );
% figure( 'Name', 'prime' ); imdisp( sinoPrime, [0 10] );
%%

close all;


p = mean(sinoRaw(sinoPrime(:) > 2 )) / mean( sinoPrime(sinoPrime(:) > 2 ));

sinoScaledPrime         = p * sinoPrime;
sinoScaledScatter       = p * sinoScatter;
sinoScaledScatterGrid   = p * sinoScatterGrid;


figure; imdisp( log( sinoRaw ));
figure; imdisp( log( sinoScaledPrime ));
figure; imdisp( ( sinoScaledScatter ), [0 10000]  );
figure; imdisp( ( sinoScaledScatterGrid ) );
figure; imdisp(  abs( sinoScaledPrime -  sinoRaw ) ./ sinoRaw, [0 1] );

%% This is an example of testing the alignment with erica's processed scatter data 

% load( 'E:\MATLAB\CTData\ScatterSimulation\scatter_simulation_xcat_4\scatter_simulations_4.mat' );

%load( 'E:\MATLAB\CTData\ScatterSimulation\scatter_simulation_xcat_5\scatter_simulations_5.mat' );

%load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Pelvis_Cone_Scout.mat');
%load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Helical_Scout.mat');
%load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\DistortedPelvisData.mat');
load( 'E:\MATLAB\CTData\scatter_simulation_pelvis\Pelvis5_Scout.mat');



%sinoPrime = rotateSinogram( primary_bowtie, 0, 1 );
%sinoPrime = rotateSinogram( primary, 0, 1, 0, 0);
%sinoPrime = rotateSinogram( Pelvis1_primary, 0, 1);
sinoPrime = rotateSinogram( all_primary, 0, 1);
sinoPrime = reverseSingoram( sinoPrime );

%sinoScatter = rotateSinogram( scatter_bowtie, 0, 1 );
%sinoScatter = rotateSinogram( Pelvis1_scatter, 0, 1);
sinoScatter = rotateSinogram( all_scatter, 0, 1);
sinoScatter = reverseSingoram( sinoScatter );

p1 = mean(sinoRaw(sinoPrime(:) > 2 )) / mean( sinoPrime(sinoPrime(:) > 2 ));
%p1 = mean(sinoRaw(sinoPrime(35:60,:,:) > 2))/mean( sinoPrime(sinoPrime(35:60,:,:) > 2));

sinoPrime         = p1 * sinoPrime;
sinoScatter       = p1 * sinoScatter;
%
figure; imdisp( (log10( sinoRaw )), [0 8]);
figure; imdisp( (log10( sinoPrime )), [0 8]);
% %figure; imdisp( sinoScatter,[0 10000] );

% p1 = mean(sinoRaw(sinoPrime(35:60,:,:) > 2))/mean( sinoPrime(sinoPrime(35:60,:,:) > 2));
% sinoPrime = p1*sinoPrime;
% figure; imdisp( (log10( sinoRaw(35:60,:,:) )), [0 8]);
% figure; imdisp( (log10( sinoPrime(35:60,:,:) )), [0 8]);