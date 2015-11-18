load 'temp.mat'

[phan, map] = loadXCATPhantom(p);

%% first of all generate prime simulation data, here is the helical case

turns = 3;
pitch = 1/2;

[ geom ] = loadProjectionGeometryHelicalCT( p, turns, pitch );

spectrum = loadSpectraCT(p, geom, 1e6);

sinosHelicalDir = [sinosDirKeV sprintf( '%i-turns-%i-pitch', turns, round( pitch * 100 ) ), '/'];

sinoRaw = simulateCTRawData( phan, geom, spectrum, sinosHelicalDir  );

%% circular scan case

[ geom ] = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 1e6);

sinoRaw = simulateCTRawData( phan, geom, spectrum, sinosDirKeV  );



%% example of reading erica scatter simulation file

% create space to save image
dims = [96 512 25];
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

load( 'E:\MATLAB\CTData\ScatterSimulation\scatter_simulation_xcat_5\scatter_simulations_5.mat' );

sinoPrime = rotateSinogram( primary_bowtie, 0, 1 );
sinoPrime = reverseSingoram( sinoPrime );

sinoScatter = rotateSinogram( scatter_bowtie, 0, 1 );
sinoScatter = reverseSingoram( sinoScatter );

p = mean(sinoRaw(sinoPrime(:) > 2 )) / mean( sinoPrime(sinoPrime(:) > 2 ));

sinoPrime         = p * sinoPrime;
sinoScatter       = p * sinoScatter;

figure; imdisp( (log10( sinoRaw )), [0 6]);
figure; imdisp( (log10( sinoPrime )), [0 6]);
figure; imdisp( sinoScatter,[0 10000] );