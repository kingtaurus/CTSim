airIndexX = validIndexX;
airIndexY = [37 72];
plot( airIndexY, airIndexX, 'r' );
plot( airIndexY(2:-1:1) , airIndexX, 'r' );
plot( airIndexY([1 1]), airIndexX, 'r' );
plot( airIndexY([2 2]), airIndexX, 'r' );
plot( airIndexY, airIndexX([1 1]), 'r' );
plot( airIndexY, airIndexX([2 2]), 'r' );


if calibrateLag
    
    regionIdx1 = [180 210 290 460];
    regionIdx2 = [185 215 055 220];
    regionIdx3 = [180 220 240 275];
    
    regionCount1 = zeros( noViews, 1 );
    regionCount2 = zeros( noViews, 1 );
    regionCount3 = zeros( noViews, 1 );
    
    for iv = firstValidViewNum : lastValidViewNum
        
        view = squeeze( sinoPhotonCount( :, :, iv) );
        
        region = view( regionIdx1(1):regionIdx1(2), regionIdx1(3):regionIdx1(4));
        regionCount1( iv ) = mean( region(:) );
        
        region = view( regionIdx2(1):regionIdx2(2), regionIdx2(3):regionIdx2(4));
        regionCount2( iv ) = mean( region(:) );
        
        region = view( regionIdx3(1):regionIdx3(2), regionIdx3(3):regionIdx3(4));
        regionCount3( iv ) = mean( region(:) );
        
    end
    
    
    ave1 = mean( regionCount1( 10:40 ) );
    ave2 = mean( regionCount2( 10:40 ) );
    ave3 = mean( regionCount3( 10:40 ) );
    
    start = 50;
    stop = 60;
    
    
    figure;
    plot( regionCount1 ); hold on;
    plot( regionCount2, '-.' ); hold on;
    plot( regionCount3, ':' ); hold on;
    title( ['Detector step response 6 MV']);
    xlabel 'No. of frames';
    legend( num2str(ave1), num2str(ave2), num2str(ave3) );
    
    
    figure;
    plot( regionCount1/ave1 ); hold on;
    plot( regionCount2/ave2, '-.' ); hold on;
    plot( regionCount3/ave3, ':' ); hold on;
    title( ['Detector step response (Normalized) 6 MV']);
    xlabel 'No. of frames';
    legend( num2str(ave1), num2str(ave2), num2str(ave3) );
    
    figure;
    semilogy( regionCount1(start:stop)/ave1 ); hold on;
    semilogy( regionCount2(start:stop)/ave2, '-.' ); hold on;
    semilogy( regionCount3(start:stop)/ave3, ':' ); hold on;
    title( ['Falling edge (Normalized) 6 MV']);
    xlabel 'No. of frames';
    legend( num2str(ave1), num2str(ave2), num2str(ave3) );
    
    
end

%% linear lag fitting

if 0
    
    
    aveCounts = mean( refCounts(firstValidViewNum:lastValidViewNum));
    
    figure;
    plot( refCounts / aveCounts  );
    title( ['Detector step response ( mean =' num2str(aveCounts) ' )']);
    xlabel 'No. of frames';
    
    
    figure;
    plot( refCounts( lastValidViewNum:end)  /  aveCounts );
    title( ['Detector falling edge ( mean =' num2str( aveCounts ) ' )']);
    xlabel 'No. of frames';
    
    
    figure;
    semilogy( refCounts( lastValidViewNum:end)  /  aveCounts );
    title( ['Detector falling edge ( mean =' num2str( aveCounts ) ' )']);
    xlabel 'No. of frames';
    
    
    
    y = refCounts( lastValidViewNum+1:end) / refCounts( lastValidViewNum );
    [b, a, f] = fallingLagLeastSquaresFittingGaussNewton( y,  y(1), 1 );
    
    sinoPhotonCount = linearRecursiveLagCorrection( sinoPhotonCount, b, a );
    
    
    % sinoPhotonCount(:,:,2:end) = sinoPhotonCount(:,:,2:end) - 0.1 * sinoPhotonCount(:,:,1:end-1);
    
end
