load 'temp.mat';
displayResult = true;
close all;


% Load simulation parameters and datas
[phan, map, roi] = loadMaterialsPhantom(p);

[geom ] = loadProjectionGeometryCarm( p, 200, 0 );

spectrum = loadSpectra(p, geom);


%%
MaxConcentration = 0.0013;
noViews = geom.noViews;

TempContrast0 = MaxConcentration * ones(1, noViews);


TempContrast1 = (0 : noViews-1) * MaxConcentration / (noViews-1);

TempContrast2 = zeros(1, noViews );

halfViews = round(noViews/2);

TempContrast2(end-halfViews+1:end) = ( 1 : halfViews) * MaxConcentration / halfViews;

%%

sinoPhotonCounts0 = computePhotonCountDynamicSinogram( phan, geom, ...
    spectrum, sinosDirKeV, TempContrast0 /2 );


sinoPhotonCounts1 = computePhotonCountDynamicSinogram( phan, geom, ...
    spectrum, sinosDirKeV, TempContrast1 );

sinoPhotonCounts2 = computePhotonCountDynamicSinogram( phan, geom, ...
    spectrum, sinosDirKeV, TempContrast2 );

sinoAtt0 = computeSinogramAttunation( sinoPhotonCounts0, spectrum );
sinoAtt1 = computeSinogramAttunation( sinoPhotonCounts1, spectrum );
sinoAtt2 = computeSinogramAttunation( sinoPhotonCounts2, spectrum );


%% FBP reconstrunction
imgAtt0 = reconFBPShortScan( sinoAtt0, geom);
imgHu0 = convertMonoAttToHu( imgAtt0, spectrum);


imgAtt1 = reconFBPShortScan( sinoAtt1, geom);
imgHu1 = convertMonoAttToHu( imgAtt1, spectrum);

imgAtt2 = reconFBPShortScan( sinoAtt2, geom);
imgHu2 = convertMonoAttToHu( imgAtt2, spectrum);

%%

if displayResult 
    figure('Name', 'FBP KeV Image', 'NumberTitle', 'off');
    imshow( imgHu1, [-500 500]);
    saveas(gcf,  'FBP1_3.eps');
end

if displayResult 
    figure('Name', 'FBP KeV Image', 'NumberTitle', 'off');
    imshow( imgHu2, [-500 500]);
    saveas(gcf,  'FBP2_3.eps');
end

%%

if displayResult 
    figure('Name', 'FBP KeV Image', 'NumberTitle', 'off');
    imshow( imgHu1-imgHu0, [-100 100]);
    saveas(gcf,  'FBP10_3.eps');
end


if displayResult 
    figure('Name', 'FBP KeV Image', 'NumberTitle', 'off');
    imshow( imgHu2-imgHu0, [-200 100]);
    saveas(gcf,  'FBP20_3.eps');
end

%%
figure;
plot( TempContrast1 );
axis([1 noViews 0 MaxConcentration]);
xlabel( 'Projection View', 'fontSize', 16);
ylabel( 'Contrast Concentration', 'fontSize', 16);
saveas(gcf,  'Contrast1.eps');