function downsampledView = detectorPixelBinning( view, binSize )
% function downsampleView = detectorPixelBinning( view, binSize )
% 
% Meng Wu

if binSize(1) == 1 && binSize(2) == 1
    downsampledView = view;
    return;
end

[ny, nx] = size( view );
downsampledView = zeros( [ny/binSize(1) nx/binSize(2)], class( view) );

for i = 1:binSize(1)
    for j = 1:binSize(2)
        
        downsampledView = downsampledView + view(i:binSize(1):end, j:binSize(2):end );

    end
end


end