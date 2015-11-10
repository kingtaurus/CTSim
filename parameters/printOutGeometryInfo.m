function printOutGeometryInfo( geom )

fprintf('Geometry information:\n');
fprintf('\tSDD: %.0fmm, SAD: %.0fmm\n', geom.SDD, geom.SAD);

if geom.helicalScan
    fprintf('\tProjections per Turn: %i, number of turns: %i with pitch %f\n', geom.noViewsTurn, geom.noTurns, geom.pitch);
else
    fprintf('\tNumber of projections: %i \n', geom.noViews);
end

fprintf('\tDetector size: %ipx X %ipx\n',  geom.detSize(1), geom.detSize(end));
fprintf('\tPixel size: %5.3fmm X %2.3f mm  \n',  geom.detSpacing(1), geom.detSpacing(end) );
fprintf('\tReconstruction size: %i  X %i  X %i  \n',  geom.reconSize(1),  geom.reconSize(2), geom.reconSize(end) );
fprintf('\tReconstruction spacing: %2.3fmm X %2.3f mm X %2.3f mm  \n',  geom.reconSpacing(1), geom.reconSpacing(2), geom.reconSpacing(end) );
