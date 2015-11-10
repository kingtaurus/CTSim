function sinoOut = bilateralFilterProjectionDomain( sinoAtt, geom, w, sd, ss )

fprintf('Bilateral filter in projection domain with sd = %2.2f, ss = %2.2f ... \n', sd, ss);

sinoOut = zeros( size(sinoAtt), 'single');

if length( geom.detSize ) == 2
    
    for iv = 1 : geom.noViews
        sinoOut(:,:,iv) = BilateralFilter( squeeze( sinoAtt(:,:,iv)), w, sd, ss);
    end
else
    
    for iv = 1 : geom.noViews
        sinoOut(:,iv) = BilateralFilter( squeeze( sinoAtt(:,:,iv)), w, sd, ss );
    end
end

fprintf('\t done.\n')