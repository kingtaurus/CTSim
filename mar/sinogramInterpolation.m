function sinoInterp = sinogramInterpolation( sino, sinoMapUnknown, interpDim, smoothMargin, sinoMapKnown )

if nargin < 3
    interpDim = 1;
end

if nargin < 4
    smoothMargin = 10;
end

if nargin < 5
    sinoMapKnown = imdilate3( sinoMapUnknown, ones(5) ) & ~ sinoMapUnknown;
end

h = fspecial( 'Gaussian', [1 8], 2);

sinoInterp = sino;

if ismatrix(sino)
    
    % sinoMapKnown  =  ~ sinoMapUnknown;
    % interpolate missig pixels using linear interpolation
    [indXSinoMapKnown, indYSinoMapKnown] = find(sinoMapKnown);
    [indXSinoMapUnknown, indYSinoMapUnknown] = find(sinoMapUnknown);
    
    sinoUnknown = griddata(indXSinoMapKnown, indYSinoMapKnown, ...
        double( sino(sinoMapKnown) ), indXSinoMapUnknown, indYSinoMapUnknown, 'linear');
    
    sinoInterp(sinoMapUnknown(:)) = sinoUnknown ;
    
elseif ndims(sino) == 3
    
    for iv = 1:size( sino, 3)
        
        sinoView            = squeeze(sino(:,:,iv));
        sinoMapUnknownView  = squeeze(sinoMapUnknown(:,:,iv) );
        sinoMapKnownView    =  squeeze(sinoMapKnown(:,:,iv) );
        
        if interpDim == 1 % 1D interpolation
            
            for iz = 1 : size( sino, 1)
                
                slice       = squeeze( sinoView(iz, : ) );
                slcKnown    = squeeze( sinoMapKnownView(iz, : ) );
                slcUnknown  = squeeze( sinoMapUnknownView(iz, : ) ) ;
                
                if sum( slcUnknown ) == 0
                    continue;
                end
                
                indxKnown = find(  slcKnown );
                indxUnknown = find( slcUnknown );
                
                slice(slcUnknown) = interp1( indxKnown, slice(slcKnown), indxUnknown, 'pchip' );
                
                sinoView(iz,:) = slice';
                
            end
            
        elseif interpDim == 2 % 2D interpolation using Matlab griddata function
            
            % sinoView            = squeeze(sino(:,:,iv));
            % sinoMapUnknownView  = squeeze(sinoMapUnknown(:,:,iv) );
            % sinoMapKnownView    =  squeeze(sinoMapKnown(:,:,iv) );
            
            if sum( sinoMapUnknownView(:) ) == 0
                continue;
            end
            
            [indXSinoMapKnown, indYSinoMapKnown] = find(sinoMapKnownView);
            [indXSinoMapUnknown, indYSinoMapUnknown] = find(sinoMapUnknownView);
            
            % interpolation methods
            sinoUnknown = griddata(indXSinoMapKnown, indYSinoMapKnown, ...
                double( sinoView(sinoMapKnownView) ), indXSinoMapUnknown, indYSinoMapUnknown, 'linear');
            
            sinoView(sinoMapUnknownView(:)) = sinoUnknown ;
            
            
        end
        
        sinoViewFilt = sinoView;
        
        if smoothMargin > 0;
            % find the boundary pixel
            sinoMapBoudaries = imdilate( sinoMapUnknownView, ones(1, smoothMargin ) ) & ~ sinoMapUnknownView;
            
            % smooth the interpolated sinogram
            sinoViewFilt = imfilter( sinoViewFilt, h );
            
            % find the smooth parts in the boundary region
            sinoMapSmooth = abs( sinoView - sinoViewFilt ) < 0.1 & sinoMapBoudaries;
            sinoMapSmooth = imerode( sinoMapSmooth, ones(2) );
            
            % use the original values in these smooth pars
            sinoViewFilt( sinoMapSmooth ) = sinoView( sinoMapSmooth );
            
            % use the original values for the rest of parts
            sinoViewFilt( ~sinoMapUnknownView & ~sinoMapBoudaries  ) = sinoView( ~sinoMapUnknownView & ~sinoMapBoudaries );
            
        end
        
        sinoInterp(:,:,iv) = sinoViewFilt;
    end
end
