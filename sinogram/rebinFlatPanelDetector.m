function sinoOut = rebinFlatPanelDetector( sinoIn, geomIn, geomOut )

if length( geomIn.detSize ) == 1 && length( geomOut.detSize ) == 1
    
    sinoOut = zeros( [ geomOut.detSize(1), geomIn.noViews ], 'single' );
    
    nu = geomIn.detSize(1);
    du = geomIn.detSpacing(1);
    ou = geomIn.detOffset(1);
    
    u = ( ( -(nu-1)/2:(nu-1)/2)+ ou ) * du ;
    
    nu = geomOut.detSize(1);
    du = geomOut.detSpacing(1);
    ou = geomOut.detOffset(1);
    
    ui = ( ( -(nu-1)/2:(nu-1)/2)+ ou ) * du ;
    
    sinoOut = interp1(u, sinoIn, ui);
    
    sinoOut( isnan(sinoOut) ) = 0;
    
elseif  length( geomIn.detSize ) == 2 && length( geomOut.detSize ) == 2
    
    sinoOut = zeros( [geomOut.detSize(2), geomOut.detSize(1), geomIn.noViews ], 'single' );
    
    nu = geomIn.detSize(1);
    nv = geomIn.detSize(2);
    du = geomIn.detSpacing(1);
    dv = geomIn.detSpacing(2);
    ou = geomIn.detOffset(1);
    ov = geomIn.detOffset(2);
    
    u = ( ( -(nu-1)/2:(nu-1)/2)+ ou ) * du ;
    v = ( ( -(nv-1)/2:(nv-1)/2)+ ov ) * dv ;
    
    [uu, vv ] = meshgrid( u, v);
    
    nu = geomOut.detSize(1);
    nv = geomOut.detSize(2);
    du = geomOut.detSpacing(1);
    dv = geomOut.detSpacing(2);
    ou = geomOut.detOffset(1);
    ov = geomOut.detOffset(2);
    
    ui = ( ( -(nu-1)/2:(nu-1)/2)+ ou ) * du ;
    vi = ( ( -(nv-1)/2:(nv-1)/2)+ ov ) * dv ;
    
    
    [uui, vvi ] = meshgrid( ui, vi);
    
    for iv = 1:geomIn.noViews
        
        sinoOutView = interp2(uu, vv, squeeze(sinoIn(:,:,iv)), uui, vvi);
        
        sinoOutView( isnan(sinoOutView) ) = 0;
        
        sinoOut(:,:,iv) = sinoOutView;
        
    end
    
end




end