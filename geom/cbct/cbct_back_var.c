bool cbct_back_var(
        int nx,
        int ny,
        int nz,
        float dx,
        float dy,
        float dz,
        float offset_x,
        float offset_y,
        float offset_z,
        int nu,
        int nv,
        float du,
        float dv,
        float offset_u,
        float offset_v,
        float sad,
        float add,
        float sdd,
        int noviews,
        double *betas,
        float* couchz,
        int detector_type,
        float *image,
        float *proj,
        bool *map,
        bool offBoundary ){
    
    float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
    float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
    
    float wu = - ( nu - 1.0 ) / 2.0 + offset_u;
    float wv = - ( nv - 1.0 ) / 2.0 + offset_v;
    
    float xc, yc, zc, xRot, yRot, rOff;
    float w2, w2s, u, v, wul, wvl;
    int ib, ix, iy, iz, iu, iv, npixel;
    int	nxy = nx * ny;
    for ( ib = 0; ib < noviews ; ib++) {
        
        float wz = - ( nz - 1.0 ) / 2.0 + offset_z + ( couchz[ib] / dz );
        
        double beta = betas[ib];
        float sin_a = (float) sin( beta );
        float cos_a = (float) cos( beta );
        
        
        for ( iy = 0; iy < ny ; iy ++ ) {
            
            for ( ix = 0; ix < nx ; ix++ ) {
                
                if( ! map[iy*nx+ix] ){
                    continue;
                }
                
                xc = ( ix + wx ) * dx;
                yc = ( iy + wy ) * dy;
                
                xRot = xc * cos_a - yc * sin_a ;
                yRot = xc * sin_a + yc * cos_a ;
                
                if( detector_type == 1 ){
                    rOff = 0;
                    w2 = sdd / ( sad + yRot );
                    w2s = w2 * w2;
                    
                    u = xRot * w2 / du - wu;
                }else{
                    
                    rOff = offset_u * sad / sdd ;
                    
                    u = atan( ( xRot - rOff  ) / ( sad + yRot  ) ) * sdd / du - wu;
                    
                    w2s = sdd * sdd  / ( ( sad + yRot )*( sad + yRot ) + (xRot - rOff) * (xRot - rOff) );
                    w2 = sdd / ( sad + yRot );
                    
                }
                
                if ( u < 1E-6 || u > nu - 1 - 1E-6 ) {
                    continue;
                }
                
                
                iu = floorf( u );
                wul = u - iu;
                
                for ( iz = 0; iz < nz ; iz ++ ){
                    
                    zc = ( iz + wz ) * dz;
                    
                    v = zc * w2 / dv - wv;
                    
                    if ( v < 1E-6 ) {
                        if ( offBoundary ){
                            continue;
                        } else {
                            v = 0.5;
                        }
                    }
                    
                    
                    if ( v > nv - 1 - 1E-6 ) {
                        if( offBoundary ){
                            break;
                        } else {
                            v = nv - 1.5;
                        }
                    }
                    
                    iv = floorf( v );
                    wvl = v - iv;
                    
                    npixel = ib * nu * nv + iu * nv + iv;
                    
                    image[iz*nxy+iy*nx+ix] = image[iz*nxy+iy*nx+ix]
                            + w2s * w2s * ( (1-wul) * ( (1-wvl)*proj[npixel] + wvl * proj[npixel+1] )
                            + wul * ( (1-wvl)*proj[npixel+nv] + wvl * proj[npixel+nv+1] ) );
                    
                    
                } //iz
                
            } //ix
            
        } //iy
        
    } //ib
    
    return true;
}

