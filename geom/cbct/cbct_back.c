/*
 *
 *pixel driven backprojection
 *
 *
 */
bool cbct_back_pd(
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
    int ib, ix, iy, iz, iu, iv, npixel, couch_z;
    int	nxy = nx * ny;
    for ( ib = 0; ib < noviews ; ib++) {
        
        float wz = - ( nz - 1.0 ) / 2.0 + offset_z;
        
        double beta = betas[ib];
        float sin_a = (float) sin( beta );
        float cos_a = (float) cos( beta );
        couch_z = couchz[ib];
        
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
                    
                    rOff = offset_u * du * sad / sdd ;
                    
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
                    
                    zc = ( iz + wz ) * dz + couch_z;
                    
                    v = zc * w2 / dv - wv;
                    
                    if( offBoundary ){
                        if ( v < 1E-6 )
                            continue;
                         if ( v > nv - 1 - 1E-6 )
                             break;
                    }else{
                        v = max( v, 0.5 );
                        v = min( v, nv - 1.5 );
                    }
                        
                    iv = floorf( v );
                    wvl = v - iv;
                    
                    npixel = ib * nu * nv + iu * nv + iv;
                    
                    image[iz*nxy+iy*nx+ix] = image[iz*nxy+iy*nx+ix]
                            + w2s * ( (1-wul) * ( (1-wvl)*proj[npixel] + wvl * proj[npixel+1] )
                            + wul * ( (1-wvl)*proj[npixel+nv] + wvl * proj[npixel+nv+1] ) );
                    
                    
                } //iz
                
            } //ix
            
        } //iy
        
    } //ib
    
    return true;
}


bool cbct_back_dd(
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
        bool *map  ){
    
    float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
    float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
    
    float wu = - ( nu - 1.0 ) / 2.0 + offset_u;
    float wv = - ( nv - 1.0 ) / 2.0 + offset_v;
    float sdd_du = sdd / du;
    
    int ib, ix, iy, iz, iu, iv, iz_start, iz_stop;
    int nxy = nx *ny;
    float u, v, uu;
    
    
    for ( ib = 0; ib < noviews ; ib++) {
        
        float wz = - ( nz - 1.0 ) / 2.0 + offset_z + ( couchz[ib] / dz );
        double beta = betas[ib];
        float sin_a = (float) sin( beta );
        float cos_a = (float) cos( beta );
        
        float dx_sin_2 = dx * sin_a / 2.0;
        float dx_cos_2 = dx * cos_a / 2.0;
        float dy_sin_2 = dy * sin_a / 2.0;
        float n_dy_cos_2 = - dy * cos_a / 2.0;
        
        float xc, yc, tproj0, dsxy0, tproj0_dsxy0, zinc, vstart;
        float cos_gamma, sin_gamma, cos_phi, sin_phi, ampu, ampv,voxel;
        float umin, umax, vmin, vmax;
        int iumin, iumax, ivmin, ivmax;
        float weights_u[FOOTPRINT_SIZE_U], weights_v[FOOTPRINT_SIZE_V];
        
        
        float *proj_view = proj + nu * nv * ib;
        // Scale the projection view with z amplitude function
        for( iu = 0; iu < nu; iu ++ ){
            
            u = ( iu +  wu ) * du;
            uu = sdd * sdd + u * u;
            
            for( iv = 0; iv < nv; iv ++){
                
                v = ( iv +  wv ) * dv;
                
                if( detector_type == 1 ){ //flat panel detector
                    ampv = sqrt( ( uu + v * v ) / uu );
                }else{ //curved detector
                    ampv = sqrt( sdd * sdd  + v * v ) / sdd;
                }
                
                proj_view[iu*nv+iv] = proj_view[iu*nv+iv] * ampv;
                
            }
        }
        
        
        for ( iy = 0; iy < ny ; iy ++ ) {
            
            for ( ix = 0; ix < nx ; ix++ ) {
                
                if( ! map[iy*nx+ix] ){
                    continue;
                }
                
                xc = ( ix + wx ) * dx;
                yc = ( iy + wy ) * dy;
                
                tproj0 = xc * cos_a - yc * sin_a ;
                dsxy0  = sad + xc * sin_a + yc * cos_a ;
                tproj0_dsxy0 = sqrt( tproj0 * tproj0 + dsxy0 * dsxy0 );
                
                cos_gamma = dsxy0 / tproj0_dsxy0;
                sin_gamma = tproj0 / tproj0_dsxy0 ;
                
                cos_phi = fabs( cos_a * cos_gamma - sin_a * sin_gamma );
                sin_phi = fabs( sin_a * cos_gamma + cos_a * sin_gamma );
                
                // compute amplitude function and footprint boundaries
                
                if( detector_type == 1 ){
                    if ( cos_phi > sin_phi) {
                        
                        ampu = dx / cos_phi / 10.0;
                        umin = sdd_du * ( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 ) - wu;
                        umax = sdd_du * ( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 ) - wu;
                        
                    } else {
                        
                        ampu = dx / sin_phi / 10.0;
                        umin = sdd_du * ( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 ) - wu;
                        umax = sdd_du * ( tproj0 - dy_sin_2 ) / ( dsxy0 - n_dy_cos_2 ) - wu;
                        
                    }
                }else{
                    if ( cos_phi > sin_phi) {
                        
                        ampu = dx / cos_phi / 10.0;
                        umin = sdd_du * atan(( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 )) - wu;
                        umax = sdd_du * atan(( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 )) - wu;
                        
                    } else {
                        
                        ampu = dx / sin_phi / 10.0;
                        umin = sdd_du * atan(( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 )) - wu;
                        umax = sdd_du * atan(( tproj0 - dy_sin_2 ) / ( dsxy0 - n_dy_cos_2 )) - wu;
                    }
                }
                if ( umin >= umax ) {
                    float temp = umax;
                    umax = umin;
                    umin = temp;
                }
                
                // decide if it outside the detector
                umin = max( umin + 0.5, 0  );
                umax = min( umax + 0.5, nu );
                
                if ( umin > umax ){
                    continue;
                }
                
                iumin = floorf( umin );
                iumax = floorf( umax );
                iumax = min( iumax, nu-1 );
                iumax = min( iumax, iumin+15 );
                
                // projection
                if (iumin == iumax) {
                    weights_u[0] =  umax - umin;
                }
                else{
                    for ( iu = iumin ; iu <= iumax ; iu++) {
                        if ( iu == iumin ) {
                            weights_u[0] = (iumin + 1.0 - umin) ;
                        } else if (iu == iumax) {
                            weights_u[iu-iumin] =  (float)( umax - iumax) ;
                        } else {
                            weights_u[iu-iumin] = 1.0;
                        }
                    }
                }
                
                zinc = dz / dv * sdd / dsxy0;
                vstart = ( wz - 0.5 ) * zinc - wv;
                
                iz_start = floorf( ( - 1 - vstart ) / zinc );
                iz_stop  = ceilf( (nv + 1 - vstart ) / zinc );
                
                iz_start = max( iz_start, 0  );
                iz_stop = min( iz_stop, nz );
                
                if ( iz_start > iz_stop ){
                    continue;
                }
                
                for ( iz = iz_start ; iz < iz_stop ; iz++ ){
                                     
                    vmin = vstart + iz * zinc;
                    vmax = vmin + zinc;
                    
                    vmin = max( vmin + 0.5, 0  );
                    vmax = min( vmax + 0.5, nv );
                    
                    if ( vmin > vmax ){
                        continue;
                    }
                    
                    ivmin = floorf( vmin );
                    ivmax = floorf( vmax );
                    ivmax = min( ivmax, nv-1 );
                    ivmax = min( ivmax, ivmin + FOOTPRINT_SIZE_V - 1 );
                    
                    //vertical footprint
                    if (ivmin == ivmax) {
                        weights_v[0] = vmax - vmin;
                    }
                    else{
                        for ( iv = ivmin ; iv <= ivmax ; iv++) {
                            if ( iv == ivmin ) {
                                weights_v[0] = ((float)ivmin + 1.0 - vmin) ;
                            } else if (iv == ivmax) {
                                weights_v[iv-ivmin] =  ( vmax - (float)ivmax) ;
                            } else {
                                weights_v[iv-ivmin] = 1.0f;
                            }
                        }
                    }
                    
                    
                    voxel = 0.0;
                    
                    for ( iu = iumin; iu <= iumax; iu++ ){
                        for ( iv = ivmin; iv <= ivmax; iv++ ){
                            voxel += proj_view[iu*nv+iv] * weights_u[iu-iumin] * weights_v[iv-ivmin];
                        }//iv
                        
                    }//iu
                    
                    image[iz*nxy+iy*nx+ix] +=  ampu * voxel;
                    
                    
                } //iz
                
                
            } //ix
            
        } //iy
        
    } //ib
    
    return true;
}


bool cbct_back_tf(
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
        bool *map  ){
    
    float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
    float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
    
    float wu = - ( nu - 1.0 ) / 2.0 + offset_u;
    float wv = - ( nv - 1.0 ) / 2.0 + offset_v;
    float sdd_du = sdd / du;
    
    int ib, ix, iy, iz, iu, iv, iz_start, iz_stop;
    int nxy = nx *ny;
    float u, v, uu;
    
    
    for ( ib = 0; ib < noviews ; ib++) {
        
        float wz = - ( nz - 1.0 ) / 2.0 + offset_z + ( couchz[ib] / dz );
        double beta = betas[ib];
        float sin_a = (float) sin( beta );
        float cos_a = (float) cos( beta );
        
        float dx_sin_2 = dx * sin_a / 2.0;
        float dx_cos_2 = dx * cos_a / 2.0;
        float dy_sin_2 = dy * sin_a / 2.0;
        float n_dy_cos_2 = - dy * cos_a / 2.0;
        
        float xc, yc, tproj0, dsxy0, tproj0_dsxy0, zinc, vstart;
        float cos_gamma, sin_gamma, cos_phi, sin_phi, ampu, ampv,voxel;
        float tau[4];
        float bl, br, xl, xr;
        float umin, umax, vmin, vmax;
        int iumin, iumax, ivmin, ivmax;
        float weights_u[FOOTPRINT_SIZE_U], weights_v[FOOTPRINT_SIZE_V];
        
        
        float *proj_view = proj + nu * nv * ib;
        // Scale the projection view with z amplitude function
        for( iu = 0; iu < nu; iu ++ ){
            
            u = ( iu +  wu ) * du;
            uu = sdd * sdd + u * u;
            
            for( iv = 0; iv < nv; iv ++){
                
                v = ( iv +  wv ) * dv;
                
                if( detector_type == 1 ){ //flat panel detector
                    ampv = sqrt( ( uu + v * v ) / uu );
                }else{ //curved detector
                    ampv = sqrt( sdd * sdd  + v * v ) / sdd;
                }
                
                proj_view[iu*nv+iv] = proj_view[iu*nv+iv] * ampv;
                
            }
        }
        
        
        for ( iy = 0; iy < ny ; iy ++ ) {
            
            for ( ix = 0; ix < nx ; ix++ ) {
                
                if( ! map[iy*nx+ix] ){
                    continue;
                }
                
                xc = ( ix + wx ) * dx;
                yc = ( iy + wy ) * dy;
                
                tproj0 = xc * cos_a - yc * sin_a ;
                dsxy0  = sad + xc * sin_a + yc * cos_a ;
                tproj0_dsxy0 = sqrt( tproj0 * tproj0 + dsxy0 * dsxy0 );
                
                cos_gamma = dsxy0 / tproj0_dsxy0;
                sin_gamma = tproj0 / tproj0_dsxy0 ;
                
                cos_phi = fabs( cos_a * cos_gamma - sin_a * sin_gamma );
                sin_phi = fabs( sin_a * cos_gamma + cos_a * sin_gamma );
                
                // compute amplitude function and footprint boundaries
                if ( cos_phi > sin_phi) {
                    ampu = dx / cos_phi / 10.0;
                } else {
                    ampu = dx / sin_phi / 10.0;
                }
                
                
                if( detector_type == 1 ){
                    
                    tau[0] = sdd_du * ( tproj0 + dx_cos_2 + dy_sin_2 ) / ( dsxy0 + dx_sin_2 + n_dy_cos_2 ) - wu;
                    tau[1] = sdd_du * ( tproj0 - dx_cos_2 + dy_sin_2 ) / ( dsxy0 - dx_sin_2 + n_dy_cos_2 ) - wu;
                    tau[2] = sdd_du * ( tproj0 + dx_cos_2 - dy_sin_2 ) / ( dsxy0 + dx_sin_2 - n_dy_cos_2 ) - wu;
                    tau[3] = sdd_du * ( tproj0 - dx_cos_2 - dy_sin_2 ) / ( dsxy0 - dx_sin_2 - n_dy_cos_2 ) - wu;
                    
                }else{
                    
                    tau[0] = sdd_du * atan(( tproj0 + dx_cos_2 + dy_sin_2 ) / ( dsxy0 + dx_sin_2 + n_dy_cos_2 )) - wu;
                    tau[1] = sdd_du * atan(( tproj0 - dx_cos_2 + dy_sin_2 ) / ( dsxy0 - dx_sin_2 + n_dy_cos_2 )) - wu;
                    tau[2] = sdd_du * atan(( tproj0 + dx_cos_2 - dy_sin_2 ) / ( dsxy0 + dx_sin_2 - n_dy_cos_2 )) - wu;
                    tau[3] = sdd_du * atan(( tproj0 - dx_cos_2 - dy_sin_2 ) / ( dsxy0 - dx_sin_2 - n_dy_cos_2 )) - wu;
                    
                }
                
                //sort taus
                if (tau[0] > tau[1]) {
                    float temp = tau[0];
                    tau[0] = tau[1];
                    tau[1] = temp;
                }
                
                
                if (tau[2] > tau[3]) {
                    float temp = tau[2];
                    tau[2] = tau[3];
                    tau[3] = temp;
                }
                
                if (tau[0] > tau[2]) {
                    float temp = tau[0];
                    tau[0] = tau[2];
                    tau[2] = temp;
                }
                
                
                if (tau[1] > tau[3]) {
                    float temp = tau[1];
                    tau[1] = tau[3];
                    tau[3] = temp;
                }
                
                if (tau[1] > tau[2]) {
                    float temp = tau[1];
                    tau[1] = tau[2];
                    tau[2] = temp;
                }
                
                // decide if it outside the detector
                umin = max( tau[0] + 0.5, 0  );
                umax = min( tau[3] + 0.5, nu );
                
                if ( umin > umax ){
                    continue;
                }
                
                iumin = floorf( umin );
                iumax = floorf( umax );
                iumax = min( iumax, nu - 1 );
                iumax = min( iumax, iumin + FOOTPRINT_SIZE_U - 1 );
                
                for ( iu = iumin ; iu <= iumax ; iu++) {
                    
                    float weight = 0;
                    bl = iu - 0.5;
                    br = iu + 0.5;
                    
                    //left piece
                    xl = max(bl, tau[0]);
                    xr = min(br, tau[1]);
                    if (xr > xl) {
                        weight = weight + (( (xr-tau[0]) * (xr-tau[0]) - (xl-tau[0]) * (xl-tau[0]) ) / (tau[1] - tau[0])) /2;
                    }
                    
                    //middle piece
                    xl = max(bl, tau[1]);
                    xr = min(br, tau[2]);
                    if (xr > xl) {
                        weight = weight + ( xr - xl );
                    }
                    
                    //right piece
                    xl = max(bl, tau[2]);
                    xr = min(br, tau[3]);
                    if (xr > xl) {
                        weight = weight + (( (xl-tau[3]) * (xl-tau[3]) - (xr-tau[3]) * (xr-tau[3]) ) / (tau[3] - tau[2])) /2;
                    }
                    
                    weights_u[iu-iumin] = weight;
                    
                } //iu
                
                zinc = dz / dv * sdd / dsxy0;
                vstart = ( wz - 0.5 ) * zinc - wv;
                
                iz_start = floorf( ( - 1 - vstart ) / zinc );
                iz_stop  = ceilf( ( nv + 1 - vstart ) / zinc );
                
                iz_start = max( iz_start, 0  );
                iz_stop = min( iz_stop, nz );
                
                if ( iz_start > iz_stop ){
                    continue;
                }
                
                for ( iz = iz_start ; iz < iz_stop ; iz++ ){
                    
                    
                    // dsz0 = ( iz + wz ) * dz ;
                    //ampv = sqrt(tproj0_dsxy0 * tproj0_dsxy0 + dsz0 * dsz0) / tproj0_dsxy0;
                    
                    vmin = vstart + iz * zinc;
                    vmax = vmin + zinc;
                    
                    vmin = max( vmin + 0.5, 0  );
                    vmax = min( vmax + 0.5, nv );
                    
                    if ( vmin > vmax ){
                        continue;
                    }
                    
                    ivmin = floorf( vmin );
                    ivmax = floorf( vmax );
                    ivmax = min( ivmax, nv-1 );
                    ivmax = min( ivmax, ivmin + FOOTPRINT_SIZE_V - 1 );
                    
                    //vertical footprint
                    if (ivmin == ivmax) {
                        weights_v[0] = vmax - vmin;
                    }
                    else{
                        for ( iv = ivmin ; iv <= ivmax ; iv++) {
                            if ( iv == ivmin ) {
                                weights_v[0] = ((float)ivmin + 1.0 - vmin) ;
                            } else if (iv == ivmax) {
                                weights_v[iv-ivmin] =  ( vmax - (float)ivmax) ;
                            } else {
                                weights_v[iv-ivmin] = 1.0f;
                            }
                        }
                    }
                    
                    
                    voxel = 0.0;
                    
                    for ( iu = iumin; iu <= iumax; iu++ ){
                        for ( iv = ivmin; iv <= ivmax; iv++ ){
                            voxel += proj_view[iu*nv+iv] * weights_u[iu-iumin] * weights_v[iv-ivmin];
                        }//iv
                        
                    }//iu
                    
                    image[iz*nxy+iy*nx+ix] +=  ampu * voxel;
                    
                    
                } //iz
                
                
            } //ix
            
        } //iy
        
    } //ib
    
    return true;
}




bool cbct_back_tf_old(
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
        bool *map  ){
    
    float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
    float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
    
    float wu = - ( nu - 1.0 ) / 2.0 + offset_u;
    float wv = - ( nv - 1.0 ) / 2.0 + offset_v;
    float sdd_du = sdd / du;
    
    int ib, ix, iy, iz, iu, iv, iz_start, iz_stop;
    int nxy = nx *ny;
    for ( ib = 0; ib < noviews ; ib++) {
        
        float wz = - ( nz - 1.0 ) / 2.0 + offset_z + ( couchz[ib] / dz );
        double beta = betas[ib];
        float sin_a = (float) sin( beta );
        float cos_a = (float) cos( beta );
        
        float dx_sin_2 = dx * sin_a / 2.0;
        float dx_cos_2 = dx * cos_a / 2.0;
        float dy_sin_2 = dy * sin_a / 2.0;
        float n_dy_cos_2 = - dy * cos_a / 2.0;
        
        float xc, yc, tproj0, dsxy0, tproj0_dsxy0, dsz0, zinc, vstart;
        float cos_gamma, sin_gamma, cos_phi, sin_phi, ampu, ampv,voxel;
        float tau[4];
        float bl, br, xl, xr;
        float umin, umax, vmin, vmax;
        int iumin, iumax, ivmin, ivmax;
        float weights_u[16], weights_v[16];
        
        
        for ( iy = 0; iy < ny ; iy ++ ) {
            
            for ( ix = 0; ix < nx ; ix++ ) {
                
                if( ! map[iy*nx+ix] ){
                    continue;
                }
                
                xc = ( ix + wx ) * dx;
                yc = ( iy + wy ) * dy;
                
                tproj0 = xc * cos_a - yc * sin_a ;
                dsxy0  = sad + xc * sin_a + yc * cos_a ;
                tproj0_dsxy0 = sqrt( tproj0 * tproj0 + dsxy0 * dsxy0 );
                
                cos_gamma = dsxy0 / tproj0_dsxy0;
                sin_gamma = tproj0 / tproj0_dsxy0 ;
                
                cos_phi = fabs( cos_a * cos_gamma - sin_a * sin_gamma );
                sin_phi = fabs( sin_a * cos_gamma + cos_a * sin_gamma );
                
                // compute amplitude function and footprint boundaries
                if ( cos_phi > sin_phi) {
                    ampu = dx / cos_phi / 10.0;
                } else {
                    ampu = dx / sin_phi / 10.0;
                }
                
                
                if( detector_type == 1 ){
                    
                    tau[0] = sdd_du * ( tproj0 + dx_cos_2 + dy_sin_2 ) / ( dsxy0 + dx_sin_2 + n_dy_cos_2 ) - wu;
                    tau[1] = sdd_du * ( tproj0 - dx_cos_2 + dy_sin_2 ) / ( dsxy0 - dx_sin_2 + n_dy_cos_2 ) - wu;
                    tau[2] = sdd_du * ( tproj0 + dx_cos_2 - dy_sin_2 ) / ( dsxy0 + dx_sin_2 - n_dy_cos_2 ) - wu;
                    tau[3] = sdd_du * ( tproj0 - dx_cos_2 - dy_sin_2 ) / ( dsxy0 - dx_sin_2 - n_dy_cos_2 ) - wu;
                    
                }else{
                    
                    tau[0] = sdd_du * atan(( tproj0 + dx_cos_2 + dy_sin_2 ) / ( dsxy0 + dx_sin_2 + n_dy_cos_2 )) - wu;
                    tau[1] = sdd_du * atan(( tproj0 - dx_cos_2 + dy_sin_2 ) / ( dsxy0 - dx_sin_2 + n_dy_cos_2 )) - wu;
                    tau[2] = sdd_du * atan(( tproj0 + dx_cos_2 - dy_sin_2 ) / ( dsxy0 + dx_sin_2 - n_dy_cos_2 )) - wu;
                    tau[3] = sdd_du * atan(( tproj0 - dx_cos_2 - dy_sin_2 ) / ( dsxy0 - dx_sin_2 - n_dy_cos_2 )) - wu;
                    
                }
                
                //sort taus
                if (tau[0] > tau[1]) {
                    float temp = tau[0];
                    tau[0] = tau[1];
                    tau[1] = temp;
                }
                
                
                if (tau[2] > tau[3]) {
                    float temp = tau[2];
                    tau[2] = tau[3];
                    tau[3] = temp;
                }
                
                if (tau[0] > tau[2]) {
                    float temp = tau[0];
                    tau[0] = tau[2];
                    tau[2] = temp;
                }
                
                
                if (tau[1] > tau[3]) {
                    float temp = tau[1];
                    tau[1] = tau[3];
                    tau[3] = temp;
                }
                
                if (tau[1] > tau[2]) {
                    float temp = tau[1];
                    tau[1] = tau[2];
                    tau[2] = temp;
                }
                
                // decide if it outside the detector
                umin = max( tau[0] + 0.5, 0  );
                umax = min( tau[3] + 0.5, nu );
                
                if ( umin > umax ){
                    continue;
                }
                
                iumin = floorf( umin );
                iumax = floorf( umax );
                iumax = min( iumax, nu - 1 );
                iumax = min( iumax, iumin + 15 );
                
                for ( iu = iumin ; iu <= iumax ; iu++) {
                    
                    float weight = 0;
                    bl = iu - 0.5;
                    br = iu + 0.5;
                    
                    //left piece
                    xl = max(bl, tau[0]);
                    xr = min(br, tau[1]);
                    if (xr > xl) {
                        weight = weight + (( (xr-tau[0]) * (xr-tau[0]) - (xl-tau[0]) * (xl-tau[0]) ) / (tau[1] - tau[0])) /2;
                    }
                    
                    //middle piece
                    xl = max(bl, tau[1]);
                    xr = min(br, tau[2]);
                    if (xr > xl) {
                        weight = weight + ( xr - xl );
                    }
                    
                    //right piece
                    xl = max(bl, tau[2]);
                    xr = min(br, tau[3]);
                    if (xr > xl) {
                        weight = weight + (( (xl-tau[3]) * (xl-tau[3]) - (xr-tau[3]) * (xr-tau[3]) ) / (tau[3] - tau[2])) /2;
                    }
                    
                    weights_u[iu-iumin] = weight;
                    
                } //iu
                
                zinc = dz / dv * sdd / dsxy0;
                vstart = ( wz - 0.5 ) * zinc - wv;
                
                iz_start = floorf( ( - 1 - vstart ) / zinc );
                iz_stop  = ceilf( ( nv + 1 - vstart ) / zinc );
                
                iz_start = max( iz_start, 0  );
                iz_stop = min( iz_stop, nz );
                
                if ( iz_start > iz_stop ){
                    continue;
                }
                
                for ( iz = iz_start ; iz < iz_stop ; iz++ ){
                    
                    
                    dsz0 = ( iz + wz ) * dz ;
                    
                    ampv = sqrt(tproj0_dsxy0 * tproj0_dsxy0 + dsz0 * dsz0) / tproj0_dsxy0;
                    
                    vmin = vstart + iz * zinc;
                    vmax = vmin + zinc;
                    
                    vmin = max( vmin + 0.5, 0  );
                    vmax = min( vmax + 0.5, nv );
                    
                    if ( vmin > vmax ){
                        continue;
                    }
                    
                    ivmin = floorf( vmin );
                    ivmax = floorf( vmax );
                    ivmax = min( ivmax, nv-1 );
                    ivmax = min( ivmax, ivmin+15 );
                    
                    //vertical footprint
                    if (ivmin == ivmax) {
                        weights_v[0] = vmax - vmin;
                    }
                    else{
                        for ( iv = ivmin ; iv <= ivmax ; iv++) {
                            if ( iv == ivmin ) {
                                weights_v[0] = (ivmin + 1.0 - vmin) ;
                            } else if (iv == ivmax) {
                                weights_v[iv-ivmin] =  (float)( vmax - ivmax) ;
                            } else {
                                weights_v[iv-ivmin] = 1.0;
                            }
                        }
                    }
                    
                    
                    voxel = 0.0;
                    
                    for ( iu = iumin; iu <= iumax; iu++ ){
                        for ( iv = ivmin; iv <= ivmax; iv++ ){
                            voxel = voxel + proj[ib*nu*nv+iu*nv+iv] * weights_u[iu-iumin] * weights_v[iv-ivmin];
                        }//iv
                        
                    }//iu
                    
                    image[iz*nxy+iy*nx+ix] = image[iz*nxy+iy*nx+ix]  +  ampu * ampv * voxel;
                    
                    
                } //iz
                
                
            } //ix
            
        } //iy
        
    } //ib
    
    return true;
}




