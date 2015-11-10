#include "ct_def.h";
/*
 *
 *pixel driven backprojection
 *
 *
 */
bool cbct_proj_dd(
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
        int noviews,
        float * pmat,
        float * cameras,
        float * sdds,
        float * offset_uv,
        float * image,
        float * proj,
        bool * map,
        bool offBoundary){
    
    float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
    float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
    float wz = - ( nz - 1.0 ) / 2.0 + offset_z;
    float wu = - ( nu - 1.0 ) / 2.0 ;
    float wv = - ( nv - 1.0 ) / 2.0 ;
    
    float *pm;
    int	nxy = nx * ny;
    
    float u, v, uu;
    
    float xc, yc, sdd, camera_x, camera_y, dsp_x, dsp_y;
    float u0, v0, f0, dudz, dvdz, umin, umax, vmin, vmax, ampu, ampv, voxel;
    int iumin, iumax, ivmin, ivmax;
    float weights_u[FOOTPRINT_SIZE_U], weights_v[FOOTPRINT_SIZE_V];
    
    int ib, ix, iy, iz, iu, iv, npixel;
    
    for ( ib = 0; ib < noviews ; ib++) {
        
        float *proj_view = proj + nu * nv * ib;
        
        pm = pmat + 12 * ib;
        
        // source-to-detector distance
        sdd = sdds[ ib ];
        
        // camera position
        camera_x = cameras[ 3 * ib ];
        camera_y = cameras[ 3 * ib + 1 ];
        
        for ( iy = 0; iy < ny ; iy ++ ) {
            
            for ( ix = 0; ix < nx ; ix++ ) {
                
                if( ! map[iy*nx+ix] ){
                    continue;
                }
                
                xc = ( ix + wx ) * dx;
                yc = ( iy + wy ) * dy;
                
                u0 = xc * pm[0] + yc * pm[3] + wz * dz * pm[6] + pm[9];
                v0 = xc * pm[1] + yc * pm[4] + wz * dz * pm[7] + pm[10];
                f0 = xc * pm[2] + yc * pm[5]  + pm[11];
                
                // source to voxel center distance
                dsp_x = fabs( camera_x - xc );
                dsp_y = fabs( camera_y - yc );
                
                
                // two cases to draw rectangular footprints
                if( fabs(dsp_x) > fabs( dsp_y ) ){
                    umin = ( u0 - dy * pm[3] / 2 ) / ( f0 - dy * pm[5] / 2 ) ;
                    umax = ( u0 + dy * pm[3] / 2 ) / ( f0 + dy * pm[5] / 2 );
                    ampu = dx * sqrt( dsp_x * dsp_x + dsp_y * dsp_y ) /  dsp_x / 10.0;
                }else{
                    umin = ( u0 - dx * pm[0] / 2 ) / ( f0 - dx * pm[2] / 2 );
                    umax = ( u0 + dx * pm[0] / 2 ) / ( f0 + dx * pm[2] / 2 );
                    ampu = dx * sqrt( dsp_x * dsp_x + dsp_y * dsp_y ) /  dsp_y / 10.0;
                }
                
                // swap if u0min is greater then u0max
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
                
                v0 = v0 / f0;
                dudz = pm[6] / f0;
                dvdz = pm[7] / f0;
                
                
                
                for ( iz = 0; iz < nz ; iz ++ ){
                    
                    if ( image[iz*nxy+iy*nx+ix] == 0 ){
                        continue;
                    }
                    
                    voxel = ampu * image[iz*nxy+iy*nx+ix];
                    
                    if( dvdz > 0){
                        vmin = v0 + ( (float) iz - 0.5 ) * dz * dvdz;
                        vmax = v0 + ( (float) iz + 0.5 ) * dz * dvdz;
                    }else{
                        vmin = v0 + ( (float) iz + 0.5 ) * dz * dvdz;
                        vmax = v0 + ( (float) iz - 0.5 ) * dz * dvdz;
                    }
                    
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
                                weights_v[iv-ivmin] =  ( (float)vmax - ivmax) ;
                            } else {
                                weights_v[iv-ivmin] = 1.0;
                            }
                        }
                    }
                    
                    // projection
                    for ( iu = iumin; iu <= iumax; iu++ ){
                        for ( iv = ivmin; iv <= ivmax; iv++ ){
                            proj_view[iu*nv+iv] += weights_u[iu-iumin] * weights_v[iv-ivmin] * voxel;
                        }//iv
                        
                    }//iu
                    
                } //iz
                
            } //ix
            
        } //iy
        
        
        // Scale the projection view with z amplitude function
        for( iu = 0; iu < nu; iu ++ ){
            
            u = ( iu +  wu + offset_uv[ 2 * ib] ) * du;
            uu = sdd * sdd + u * u;
            
            for( iv = 0; iv < nv; iv ++){
                
                v = ( iv +  wv + offset_uv[ 2 * ib] ) * dv;
                ampv = sqrt( ( uu + v * v ) / uu );
                
                proj_view[iu*nv+iv] = proj_view[iu*nv+iv] * ampv;
                
            }
        }
        
    } //ib
    
    return true;
}

