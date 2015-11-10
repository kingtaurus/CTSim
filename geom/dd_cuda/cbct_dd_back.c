#define FOOTPRINT_SIZE_U 8
#define FOOTPRINT_SIZE_V 8

////////////////////////////////////////////////
///         kernel_1
////////////////////////////////////////////////
static __global__ void back_dd_kernel_1(
        double beta,
        int nx,
        int ny,
        float dx,
        float dy,
        float wx,
        float wy,
        float sad,
        float add,
        float sdd,
        int nu,
        float du,
        float wu,
        float offset_u,
        float *dev_fpu,
        int *dev_iumins,
        int *dev_iumaxs,
        float *dev_magns,
        int detector_type)
{
    
    
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    
    // if index is out of bound
    if (ix >= nx || iy >= ny )
        return;
    
    float sin_a = (float) sin( beta );
    float cos_a = (float) cos( beta );
    
    float sdd_du = sdd / du;
    
    float dx_sin_2 = dx * sin_a / 2.0;
    float dx_cos_2 = dx * cos_a / 2.0;
    float dy_sin_2 = dy * sin_a / 2.0;
    float n_dy_cos_2 = - dy * cos_a / 2.0;
    
    float xc, yc, tproj0, dsxy0, tproj0_dsxy0;
    float cos_gamma, sin_gamma, cos_phi, sin_phi, ampu;
    float umin, umax;
    int iu, iumin, iumax;
    float weights_u[FOOTPRINT_SIZE_U];
    
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
    
    iumin = floorf( umin );
    iumax = floorf( umax );
    iumax = min( iumax, nu - 1 );
    iumax = min( iumax, iumin + FOOTPRINT_SIZE_U - 1 );
    
    // set as out side fov
    if ( umin > umax ){
        iumin = nu;
    }
    
    dev_iumins[ ix + iy * nx ] = iumin;
    dev_iumaxs[ ix + iy * nx ] = iumax;
    dev_magns[ix + iy * nx  ] = sdd / dsxy0;
    
    // set as out side fov
    if ( umin > umax ){
        return;
    }
    
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
    
    // store the data
    int ncache = FOOTPRINT_SIZE_U * ( ix + iy * nx );
    //int ncache = ( ix + iy * nx );
    
    for( iu = 0; iu < FOOTPRINT_SIZE_U; iu ++ ){
        dev_fpu[ ncache + iu ] = ampu * weights_u[iu];
    }
    
    return;
    
}


////////////////////////////////////////////////
///         kernel_2
////////////////////////////////////////////////
static __global__ void back_dd_kernel_2(
        int nx,
        int ny,
        int nz,
        int nu,
        int nv,
        float dv,
        float wv,
        float dz,
        float wz,
        int *dev_iumins,
        int *dev_iumaxs,
        float *dev_magns,
        float *dev_fpu,
        float *dev_image,
        float *dev_proj,
        bool *dev_map){
    
    int iz = threadIdx.x;
    int ix = blockIdx.x;
    int iy = blockIdx.y;
    
    int index = iy * nx + ix;
    
    if( !dev_map[index] ){
     return;   
    }
    
    //se shared memory
    __shared__ int iumin;
    __shared__ int iumax;
    __shared__ float mag;
    __shared__ float weights_u[FOOTPRINT_SIZE_U];
    
    
    float zinc, vstart, vmin, vmax;
    int ivmin, ivmax, iv, iu;
    float voxel;
    
    
    if( iz == 0 ){
        iumin   = dev_iumins[ index ];
        iumax   = dev_iumaxs[ index ];
        mag     = dev_magns[ index ];
        
        for( iu = 0; iu  <= iumax - iumin ; iu++  ){
            weights_u[iu] = dev_fpu[ index * FOOTPRINT_SIZE_U + iu];
        }
    }
    
    __syncthreads();
    
    zinc = dz / dv * mag;
    
    vstart = ( wz - 0.5 ) * zinc - wv;
    
    vmin = vstart + iz * zinc;
    vmax = vmin + zinc;
    
    vmin = max( vmin + 0.5, 0  );
    vmax = min( vmax + 0.5, nv );
    
    if ( vmin > vmax ){
        return;
    }
    
    ivmin = floorf( vmin );
    ivmax = floorf( vmax );
    ivmax = min( ivmax, nv-1 );
    ivmax = min( ivmax, ivmin + FOOTPRINT_SIZE_V - 1 );
    
    voxel = 0;
    
    for ( iu = iumin; iu <= iumax; iu++ ){
        
        float voxel_v = 0;
        
        //vertical footprint
        if (ivmin == ivmax) {
            voxel_v += dev_proj[ iu * nv + ivmin ];
        }
        else if(ivmin == ivmax - 1){
            voxel_v += dev_proj[ iu * nv + ivmin ] * ( (float)ivmin + 1.0f - vmin) + dev_proj[ iu * nv + ivmax ] * ( vmax - (float)ivmax );
        }
        else{
            for ( iv = ivmin ; iv <= ivmax ; iv++) {
                if ( iv == ivmin ) {
                    voxel_v += dev_proj[ iu * nv + ivmin ] * ( (float)ivmin + 1.0f - vmin) ;
                } else if (iv == ivmax) {
                    voxel_v += dev_proj[ iu * nv + ivmax ] * ( vmax - (float)ivmax );
                } else {
                    voxel_v += dev_proj[ iu * nv + iv ];
                }
            }
        }
        
        voxel += weights_u[iu - iumin] * voxel_v;
        
        
    }//iu
    
    dev_image[ index * nz + iz ] += voxel;
    
    return;
}


////////////////////////////////////////////////
///         kernel_3
////////////////////////////////////////////////
static __global__ void back_dd_kernel_3(
        int nv,
        float du,
        float dv,
        float wu,
        float wv,
        float sdd,
        float * dev_proj,
        bool detector_type )
{
    
    int iu = blockIdx.x;
    int iv = threadIdx.x;
    
    float   u = ( iu +  wu ) * du;
    float   uu = sdd * sdd + u * u;
    float   v = ( iv +  wv ) * dv;
    float   ampv;
    
    if( detector_type == 1 ){ //flat panel detector
        ampv = sqrt( ( uu + v * v ) / uu );
    }else{ //curved detector
        ampv = sqrt( sdd * sdd  + v * v ) / sdd;
    }
    
    dev_proj[iu*nv+iv] = dev_proj[iu*nv+iv] * ampv;
    
    return;
    
}


bool cbct_back_dd_gpu(
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
    
    int	nxy     = nx * ny;
    int nxyz    = nx * ny * nz;
    int nuv     = nu * nv;
    
    //decide the block and grid sturcture
    dim3 dimBlock1(16, 16);
    dim3 dimGrid1( iDivUp( nx, dimBlock1.x), iDivUp( ny, dimBlock1.y));
    
    dim3 dimBlock2( nz, 1, 1);
    dim3 dimGrid2(nx, ny);
    
    dim3 dimBlock3(nv, 1);
    dim3 dimGrid3( nu, 1);
    
    
    float *dev_image;
    bool *dev_map;
    float *dev_fpu;
    int *dev_iumins;
    int *dev_iumaxs;
    float *dev_proj;
    float *dev_magns;;
    
    // allocate the memory on the GPU
    HANDLE_ERROR( cudaMalloc( (void**) & dev_image, nxyz * sizeof( float ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_map, nxy * sizeof( bool ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_fpu, FOOTPRINT_SIZE_U * nxy * sizeof( float ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_iumins, nxy * sizeof( int ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_iumaxs, nxy * sizeof( int ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_proj, nuv * sizeof( float ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_magns, nxy * sizeof( float ) ) );
    
    HANDLE_ERROR( cudaMemcpy( dev_map, map , nxy * sizeof(bool), cudaMemcpyHostToDevice ) );
    HANDLE_ERROR( cudaMemset( dev_image, 0, nxyz * sizeof( float )  ) );
    
    int ib;
    
    float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
    float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
    
    float wu = - ( nu - 1.0 ) / 2.0 + offset_u;
    float wv = - ( nv - 1.0 ) / 2.0 + offset_v;
    
    for ( ib = 0; ib < noviews ; ib++ ) {
        
        float wz = - ( nz - 1.0 ) / 2.0 + offset_z + ( couchz[ib] / dz );
        double beta = betas[ib];
  
        // GPU kernel for computing the u foot print
        back_dd_kernel_1<<<dimGrid1, dimBlock1>>>(
                beta,
                nx,
                ny,
                dx,
                dy,
                wx,
                wy,
                sad,
                add,
                sdd,
                nu,
                du,
                wu,
                offset_u,
                dev_fpu,
                dev_iumins,
                dev_iumaxs,
                dev_magns,
                detector_type);
        
        float *proj_view = proj + nuv * ib;
        HANDLE_ERROR( cudaMemcpy( dev_proj, proj_view, nuv * sizeof(float), cudaMemcpyHostToDevice ) );
        
        // Scale the projection view with z amplitude function
        back_dd_kernel_3<<<dimGrid3, dimBlock3>>>(
                nv,
                du,
                dv,
                wu,
                wv,
                sdd,
                dev_proj,
                detector_type );

        back_dd_kernel_2<<<dimGrid2, dimBlock2>>>(
                nx,
                ny,
                nz,
                nu,
                nv,
                dv,
                wv,
                dz,
                wz,
                dev_iumins,
                dev_iumaxs,
                dev_magns,
                dev_fpu,
                dev_image,
                dev_proj,
                dev_map);
        
    } //ib
    
    float *image_temp = (float *) malloc( nxyz * sizeof(float) );
    
    // copy the arrays ‘image' back to the CPU
    HANDLE_ERROR( cudaMemcpy( image_temp, dev_image, nxyz * sizeof(float), cudaMemcpyDeviceToHost ) );
    
    int ix, iy, iz;
    for ( iz = 0; iz < nz ; iz ++ ) {
        
        for ( iy = 0; iy < ny ; iy ++ ) {
            
            for ( ix = 0; ix < nx ; ix++ ) {
                
                image[iz * nxy + iy * nx + ix] = image_temp[iz +  nz * ( iy * nx + ix )];
                
            }
        }
    }
    
    free( image_temp );
    cudaFree( dev_image );
    cudaFree( dev_map );
    cudaFree( dev_fpu );
    cudaFree( dev_iumins );
    cudaFree( dev_iumaxs );
    cudaFree( dev_proj );
    cudaFree( dev_magns );
    
    
    
    return true;
}



bool cbct_back_dd_cpu(
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




