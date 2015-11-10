#define FOOTPRINT_SIZE_U 8
#define FOOTPRINT_SIZE_V 8
#define SUBSET_DISTANCE 4
#define CACHE_COLUMN FALSE
#define CACHE_COLUMN_LENGTH 256

////////////////////////////////////////////////
///         kernel_1
////////////////////////////////////////////////
static __global__ void proj_sf_kernel_1(
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
    float tau[4];
    float bl, br, xl, xr;
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
    
    // store the data
    int ncache = FOOTPRINT_SIZE_U * ( ix + iy * nx );
    //int ncache = ( ix + iy * nx );
    
    for( iu = 0; iu < FOOTPRINT_SIZE_U; iu ++ ){
        dev_fpu[ ncache + iu ] = ampu * weights_u[iu];
    }
    
    
}

////////////////////////////////////////////////
///         kernel_2
////////////////////////////////////////////////
static __global__ void proj_sf_kernel_2(
        int is,
        int i,
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
        bool xset){
    
    int iv = threadIdx.x;
    
    //se shared memory
    __shared__ int iumin;
    __shared__ int iumax;
    __shared__ float mag;
    __shared__ float weights_u[FOOTPRINT_SIZE_U];
    
    int iz, index, iu;
    int izmin, izmax;
    float zmax, zmin;
    float zinc, vinc, voxel;
    
    if( xset ){
        int ix = blockIdx.x * SUBSET_DISTANCE + is;
        if( ix >= nx ){
            return;
        }
        index = ix + i * nx;
    }else{
        int iy = blockIdx.x * SUBSET_DISTANCE + is;
        if( iy >= ny ){
            return;
        }
        index = i + iy * nx;
    }
    
    if( iv == 0 ){
        iumin   = dev_iumins[ index ];
        iumax   = dev_iumaxs[ index ];
        mag     = dev_magns[ index ];
        
        for( iu = 0; iu  <= iumax - iumin ; iu++  ){
            weights_u[iu] = dev_fpu[ index * FOOTPRINT_SIZE_U + iu];
        }
    }
    
    __syncthreads();
    
    zinc = dz / dv * mag;
    vinc = 1 / zinc ;
    
    zmin = ( iv + wv - 0.5 ) * vinc - wz + 0.5;
    zmax  = zmin + vinc;
    
    if( zmin > zmax )
        return;
    
    //zmin = max( zmin + 0.5, 0  );
    //zmax = min( zmax + 0.5, nz - 1E-6 );
    
    izmin = max( floorf( zmin ), 0) ;
    izmax = min( floorf( zmax ), nz - 1);
    
    //vertical footprint
    if (izmin == izmax) {
        voxel = dev_image[ izmin + index * nz ];
    }
    else if(izmin == izmax - 1){
        voxel = dev_image[ izmin + index * nz ] * ( (float)izmin + 1.0f - zmin) + dev_image[ izmax + index * nz ] * ( zmax - (float)izmax );
        voxel = voxel / vinc;
    }
    else{
        voxel = 0;
        for ( iz = izmin ; iz <= izmax ; iz++) {
            if ( iz == izmin ) {
                voxel += ( (float)izmin + 1.0 - zmin) * dev_image[ iz + index * nz ];
            } else if (iz == izmax) {
                voxel += dev_image[ iz + index * nz ] * ( zmax - (float)izmax );
            } else {
                voxel += dev_image[ iz + index * nz ];
            }
        }
        voxel = voxel / vinc;
    }
    
    for( iu = iumin ; iu <= iumax ; iu ++ ){
        dev_proj[ iu * nv + iv ] += weights_u[iu - iumin] * voxel;
    }
    
    __syncthreads();
    
    return;
}

////////////////////////////////////////////////
///         kernel_3
////////////////////////////////////////////////
static __global__ void proj_sf_kernel_3(
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
    
}

bool cbct_proj_tf_gpu(
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
    
    
    dim3 dimBlock3(nv, 1);
    dim3 dimGrid3( nu, 1);
    
    float *dev_image;
    //bool *dev_map;
    float *dev_fpu;
    int *dev_iumins;
    int *dev_iumaxs;
    float *dev_proj;
    float *dev_magns;;
    
    // allocate the memory on the GPU
    HANDLE_ERROR( cudaMalloc( (void**) & dev_image, nxyz * sizeof( float ) ) );
    //HANDLE_ERROR( cudaMalloc( (void**) & dev_map, nxy * sizeof( bool ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_fpu, FOOTPRINT_SIZE_U * nxy * sizeof( float ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_iumins, nxy * sizeof( int ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_iumaxs, nxy * sizeof( int ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_proj, nuv * sizeof( float ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_magns, nxy * sizeof( float ) ) );
    
    
    // initialize image volume to 0 in the GPU
    //HANDLE_ERROR( cudaMemset( dev_image, 0, nxyz * sizeof(float) ) );
    
    
    // change the image volume structure
    float *image_temp = (float *) malloc( nxyz * sizeof(float) );
    
    int ix, iy, iz;
    for ( iz = 0; iz < nz ; iz ++ ) {
        for ( iy = 0; iy < ny ; iy ++ ) {
            for ( ix = 0; ix < nx ; ix++ ) {
                
                image_temp[iz +  nz * ( iy * nx + ix )] = image[iz * nxy + iy * nx + ix];
                
            }
        }
    }
    // copy the arrays image to the GPU
    HANDLE_ERROR( cudaMemcpy( dev_image, image_temp, nxyz * sizeof( float ), cudaMemcpyHostToDevice ) );
    free( image_temp );
    
    
    //HANDLE_ERROR( cudaMemcpy( dev_map, map , nxy * sizeof(bool), cudaMemcpyHostToDevice ) );
    //HANDLE_ERROR( cudaMemset( dev_fpu, 0, FOOTPRINT_SIZE_U * nxy * sizeof(float) ) );
    
    int ib, is;
    
    float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
    float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
    
    float wu = - ( nu - 1.0 ) / 2.0 + offset_u;
    float wv = - ( nv - 1.0 ) / 2.0 + offset_v;
    
    
    for ( ib = 0; ib < noviews ; ib++ ) {
        
        float wz = - ( nz - 1.0 ) / 2.0 + offset_z + ( couchz[ib] / dz );
        double beta = betas[ib];
        
        HANDLE_ERROR( cudaMemset( dev_proj, 0, nuv * sizeof( float )  ) );
        //HANDLE_ERROR( cudaMemset( dev_iumins, 0, nxy * sizeof( int )  ) );
        //HANDLE_ERROR( cudaMemset( dev_iuset, 0, nu * nk * sizeof(int)  ) );
        
        float *proj_view = proj + nu * nv * ib;
        
        // GPU kernel for computing the u foot print
        proj_sf_kernel_1<<<dimGrid1, dimBlock1>>>(
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
        
        
        
        // Scale the projection view with z amplitude function
        if ( abs( sin( beta ) ) < abs( cos( beta ) ) ){
            
            
            dim3 dimBlock2(nv, 1);
            dim3 dimGrid2( iDivUp( nx, SUBSET_DISTANCE), 1);
            
            for( is = 0; is < SUBSET_DISTANCE ; is ++ ){
                
                for( iy = 0; iy < ny; iy ++ ){
                    
                    proj_sf_kernel_2<<<dimGrid2, dimBlock2>>>(
                            is,
                            iy,
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
                            true);
                    
                }
            }
            
        }else{
            
            dim3 dimBlock2(nv, 1 );
            dim3 dimGrid2( iDivUp( ny, SUBSET_DISTANCE), 1);
            
            for( is = 0; is < SUBSET_DISTANCE ; is ++ ){
                
                for( ix = 0; ix < nx; ix ++ ){
                    
                    proj_sf_kernel_2<<<dimGrid2, dimBlock2>>>(
                            is,
                            ix,
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
                            false);
                    
                }
            }
        }
        
        proj_sf_kernel_3<<<dimGrid3, dimBlock3>>>(
                nv,
                du,
                dv,
                wu,
                wv,
                sdd,
                dev_proj,
                detector_type );
        
        
        HANDLE_ERROR( cudaMemcpy( proj_view, dev_proj, nuv * sizeof(float), cudaMemcpyDeviceToHost ) );
        
    } //ib
    
    
    cudaFree( dev_image );
    //cudaFree( dev_map );
    cudaFree( dev_fpu );
    cudaFree( dev_iumins );
    cudaFree( dev_iumaxs );
    cudaFree( dev_proj );
    cudaFree( dev_magns );
    
    return true;
}



bool cbct_proj_tf_cpu(
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
    int	nxy = nx * ny;
    int ib, ix, iy, iz, iu, iv;
    float u, v, uu;
    
    // forward projection per view
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
        int iumin, iumax, ivmin, ivmax, iz_start, iz_stop;
        float weights_u[FOOTPRINT_SIZE_U], weights_v[FOOTPRINT_SIZE_V];
        float *proj_view = proj + nu * nv * ib;
        
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
                    
                    if ( image[iz*nxy+iy*nx+ix] == 0 ){
                        continue;
                    }
                    
                    //dsz0 = ( iz + wz ) * dz ;
                    
                    //ampv = sqrt(tproj0_dsxy0 * tproj0_dsxy0 + dsz0 * dsz0) / tproj0_dsxy0;
                    
                    voxel = ampu * image[iz*nxy+iy*nx+ix];
                    
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
                                weights_v[0] = (ivmin + 1.0 - vmin) ;
                            } else if (iv == ivmax) {
                                weights_v[iv-ivmin] =  (float)( vmax - ivmax) ;
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
        
        
    } //ib
    
    return true;
}


