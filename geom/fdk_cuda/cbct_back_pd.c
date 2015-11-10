

#define TEXTURE_MEMORY FALSE

#if TEXTURE_MEMORY
texture< float, 2, cudaReadModeElementType> texRef;
#endif

//////////////////////////////////////////////////////////////////////////////////////
///         kernel_1
//////////////////////////////////////////////////////////////////////////////////////
static __global__ void fdk_back_kernel_1(
        float sin_a,
        float cos_a,
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
        float *sval,
        int detector_type
        )
{
    // index into image array
    // determine the index of x, y
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    
    float xc, yc, xRot, yRot, rOff;
    float w2, w2s, u;
    
    // if index is out of bound
    if (ix >= nx || iy >= ny )
        return;
    
    
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
        
        w2 = sdd / ( sad + yRot );
        
        w2s = sdd * sdd  / ( ( sad + yRot )*( sad + yRot ) + (xRot - rOff) * (xRot - rOff) );
        
        
    }
    
    int ncache = 3 * ( ix + iy * nx );
    
    sval[ ncache ] = (float)u;
    sval[ ncache + 1] = (float)w2;
    sval[ ncache + 2] = (float)w2s;
}

//////////////////////////////////////////////////////////////////////////////////////
///         kernel_2
//////////////////////////////////////////////////////////////////////////////////////
#if TEXTURE_MEMORY
__global__ void fdk_back_kernel_2(
        int nx,
        int nz,
        int nxy,
        float dz,
        float wz,
        int nu,
        int nv,
        float du,
        float dv,
        float wu,
        float wv,
        float couch_z,
        float *sval,
        float *image,
        bool *map,
        bool offBoundary )
{
    
    int ix = blockIdx.x ;
    int iy = blockIdx.y ;
    int iz = threadIdx.x;
    
    //NEW 5: use shared memory
    __shared__ float u;
    __shared__ float w2;
    __shared__ float w2s;
    
    int ncache = 3 * ( ix + iy * nx );
    if ( iz == 1 ){
        u = sval[ ncache ];
        w2 = sval[ ncache + 1];
        w2s = sval[ ncache + 2];
    }
    
    __syncthreads();
    
    
    if ( u < 1E-6 || u > nu - 1 - 1E-6 ) {
        return;
    }
    
    float zc = ( iz + wz ) * dz + couch_z;
    float v = zc * w2 / dv - wv;
    
    if( offBoundary ){
        if ( v < 1E-6 || v > nv - 1 - 1E-6 ) {
            return;
        }
    }
    
    //New: use texture memory for bilinear interpolation
    image[iz * nxy + iy * nx + ix] += w2s * tex2D(texRef, v + 0.5 , u + 0.5 );
}

#else

__global__ void fdk_back_kernel_2(
        int nx,
        int nz,
        int nxy,
        float dz,
        float wz,
        int nu,
        int nv,
        float du,
        float dv,
        float wu,
        float wv,
        float couch_z,
        float *sval,
        float *image,
        float *proj,
        bool *map,
        bool offBoundary )
{
    
    int ix = blockIdx.x ;
    int iy = blockIdx.y ;
    int iz = threadIdx.x;
    
    if( ! map[iy*nx+ix] ){
        return;
    }
    
    //NEW 5: use shared memory
    __shared__ float u;
    __shared__ float w2;
    __shared__ float w2s;
    
    int ncache = 3 * ( ix + iy * nx );
    if ( iz == 1 ){
        u = sval[ ncache ];
        w2 = sval[ ncache + 1];
        w2s = sval[ ncache + 2];
    }
    
    __syncthreads();
    
    
    if ( u < 1E-6 || u > nu - 1 - 1E-6 ) {
        return;
    }
    
    int iu = floorf( u );
    float wul = u - iu;
    
    float zc = ( iz + wz ) * dz + couch_z;
    float v = zc * w2 / dv - wv;
    
    if( offBoundary ){
        if ( v < 1E-6 || v > nv - 1 - 1E-6 ) {
            return;
        }
    }else{
        
        if( v < 0.5 )
            v = 0.5;
        
        if( v > nv - 1.5 )
            v = nv - 1.5;
    }
    
    int  iv = floorf( v );
    float wvl = v - iv;
    
    int np = iu * nv + iv;
    float p1 = (1-wvl)*proj[np] + wvl * proj[np+1];
    float p2 = (1-wvl)*proj[np+nv] + wvl * proj[np+nv+1];
    
    // image[iz * nxy + iy * nx + ix] +=  w2s * ( (1-wul) * p1 + wul * p2 );
    image[iz +  nz * ( iy * nx + ix )] +=  w2s * ( (1-wul) * p1 + wul * p2 );
}

#endif



/*
 *
 * pixel driven backprojection gpu version
 *
 *
 */
bool cbct_back_pd_gpu(
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
    
    int	nxy     = nx * ny;
    int nxyz    = nx * ny * nz;
    int nuv     = nu * nv;
    
    //decide the block and grid sturcture
    dim3 dimBlock(16, 16);
    dim3 dimGrid( iDivUp( nx, dimBlock.x), iDivUp( ny, dimBlock.y));
    dim3 dimBlock2( nz, 1, 1);
    dim3 dimGrid2(nx, ny);
    
    
    
    float *dev_image;
    bool *dev_map;
    float *dev_sval;
    
    // allocate the memory on the GPU
    HANDLE_ERROR( cudaMalloc( (void**) & dev_image, nxyz * sizeof( float ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_map, nxy * sizeof( bool ) ) );
    HANDLE_ERROR( cudaMalloc( (void**) & dev_sval, 3 * nxy * sizeof( float ) ) );
    
    // initialize image volume to 0 in the GPU
    HANDLE_ERROR( cudaMemset( dev_image, 0, nxyz * sizeof(float) ) );
    
    // initialize cache volume to 0 in the GPU
    HANDLE_ERROR( cudaMemset( dev_sval, 0, 3 * nxy * sizeof(float) ) );
    
    // copy the arrays ‘map' to the GPU
    HANDLE_ERROR( cudaMemcpy( dev_map, map , nxy * sizeof(bool), cudaMemcpyHostToDevice ) );
    
#if TEXTURE_MEMORY
    // use texture memory
    cudaArray *dev_proj;
    cudaChannelFormatDesc input_tex = cudaCreateChannelDesc<float>();
    HANDLE_ERROR(  cudaMallocArray(&dev_proj, &input_tex, nv, nu) );
    
    texRef.addressMode[0]	= cudaAddressModeClamp;
    texRef.addressMode[1]	= cudaAddressModeClamp;
    texRef.filterMode	= cudaFilterModeLinear;
    texRef.normalized	= 0;
#else
    float *dev_proj;
    HANDLE_ERROR( cudaMalloc( (void**) & dev_proj, nuv * sizeof( float ) ) );
#endif
    
    int ib;
    float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
    float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
    
    float wu = - ( nu - 1.0 ) / 2.0 + offset_u;
    float wv = - ( nv - 1.0 ) / 2.0 + offset_v;
    
    
    for ( ib = 0; ib < noviews ; ib++ ) {
        
        float wz = - ( nz - 1.0 ) / 2.0 + offset_z ;
        
        float couch_z = couchz[ib];
        double beta = betas[ib];
        float sin_a = (float) sin( beta );
        float cos_a = (float) cos( beta );
        
#if TEXTURE_MEMORY
        // copy the arrays projection on the GPU
        cudaMemcpyToArray(dev_proj, 0, 0, proj + ib * nuv, nuv*sizeof(float), cudaMemcpyHostToDevice );
        cudaBindTextureToArray(texRef, dev_proj, input_tex);
#else
        HANDLE_ERROR( cudaMemcpy( dev_proj, proj + ib * nuv, nuv * sizeof(float ), cudaMemcpyHostToDevice ) );
#endif
        
        fdk_back_kernel_1<<<dimGrid, dimBlock>>>(
                sin_a,
                cos_a,
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
                dev_sval,
                detector_type );
        
#if TEXTURE_MEMORY
        
        fdk_back_kernel_2<<<dimGrid2, dimBlock2>>>(
                nx,
                nz,
                nxy,
                dz,
                wz,
                nu,
                nv,
                du,
                dv,
                wu,
                wv,
                couch_z,
                dev_sval,
                dev_image,
                dev_map,
                offBoundary );
#else
        fdk_back_kernel_2<<<dimGrid2, dimBlock2>>>(
                nx,
                nz,
                nxy,
                dz,
                wz,
                nu,
                nv,
                du,
                dv,
                wu,
                wv,
                couch_z,
                dev_sval,
                dev_image,
                dev_proj,
                dev_map,
                offBoundary );
        
#endif
        
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
    cudaFree( dev_sval );
    cudaFree( dev_proj );
    
    
    return true;
}


/*
 *
 *pixel driven backprojection
 *
 *
 */
bool cbct_back_pd_cpu(
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
                    u = xRot * w2 / du - wu;
                    
                    w2 = sdd / ( sad + yRot );
                    w2s = w2 * w2;
                    
                    
                }else{
                    
                    rOff = offset_u * du * sad / sdd ;
                    u = atan( ( xRot - rOff  ) / ( sad + yRot  ) ) * sdd / du - wu;
                    
                    w2 = sdd / ( sad + yRot );
                    w2s = sdd * sdd  / ( ( sad + yRot )*( sad + yRot ) + (xRot - rOff) * (xRot - rOff) );
                    
                    
                    
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
                            + w2s * ( (1-wul) * ( (1-wvl)*proj[npixel] + wvl * proj[npixel+1] )
                            + wul * ( (1-wvl)*proj[npixel+nv] + wvl * proj[npixel+nv+1] ) );
                    
                    
                } //iz
                
            } //ix
            
        } //iy
        
    } //ib
    
    return true;
}
