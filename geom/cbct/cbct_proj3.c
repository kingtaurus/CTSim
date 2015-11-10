/*
 * 
 * Ray driven forward projector
 * Special arc detector included, last updated by He Yang, Dec 2014
 *
 *
 */

bool cbct_proj_rd(
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
	double* betas,
	float* couchz,
	int detector_type,
	float* image,
	float* proj ){


		float wu = ( - ( nu - 1.0 ) / 2.0 + offset_u ) * du;
		float wv = ( - ( nv - 1.0 ) / 2.0 + offset_v ) * dv;

		float cubemin_x = ( - ( nx - 1.0 ) / 2.0 + offset_x ) * dx;
		float cubemax_x = ( + ( nx - 1.0 ) / 2.0 + offset_x ) * dx;
		float cubemin_y = ( - ( ny - 1.0 ) / 2.0 + offset_y ) * dy;
		float cubemax_y = ( + ( ny - 1.0 ) / 2.0 + offset_y ) * dy;


		float delta = dx / 2.0 ;
		int ib, ix, iy, iz, iu, iv, npixel, nvoxel;
		float x, y, z, wxl, wyl, wzl, t, temp;
		int nxy = nx * ny ;

		for (  ib = 0 ; ib < noviews; ib ++ ) {   // 1st loop: # of view


			float cubemin_z = ( - ( nz - 1.0 ) / 2.0 + offset_z ) * dz + couchz[ib];
			float cubemax_z = ( + ( nz - 1.0 ) / 2.0 + offset_z ) * dz + couchz[ib];


			double beta = betas[ib];
			float sin_a = (float) sin( beta );
			float cos_a = (float) cos( beta );

			float source_x = - sad * sin_a;
			float source_y = - sad * cos_a;
			float source_z = 0.0;


			for ( iu = 0 ; iu < nu ; iu++ ) { // 2nd loop: detector coordinate u

				float u =  iu * du + wu;
				float spd, phi, raydir_x, raydir_y, tn, tf;

				// Flat Panel
				if( detector_type == 1 ){
					spd = sqrt( sdd * sdd + u * u ); // source to pixel array distance
					phi = beta + atan( (double)( u / sdd ) );
				} 
				// Standard Curved Panel
				else if( detector_type == 0 ){
					spd = sdd; // source to pixel array distance
					phi = beta + (double)( u / sdd );
				}
				// Special Arc Panel
				else if( detector_type == 2){
                      double alpha, gamma;
                      alpha = (double)(u/add);
                      gamma = (double)(add*sin(alpha)/(sad + add*cos(alpha)));
					  gamma = atan(gamma);
					  phi = beta + gamma;

					  spd=(add*sin(alpha))*(add*sin(alpha))+(sad+add*cos(alpha))*(sad+add*cos(alpha));
					  spd = sqrt(spd);   // Source to detector cell distance
				}

				raydir_x = (float) sin( phi );
				raydir_y = (float) cos( phi );

				// find the enter and leave point of ray
				tn = - 1E10;
				tf = + 1E10;

				if ( raydir_x > -1E-12 && raydir_x < 1E-12 ) {
					if ( source_x < cubemin_x || source_x > cubemax_x ) {
						continue;
					}
				}
				else {

					float t1 = ( cubemin_x - source_x ) / raydir_x ;
					float t2 = ( cubemax_x - source_x ) / raydir_x ;

					if ( t1 > t2 ) {
						float ttmp = t1;
						t1 = t2;
						t2 = ttmp;
					}

					if (t1 > tn ){
						tn = t1;
					}

					if ( t2 < tf ){
						tf = t2;
					}

					if (tn > tf || tf < 0 ) {
						continue;
					}

				}

				if ( raydir_y > -1E-12 && raydir_y < 1E-12  ) {
					if ( source_y < cubemin_y || source_y > cubemax_y ) {
						continue;
					}
				}
				else {
					float t1 = (cubemin_y - source_y ) / raydir_y ;
					float t2 = (cubemax_y - source_y ) / raydir_y ;

					if ( t1 > t2 ) {
						float ttmp = t1;
						t1 = t2;
						t2 = ttmp;
					}

					if (t1 > tn ){
						tn = t1;
					}

					if ( t2 < tf ){
						tf = t2;
					}

					if (tn > tf || tf < 0 ) {
						continue;
					}
				}

				npixel = ib*nu*nv + iu*nv;

				for ( iv = 0; iv < nv ; iv ++ ){  //3rd loop: det coordinate v

					float raydir_z = (  iv * dv + wv - source_z ) / spd;

					temp = 0;
					// go along the x-ray
					for ( t = tn ; t < tf ; t = t + delta) {

						z = ( source_z + t * raydir_z - cubemin_z ) / dz ;

						if ( z < 1 ){
							continue;
						}

						if ( z > nz - 1.1 ){
							break;
						}

						x = ( source_x + t * raydir_x - cubemin_x ) / dx ;
						y = ( source_y + t * raydir_y - cubemin_y ) / dy ;


						ix = floorf( x );
						iy = floorf( y );
						iz = floorf( z );

						wzl = z - iz;

						nvoxel = iz * nxy + iy * nx + ix;

						if( ix > 0 && ix < nx - 1 && iy > 0 &&  iy < ny - 1 ){

							wxl = x - ix;
							wyl = y - iy;

							//trilinear interpolation
							temp = temp + (1-wzl) *(  (1-wyl) * ( (1-wxl) * image[nvoxel] + wxl * image[nvoxel+1]) + wyl * ( (1-wxl) * image[nvoxel+nx] + wxl * image[nvoxel+nx+1] ) ) 
								+ wzl * ( (1-wyl) * ( (1-wxl) * image[nvoxel+nxy] + wxl * image[nvoxel+nxy+1]) + wyl * ( (1-wxl) * image[nvoxel+nxy+nx] + wxl * image[nvoxel+nxy+nx+1] ) );
						}
					} //t

					proj[npixel + iv] = temp * delta * sqrt( 1 + raydir_z * raydir_z )  / 10;


				} //iv

			} //iu

		} //ib

		return true;
}

#define FOOTPRINT_SIZE_U 10
#define FOOTPRINT_SIZE_V 10
/*
 *
 * distance driven projection with binary image map
 * 
 * Special arc detector included, last updated by He Yang, Jan 2015
 *
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
        float umin, umax, vmin, vmax;
        int iumin, iumax, ivmin, ivmax, iz_start, iz_stop;
		float gam_min, gam_max, xi_min, xi_max;
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
                
                if( detector_type == 1 ){ //Flat Panel
                    if ( cos_phi > sin_phi) {
                        
                        ampu = dx / cos_phi / 10.0;
                        umin = sdd_du * ( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 ) - wu;
                        umax = sdd_du * ( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 ) - wu;
                        
                    } else {
                        
                        ampu = dx / sin_phi / 10.0;
                        umin = sdd_du * ( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 ) - wu;
                        umax = sdd_du * ( tproj0 - dy_sin_2 ) / ( dsxy0 - n_dy_cos_2 ) - wu;
                        
                    }
                }else if( detector_type == 0){ // Standard curved panel
                    if ( cos_phi > sin_phi) {
                        
                        ampu = dx / cos_phi / 10.0;
                        umin = sdd_du * atan(( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 )) - wu;
                        umax = sdd_du * atan(( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 )) - wu;
                        
                    } else {
                        
                        ampu = dx / sin_phi / 10.0;
                        umin = sdd_du * atan(( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 )) - wu;
                        umax = sdd_du * atan(( tproj0 - dy_sin_2 ) / ( dsxy0 - n_dy_cos_2 )) - wu;
                    }
                }else if( detector_type == 2){
					if ( cos_phi > sin_phi) {

						ampu = dx / cos_phi / 10.0;

						gam_min = atan(( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 ));
						xi_min = gam_min + asin(sad/add*sin(gam_min));
						umin = xi_min*add/du  - wu;

						gam_max = atan(( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 ));
						xi_max = gam_max + asin(sad/add*sin(gam_max));
						umax = xi_max*add/du  - wu;
					} else {

						ampu = dx / sin_phi / 10.0;

						gam_min = atan(( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 ));
						xi_min = gam_min + asin(sad/add*sin(gam_min));
						umin = xi_min*add/du  - wu;

						gam_max = atan(( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 ));
						xi_max = gam_max + asin(sad/add*sin(gam_max));
						umax = xi_max*add/du  - wu;
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
                    
                    if ( image[iz*nxy+iy*nx+ix] == 0 ){
                        continue;
                    }
                    
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
            
            u = ( iu +  wu ) * du;
            uu = sdd * sdd + u * u;
            
            for( iv = 0; iv < nv; iv ++){
                
                v = ( iv +  wv ) * dv;
                
                if( detector_type == 1 ){ //flat panel detector
                    ampv = sqrt( ( uu + v * v ) / uu );
                }else if( detector_type == 0){ //curved detector
                    ampv = sqrt( sdd * sdd  + v * v ) / sdd;
                }else if( detector_type == 2){ // Special Arc Detector
					ampv = (sad+add*cos(u/add))*(sad+add*cos(u/add))+(add*sin(u/add))*(add*sin(u/add));
					ampv=sqrt((ampv+v*v)/ampv);
				}
                
                proj_view[iu*nv+iv] = proj_view[iu*nv+iv] * ampv;
                
            }
        }
        
        
    } //ib
    
    return true;
}

/*
 *
 * trapesoid footprint projection with binary image map
 *
 *
 */
bool cbct_proj_tf(
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



