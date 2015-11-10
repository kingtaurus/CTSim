
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
	float *image,
	float *proj ){

		float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
		float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
		float wz = - ( nz - 1.0 ) / 2.0 + offset_z;
		float wu = - ( nu - 1.0 ) / 2.0 + offset_u;
		float wv = - ( nv - 1.0 ) / 2.0 + offset_v;
		float sdd_du = sdd / du;

		int ib, ix, iy, iz, iu, iv;

		for ( ib = 0; ib < noviews ; ib++) {

			double beta = betas[ib];
			float sin_a = (float) sin( beta );
			float cos_a = (float) cos( beta );

			float dx_sin_2 = dx * sin_a / 2.0;
			float dx_cos_2 = dx * cos_a / 2.0;
			float dy_sin_2 = dy * sin_a / 2.0;
			float n_dy_cos_2 = - dy * cos_a / 2.0;

			float xc, yc, tproj0, dsxy0, tproj0_dsxy0, dsz0, zinc, vstart;
			float cos_gamma, sin_gamma, cos_phi, sin_phi, ampu, ampv, voxel;
			float umin, umax, vmin, vmax;
			int iumin, iumax, ivmin, ivmax;
			float weights_u[8], weights_v[8];
			int nxy = nx * ny;
			for ( iy = 0; iy < ny ; iy ++ ) {

				for ( ix = 0; ix < nx ; ix++ ) {

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
						umin = sdd_du * ( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 ) - wu;
						umax = sdd_du * ( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 ) - wu;

					} else {

						ampu = dx / sin_phi / 10.0;
						umin = sdd_du * ( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 ) - wu;
						umax = sdd_du * ( tproj0 - dy_sin_2 ) / ( dsxy0 - n_dy_cos_2 ) - wu;

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
					iumax = min( iumax, iumin+7 );

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

					for ( iz = 0 ; iz < nz ; iz++ ){

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
						ivmax = min( ivmax, ivmin+7 );

						// projection
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

						image[iz*nxy+iy*nx+ix] = image[iz*nxy+iy*nx+ix] +  ampu * ampv * voxel;


					} //iz

				} //ix

			} //iy

		} //ib

		return true;
}


