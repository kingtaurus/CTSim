bool fbct_back_pd( 
	int nx,
	int ny,
	float dx,
	float dy,
	float offset_x,
	float offset_y,
	int nu,
	float du,
	float offset_u,
	float sad,
	float add,
	float sdd,
	int noviews,
	double *betas,
	int detector_type,
	float *image,
	float *proj ){

		float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
		float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
		float wu = - ( nu - 1.0 ) / 2.0 + offset_u;


		float xc, yc, xRot, yRot, rOff;
		float w2, w2s, u, wul;
		int ib, ix, iy, iu;

		for ( ib = 0; ib < noviews ; ib++) {

			double beta = betas[ib];
			float sin_a = (float) sin( beta );
			float cos_a = (float) cos( beta );

			for ( ix = 0; ix < nx ; ix++ ) {

				for ( iy = 0; iy < ny ; iy ++ ) {

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

					image[iy*nx+ix] = image[iy*nx+ix] + w2s * ( (1-wul) * proj[ib*nu+iu] + wul * proj[ib*nu+iu+1] );

				} //ix

			} //iy

		} //ib

		return true;
}

// Distance-Driven Backprojector
// Special arc detector included, last updated by He Yang, Dec 2014
bool fbct_back_dd(
	int nx,
	int ny,
	float dx,
	float dy,
	float offset_x,
	float offset_y,
	int nu,
	float du,
	float offset_u,
	float sad,
	float add,
	float sdd,
	int noviews,
	double *betas,
	int detector_type,
	float *image,
	float *proj  ){


		float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
		float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
		float wu = - ( nu - 1.0 ) / 2.0 + offset_u;
		float sdd_du = sdd / du;
		int ib, ix, iy, iu;

		for ( ib = 0; ib < noviews ; ib++) {

			double beta = betas[ib];
			float sin_a = (float) sin( beta );
			float cos_a = (float) cos( beta );

			float dx_sin_2 = dx * sin_a / 2.0;
			float dx_cos_2 = dx * cos_a / 2.0;
			float dy_sin_2 = dy * sin_a / 2.0;
			float n_dy_cos_2 = - dy * cos_a / 2.0;

			float xc, yc, tproj0, dsxy0, tproj0_dsxy0;
			float cos_gamma, sin_gamma, cos_phi, sin_phi, amplitude;
			float umin, umax;
			int iumin, iumax;
			float temp;
			float gam_min, gam_max, xi_min, xi_max;

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

					if( detector_type == 1 ){    // Flat panel
						if ( cos_phi > sin_phi) {

							amplitude = dx / cos_phi / 10.0;
							umin = sdd_du * ( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 ) - wu;
							umax = sdd_du * ( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 ) - wu;

						} else {

							amplitude = dx / sin_phi / 10.0;
							umin = sdd_du * ( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 ) - wu;
							umax = sdd_du * ( tproj0 - dy_sin_2 ) / ( dsxy0 - n_dy_cos_2 ) - wu;

						}
					}
					else if( detector_type == 0){  // Standard curved detector
						if ( cos_phi > sin_phi) {

							amplitude = dx / cos_phi / 10.0;
							umin = sdd_du * atan(( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 )) - wu;
							umax = sdd_du * atan(( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 )) - wu;

						} else {

							amplitude = dx / sin_phi / 10.0;
							umin = sdd_du * atan(( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 )) - wu;
							umax = sdd_du * atan(( tproj0 - dy_sin_2 ) / ( dsxy0 - n_dy_cos_2 )) - wu;

						}
					}
					else if( detector_type == 2){   // Special arc detector 
						if ( cos_phi > sin_phi) {

							amplitude = dx / cos_phi / 10.0;

							gam_min = atan(( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 ));
							xi_min = gam_min + asin(sad/add*sin(gam_min));
							umin = xi_min*add/du  - wu;

							gam_max = atan(( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 ));
							xi_max = gam_max + asin(sad/add*sin(gam_max));
							umax = xi_max*add/du  - wu;

						} else {

							amplitude = dx / sin_phi / 10.0;

							gam_min = atan(( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 ));
							xi_min = gam_min + asin(sad/add*sin(gam_min));
							umin = xi_min*add/du  - wu;

							gam_max = atan(( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 ));
							xi_max = gam_max + asin(sad/add*sin(gam_max));
							umax = xi_max*add/du  - wu;

						}
					}
//					else
//						printf("Wrong detector type! Detector_type Should be 0, 1 or 2.");
//				}

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

					temp = 0.0;

					// projection
					if (iumin == iumax) {

						temp = (umax-umin) * proj[ib * nu + iumin ];

					}
					else{

						for ( iu = iumin ; iu <= iumax ; iu++) {

							if ( iu == iumin ) {
								temp = temp + (iumin+1.0-umin) * proj[ib * nu + iu ];
							} else if (iu == iumax) {
								temp = temp + (umax-iumax) * proj[ib * nu + iu ];
							} else {
								temp = temp + proj[ib * nu + iu ];
							}

						} //iu

					}

					image[iy*nx+ix] = image[iy*nx+ix] + amplitude * temp;

				} //ix

			} //iy

		} //ib

		return true;
}

bool fbct_back_tf( 
	int nx,
	int ny,
	float dx,
	float dy,
	float offset_x,
	float offset_y,
	int nu,
	float du,
	float offset_u,
	float sad,
	float add,
	float sdd,
	int noviews,
	double *betas,
	int detector_type,
	float *image,
	float *proj  ){

		float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
		float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
		float wu = - ( nu - 1.0 ) / 2.0 + offset_u;
		float sdd_du = sdd / du;

		int ib, ix, iy, iu;
		float voxel, temp;

		for ( ib = 0; ib < noviews ; ib++) {

			double beta = betas[ib];		
			float sin_a = (float) sin( beta );
			float cos_a = (float) cos( beta );

			float dx_sin_2 = dx * sin_a / 2.0;
			float dx_cos_2 = dx * cos_a / 2.0;
			float dy_sin_2 = dy * sin_a / 2.0;
			float n_dy_cos_2 = - dy * cos_a / 2.0;

			float xc, yc, tproj0, dsxy0, tproj0_dsxy0;
			float cos_gamma, sin_gamma, cos_phi, sin_phi, amplitude;
			float tau[4];
			int iumin, iumax;
			float bl, br, xl, xr;

			for ( ix = 0; ix < nx ; ix++ ) {

				for ( iy = 0; iy < ny ; iy ++ ) {


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
						amplitude = dx / cos_phi / 10.0;
					} else {
						amplitude = dx / sin_phi / 10.0;
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

					// boundaries
					iumin = max( floorf(tau[0] + 0.5), 0 );
					iumax = min( ceilf(tau[3] + 0.5), nu - 1);

					if (iumin > iumax) {
						continue;
					}

					temp = 0.0;

					for ( iu = iumin ; iu <= iumax ; iu++) {

						float weight = 0.0;
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

						temp = temp + amplitude * weight * proj[ib * nu + iu ];

					} //iu


					image[iy*nx+ix] = image[iy*nx+ix] + temp;

				} //ix

			} //iy

		} //ib



		return true;
}

bool fbct_back_var( 
	int nx,
	int ny,
	float dx,
	float dy,
	float offset_x,
	float offset_y,
	int nu,
	float du,
	float offset_u,
	float sad,
	float add,
	float sdd,
	int noviews,
	double *betas,
	int detector_type,
	float *image,
	float *proj ){

		float wx = - ( nx - 1.0 ) / 2.0 + offset_x;
		float wy = - ( ny - 1.0 ) / 2.0 + offset_y;
		float wu = - ( nu - 1.0 ) / 2.0 + offset_u;


		float xc, yc, xRot, yRot, rOff;
		float w2, w2s, u, wul;
		int ib, ix, iy, iu;

		for ( ib = 0; ib < noviews ; ib++) {

			double beta = betas[ib];
			float sin_a = (float) sin( beta );
			float cos_a = (float) cos( beta );

			for ( ix = 0; ix < nx ; ix++ ) {

				for ( iy = 0; iy < ny ; iy ++ ) {

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

					image[iy*nx+ix] = image[iy*nx+ix] + w2s * w2s * ( (1-wul) * (1-wul)  * proj[ib*nu+iu] + wul * wul * proj[ib*nu+iu+1] );

				} //ix

			} //iy

		} //ib

		return true;
}

