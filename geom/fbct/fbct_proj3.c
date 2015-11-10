bool fbct_proj_rd(
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
/* last updated by He Yang, Nov 2014 */
		float wu = ( - ( nu - 1.0 ) / 2.0 + offset_u ) * du;
		float cubemin_x = ( - ( nx - 1.0 ) / 2.0 + offset_x ) * dx;
		float cubemax_x = ( + ( nx - 1.0 ) / 2.0 + offset_x ) * dx;
		float cubemin_y = ( - ( ny - 1.0 ) / 2.0 + offset_y ) * dy;
		float cubemax_y = ( + ( ny - 1.0 ) / 2.0 + offset_y ) * dy;
		float delta = dx / 2.0 ;

		int ib, ix, iy;
		float temp, temp1;
		float x, y, wxl, wyl, t;

//		printf("nx = %d, dx = %f, offset_x = %f\n", nx, dx, offset_x);
//		printf("cubemin_x = %f \n", cubemin_x);
		for (  ib = 0 ; ib < noviews; ib ++ ) {

			double beta = betas[ib];
			float sin_a = (float) sin( beta );
			float cos_a = (float) cos( beta );

			float source_x = - sad * sin_a;
			float source_y = - sad * cos_a;

			int iu ;
			for ( iu = 0 ; iu < nu ; iu++ ) {

				float u =  iu * du + wu;
				float spd, raydir_x, raydir_y, tn, tf;
				double phi;

                // Flat panel
				if( detector_type == 1 ){
					spd = sqrt( sdd * sdd + u * u ); // source to pixel array distance
					phi = beta + atan( (double)( u / sdd ) );
				}//endif type1


                // Standard arc panel
                if( detector_type == 0 ){
					spd = sdd; // source to pixel array distance
					phi = beta + (double)( u / sdd );
				}//endif type0


                // Arc panel, along the detector trajectory
                if( detector_type == 2 ){
                      double alpha, gamma;
                      alpha = (double)(u/(sdd-sad));
                      gamma = (double)((sdd-sad)*sin(alpha)/(sad + (sdd-sad)*cos(alpha)));
					  phi = beta + atan(gamma);
                }//endif type2


				// Direction of the ray
				raydir_x = (float) sin( phi );
				raydir_y = (float) cos( phi );


				// find the enter and leave point of ray
				tn = - 1E10;
				tf = + 1E10;

				//if (ib==0 && iu==0){
				//	printf("Before computation\n");
				//	printf("phi=%f\n",phi);
				//	printf("tn=%f\n",tn);
				//	printf("raydir_x = %f\n", raydir_x);
				//}


				// accumulate attenuation values, nov2014
				proj[ib*nu + iu] = 0.;


				if ( raydir_x > -1E-12 && raydir_x < 1E-12 ) {
					// The ray misses the box from approximately the y-direction 
					if ( source_x < cubemin_x || source_x > cubemax_x ) {
						continue;
					}
				}
				else {

					float t1 = ( cubemin_x - source_x ) / raydir_x ; // Distance to approach cubemin_x
					float t2 = ( cubemax_x - source_x ) / raydir_x ; // Distance to approach cubemax_x

					//if (ib==0 && iu==0){
				    //    printf("proj[0]=%f\n",proj[0]);
				    //    printf("tn=%f\n",tn);
				 	//    printf("tf=%f\n",tf);
					//	printf("cubemin_x = %f\n", cubemin_x);
					//	printf("cubemax_x = %f\n", cubemax_x);
					//	printf("source_x = %f\n", source_x);
					//    printf("t1=%f\n",t1);
					//    printf("t2=%f\n",t2);
				    //}


					// t1: distance from source to the close intersection, t2: at far intersection 
					if ( t1 > t2 ) {
						float ttmp = t1;
						t1 = t2;
						t2 = ttmp;
					}

					// Choose the largest distance at the close intersection
					if (t1 > tn ){
						tn = t1;
					}

					// Choose the smallest distance at the far intersection
					if ( t2 < tf ){
						tf = t2;
					}

					// When the ray misses the object
					if (tn > tf || tf < 0 ) {
						continue;
					}
				}//end of the test for x direction

				//if (ib==0 && iu==0){
				//	printf("After the discussion of x boundaries\n");
				//	printf("proj[0]=%f\n",proj[0]);
				//    printf("tn=%f\n",tn);
				//	printf("tf=%f\n",tf);
				//	printf("t1=%f\n",t1);
				//	printf("t2=%f\n",t2);
				//}

				if ( raydir_y > -1E-12 && raydir_y < 1E-12  ) {
					if ( source_y < cubemin_y || source_y > cubemax_y ) {
						continue;
					}
				}
				else {
					float t1 = (cubemin_y - source_y ) / raydir_y ; // Distance to approach cubemin_y
					float t2 = (cubemax_y - source_y ) / raydir_y ; // Distance to approach cubemax_y

					//if (ib==0 && iu==0){
					//	printf("\n t1 = %f and t2 = %f \n", t1, t2);
					//}

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
				}//end of the test for y direction

				//if (ib==0 && iu==0){
				//	printf("proj[0]=%f\n",proj[0]);
				//    printf("tn=%f\n",tn);
				//	printf("tf=%f\n",tf);
				//	printf("t1=%f\n",t1);
				//	printf("t2=%f\n",t2);
				//}

				temp = 0;
				// go along the x-ray
				for ( t = tn ; t < tf ; t = t + delta) {

					x = ( source_x + t * raydir_x - cubemin_x ) / dx ;
					y = ( source_y + t * raydir_y - cubemin_y ) / dy ;   


					ix = floorf( x );
					iy = floorf( y );

					if( ix > 0 &&  ix < nx - 1 && iy > 0 && iy < ny-1 ){

						wxl = x - ix;
						wyl = y - iy;

						//bilinear interpolation
						temp = temp + (1-wyl) * ( (1-wxl) * image[iy*nx+ix] + wxl * image[iy*nx+ix+1]) + wyl * ( (1-wxl) * image[iy*nx+nx+ix] + wxl * image[iy*nx+nx+ix+1] ) ;
					}
				}//end of for-loop for t
				proj[ib * nu + iu ] = temp*delta/10;
}//end of for-loop for iu
}//end of for-loop for ib


		return true;
}



bool fbct_proj_dd(
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

					if( detector_type == 1 ){
						if ( cos_phi > sin_phi) {

							amplitude = dx / cos_phi / 10.0;
							umin = sdd_du * ( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 ) - wu;
							umax = sdd_du * ( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 ) - wu;

						} else {

							amplitude = dx / sin_phi / 10.0;
							umin = sdd_du * ( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 ) - wu;
							umax = sdd_du * ( tproj0 - dy_sin_2 ) / ( dsxy0 - n_dy_cos_2 ) - wu;

						}
					}else{
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

					// projection
					if (iumin == iumax) {

						proj[ib * nu + iumin ] = proj[ib * nu + iumin ] + amplitude * (umax-umin) * image[ iy*nx+ix];

					}
					else{

						for ( iu = iumin ; iu <= iumax ; iu++) {

							if ( iu == iumin ) {
								proj[ib * nu + iu ] = proj[ib * nu + iu ] + amplitude * (iumin+1.0-umin) * image[iy*nx+ix];
							} else if (iu == iumax) {
								proj[ib * nu + iu ] = proj[ib * nu + iu ] + amplitude * (umax-iumax) * image[iy*nx+ix];
							} else {
								proj[ib * nu + iu ] = proj[ib * nu + iu ] + amplitude * image[iy*nx+ix];
							}

						} //iu

					}

				} //ix

			} //iy

		} //ib

		return true;
}


bool fbct_proj_tf(
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
		float sdd_du = sdd / du;

		int ib, ix, iy, iu;
		float voxel;

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
			float bl, br, xl, xr, umin, umax;


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

					// decide if it outside the detector
					umin = max( tau[0] + 0.5, 0  );
					umax = min( tau[3] + 0.5, nu );

					if ( umin > umax ){
						continue;
					}

					iumin = floorf( umin );
					iumax = floorf( umax );
					iumax = min( iumax, nu - 1 );

					voxel = amplitude * image[ iy * nx + ix ];

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

						proj[ib * nu + iu ] = proj[ib * nu + iu ] +  weight * voxel;

					} //iu

				} //ix

			} //iy

		} //ib

		return true;
}
