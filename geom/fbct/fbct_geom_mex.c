#include "stdio.h"
#include "string.h"
#include "mex.h"
#include "ct_def.h"


#include "fbct_proj.c"
#include "fbct_back.c"

/**************************************************************************
* mxu_string()
* caller must free using mxu_string_free()
**************************************************************************/
char *mxu_getstring( Cmx mx)
{
	char	*string;
	int	n = mxGetM(mx) * mxGetN(mx) + 1;

	string = (char *) mxCalloc(n, sizeof(char));

	if (mxGetString(mx, string, n)){
		mexErrMsgTxt("bug with mxGetString. \n");
	}

	return string;
}

void mxu_freestring(char *s)
{
	mxFree(s);
}

static void fbct_geom_mex_help(){
	mexPrintf("\n\
			  coming soon!\
			  \n");
}

/*
* static bool fbct_proj_all( ... )
* forward projection call function
*/
static bool fbct_proj_all
	(
	mxArray *plhs[],
	Cmx mx_nxyz,
	Cmx mx_dxyz,
	Cmx mx_cxyz,
	Cmx mx_nuv,
	Cmx mx_duv,
	Cmx mx_cuv,
	Cmx mx_sad,
	Cmx mx_add,
	Cmx mx_sdd,
	Cmx mx_noviews,
	Cmx mx_betas,
	Cmx mx_detector_type,
	Cmx mx_image,
	GEOM_TYPE geom_type){

#if Debug1
		mexPrintf("fbct_proj_all()... \n");
#endif

		bool sof = false;

		int nx = mxGetPr_cint32(mx_nxyz)[0];
		int ny = mxGetPr_cint32(mx_nxyz)[1];
		float dx = mxGetPr_cfloat(mx_dxyz)[0];
		float dy = mxGetPr_cfloat(mx_dxyz)[1];
		float offset_x = mxGetPr_cfloat(mx_cxyz)[0];
		float offset_y = mxGetPr_cfloat(mx_cxyz)[1];

		int nu = mxGetPr_cint32(mx_nuv)[0];
		float du = mxGetPr_cfloat(mx_duv)[0];
		float offset_u = mxGetPr_cfloat(mx_cuv)[0];

		float sad = mxGetPr_cfloat(mx_sad)[0];
		float add = mxGetPr_cfloat(mx_add)[0];
		float sdd = mxGetPr_cfloat(mx_sdd)[0];
		int noviews = mxGetPr_cint32(mx_noviews)[0];
		int detector_type = mxGetPr_cint32(mx_detector_type)[0];

		double *betas;
		float *image;
		float *proj;
		int *dim_image;
		mwSize dim_proj[2];

		/* get parameters for geometry */
		if (noviews != mxGetM(mx_betas) * mxGetN(mx_betas) || !mxIsRealDouble(mx_betas)){
			mexErrMsgTxt("Error: mx_betas must have noviews X 1 double array.\n");
			return false;
		}

		betas = (double *) mxGetData(mx_betas);


		/* get image data */
		dim_image = mxGetDimensions( mx_image );
		if (nx != dim_image[0] || ny != dim_image[1] ||  !mxIsSingle(mx_image) ){
			mexErrMsgTxt("Error: mx_image must have [nx ny] float array.\n");
			return false;
		}

		image = mxGetPr_cfloat(mx_image);

		/* create memory space for proj */
		dim_proj[0] = nu;
		dim_proj[1] = noviews;

		plhs[0] = mxCreateNumericArray( 2 , dim_proj, mxSINGLE_CLASS, mxREAL);
		proj = ((float *) mxGetData(plhs[0]));


		if ( geom_type == PROJ_RD ) {
			sof = fbct_proj_rd(  nx, ny, dx, dy, offset_x, offset_y, nu, du, offset_u, sad, add, sdd, noviews, betas, detector_type, image, proj );
		}else if (geom_type == PROJ_DD ) {
			sof = fbct_proj_dd(  nx, ny, dx, dy, offset_x, offset_y, nu, du, offset_u, sad, add, sdd, noviews, betas, detector_type, image, proj );
		} else if (geom_type == PROJ_TF ) {
			sof = fbct_proj_tf(  nx, ny, dx, dy, offset_x, offset_y, nu, du, offset_u, sad, add, sdd, noviews, betas, detector_type, image, proj );
		}

		return sof;
}

/*
* static bool fbct_back_all( ... )
* back projection call function
*/
static bool fbct_back_all
	( mxArray *plhs[],
	Cmx mx_nxyz,
	Cmx mx_dxyz,
	Cmx mx_cxyz,
	Cmx mx_nuv,
	Cmx mx_duv,
	Cmx mx_cuv,
	Cmx mx_sad,
	Cmx mx_add,
	Cmx mx_sdd,
	Cmx mx_noviews,
	Cmx mx_betas,
	Cmx mx_detector_type,
	Cmx mx_proj,
	GEOM_TYPE geom_type){

#if Debug1
		mexPrintf("fbct_back_all()... \n");
#endif

		bool sof = false;

		int nx = mxGetPr_cint32(mx_nxyz)[0];
		int ny = mxGetPr_cint32(mx_nxyz)[1];
		float dx = mxGetPr_cfloat(mx_dxyz)[0];
		float dy = mxGetPr_cfloat(mx_dxyz)[1];
		float offset_x = mxGetPr_cfloat(mx_cxyz)[0];
		float offset_y = mxGetPr_cfloat(mx_cxyz)[1];

		int nu = mxGetPr_cint32(mx_nuv)[0];
		float du = mxGetPr_cfloat(mx_duv)[0];
		float offset_u = mxGetPr_cfloat(mx_cuv)[0];

		float sad = mxGetPr_cfloat(mx_sad)[0];
		float add = mxGetPr_cfloat(mx_add)[0];
		float sdd = mxGetPr_cfloat(mx_sdd)[0];
		int noviews = mxGetPr_cint32(mx_noviews)[0];
		int detector_type = mxGetPr_cint32(mx_detector_type)[0];

		double *betas;
		float *image;
		float *proj;
		int *dim_proj;
		mwSize dim_image[2];

		/* get parameters for geometry */
		if ( noviews != mxGetM(mx_betas) * mxGetN(mx_betas) || !mxIsRealDouble(mx_betas)){
			mexErrMsgTxt("Error: mx_betas must have noviews X 1 float array.\n");
			return false;
		}

		betas = (double *) mxGetData(mx_betas);


		/* get image data */
		dim_proj = mxGetDimensions( mx_proj );
		if (nu != dim_proj[0] || noviews != dim_proj[1] ||  !mxIsSingle(mx_proj) ){
			mexErrMsgTxt("Error: mx_proj must have [nu noviews] float array.\n");
			return false;
		}

		proj = mxGetPr_cfloat(mx_proj);

		/* create memory space for proj */
		dim_image[0] = nx;
		dim_image[1] = ny;

		plhs[0] = mxCreateNumericArray( 2 , dim_image, mxSINGLE_CLASS, mxREAL);
		image = ((float *) mxGetData(plhs[0]));

		if ( geom_type == BACK_PD ) {
			sof = fbct_back_pd(  nx, ny, dx, dy, offset_x, offset_y, nu, du, offset_u, sad, add, sdd, noviews, betas, detector_type, image, proj );
		}else if (geom_type == BACK_DD ) {
			sof = fbct_back_dd(  nx, ny, dx, dy, offset_x, offset_y, nu, du, offset_u, sad, add, sdd, noviews, betas, detector_type, image, proj );
		} else if (geom_type == BACK_TF ) {
			sof = fbct_back_tf(  nx, ny, dx, dy, offset_x, offset_y, nu, du, offset_u, sad, add, sdd, noviews, betas, detector_type, image, proj );
		}

		return sof;
}




/*
* bool fbct_geom_mex( ... )
* decide the type of forward or back projection
*/
bool fbct_geom_mex (int nlhs, mxArray *plhs[], int nrhs,   mxArray *prhs[])
{



	bool sof = false;
	char *arg;
	GEOM_TYPE geom_type;

	if ( nlhs <= 1 && !mxIsChar(prhs[0])) {
		mexErrMsgTxt("The 1st argument is invalid! \n");
		fbct_geom_mex_help();
		return false;
	}

	arg = mxu_getstring(prhs[0]);

	if (nrhs != 14 || nlhs != 1){
		mexErrMsgTxt("Error: incorrect number of arguments. \n");
		fbct_geom_mex_help();
		return false;
	}

	if ( !strncmp( arg, "proj", 4 ) ){

		if ( !strcmp( arg, "proj,rd" ) )
			geom_type = PROJ_RD;
		else if ( !strcmp( arg, "proj,dd" ) )
			geom_type = PROJ_DD;
		else if ( !strcmp( arg, "proj,tf" ) )
			geom_type = PROJ_TF;
		else{
			mexErrMsgTxt("Error: unknown geom type. \n");
			return false;
		}

		sof = fbct_proj_all(plhs, prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], \
			prhs[6], prhs[7], prhs[8], prhs[9], prhs[10], \
			prhs[11], prhs[12], prhs[13], geom_type );
	}

	if ( strncmp( arg, "back", 4 ) == 0 ){

		if ( !strcmp( arg, "back,pd" ) )
			geom_type = BACK_PD;
		else if ( !strcmp( arg, "back,dd" ) )
			geom_type = BACK_DD;
		else if ( !strcmp( arg, "back,tf" ) )
			geom_type = BACK_TF;
		else{
			mexErrMsgTxt("Error: unknown geom type. \n");
			return false;
		}

		sof = fbct_back_all(plhs, prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], \
			prhs[6], prhs[7], prhs[8], prhs[9], prhs[10], \
			prhs[11], prhs[12], prhs[13], geom_type );
	}


#if Debug1
	mexPrintf("fbct_geom_mex() successed! \n");
#endif

	mxu_freestring(arg);
	return sof;

}

/*
* mex entry function
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,   mxArray *prhs[]){
	if (!nlhs && !nrhs) {
		fbct_geom_mex_help();
		return;
	}
	if (!fbct_geom_mex(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("Error: fbct_geom_mex() failed. \n");
	return;

}



