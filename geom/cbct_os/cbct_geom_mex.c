#include "stdio.h"
#include "string.h"

typedef unsigned int char16_t;

#include "mex.h"
#include "math.h"

#include "ct_def.h"
#include "cbct_proj.h"
#include "cbct_back.h"

#define mxGetPr_cbool(mx) ((  bool *) mxGetPr(mx))
#define mxGetPr_cint32(mx) ((  int *) mxGetPr(mx))
#define mxGetPr_cfloat(mx) ((  float *) mxGetPr(mx))
#define mxGetPr_cdouble(mx) ((double *) mxGetPr(mx))

#define mxGetInt(mx) (*((  int *) mxGetData(mx)))
#define mxGetDouble(mx) (*((double *) mxGetData(mx)))
#define mxGetSingle(mx) (*((  float *) mxGetData(mx)))

#define mxIsScalar(mx) \
( (2 == mxGetNumberOfDimensions(mx)) \
        && (1 == mxGetM(mx)) && (1 == mxGetN(mx)) )
        
#define mxIsScalarInt32(mx) \
        ( mxIsScalar(mx) && mxIsInt32(mx) )
        
#define mxIsComplexSingle(mx) \
        (mxIsSingle(mx) && mxIsComplex(mx))
        
#define mxIsComplexDouble(mx) \
        (mxIsDouble(mx) && mxIsComplex(mx))
        
#define mxIsRealSingle(mx) \
        (mxIsSingle(mx) && !mxIsComplex(mx))
        
#define mxIsRealDouble(mx) \
        (mxIsDouble(mx) && !mxIsComplex(mx))
        
#define mxIsScalarSingle(mx) \
        ( mxIsScalar(mx) && mxIsRealSingle(mx) )
        
#define mxIsScalarDouble(mx) \
        ( mxIsScalar(mx) && mxIsRealDouble(mx) )
        
        typedef const mxArray *Cmx;


static void cbct_geom_mex_help(){
    mexPrintf("\n\
            coming soon!\
            \n");
}

/**************************************************************************
 * mxu_string()
 * caller must free using mxu_string_free()
 **************************************************************************/
char *mxu_getstring( Cmx mx)
{
    char	*string;
    unsigned long	n = mxGetM(mx) * mxGetN(mx) + 1;
    
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

/*
 * static bool fbct_proj_all( ... )
 * forward projection call function
 */
static bool cbct_proj_all
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
        Cmx mx_couchz,
        Cmx mx_detector_type,
        Cmx mx_image,
        GEOM_TYPE geom_type){
    
    bool sof = false;
    
    int nx = mxGetPr_cint32(mx_nxyz)[0];
    int ny = mxGetPr_cint32(mx_nxyz)[1];
    int nz = mxGetPr_cint32(mx_nxyz)[2];
    float dx = mxGetPr_cfloat(mx_dxyz)[0];
    float dy = mxGetPr_cfloat(mx_dxyz)[1];
    float dz = mxGetPr_cfloat(mx_dxyz)[2];
    float offset_x = mxGetPr_cfloat(mx_cxyz)[0];
    float offset_y = mxGetPr_cfloat(mx_cxyz)[1];
    float offset_z = mxGetPr_cfloat(mx_cxyz)[2];
    
    int nu = mxGetPr_cint32(mx_nuv)[0];
    int nv = mxGetPr_cint32(mx_nuv)[1];
    float du = mxGetPr_cfloat(mx_duv)[0];
    float dv = mxGetPr_cfloat(mx_duv)[1];
    float offset_u = mxGetPr_cfloat(mx_cuv)[0];
    float offset_v = mxGetPr_cfloat(mx_cuv)[1];
    
    float sad = mxGetPr_cfloat(mx_sad)[0];
    float add = mxGetPr_cfloat(mx_add)[0];
    float sdd = mxGetPr_cfloat(mx_sdd)[0];
    int noviews = mxGetPr_cint32(mx_noviews)[0];
    int detector_type = mxGetPr_cint32(mx_detector_type)[0];
    
    double* betas;
    float* couchz;
    float* image;
    float* proj;
    const unsigned long* dim_image = mxGetDimensions( mx_image );
    mwSize dim_proj[3];
    
    /* get parameters for projection angle */
    if (noviews != mxGetM(mx_betas) * mxGetN(mx_betas) || !mxIsRealDouble(mx_betas)){
        mexErrMsgTxt("Error: mx_betas must have noviews X 1 double array.\n");
        return false;
    }
    
    betas = (double *) mxGetData(mx_betas);
    
    /* get parameters for couch position */
    if (noviews != mxGetM(mx_couchz) * mxGetN(mx_couchz) || !mxIsSingle(mx_couchz)){
        mexErrMsgTxt("Error: mx_couchz must have noviews X 1 double array.\n");
        return false;
    }
    
    couchz = (float *) mxGetData(mx_couchz);
    
    
    /* get image data */
    if ( nz != dim_image[2] || nx != dim_image[0] || ny != dim_image[1] ||  !mxIsSingle(mx_image) ){
        mexErrMsgTxt("Error: mx_image must have [nx ny nz] float array.\n");
        return false;
    }
    
    image = mxGetPr_cfloat(mx_image);
    
    /* create memory space for proj */
    dim_proj[0] = (unsigned long) nv;
    dim_proj[1] = (unsigned long) nu;
    dim_proj[2] = (unsigned long) noviews;
    
    plhs[0] = mxCreateNumericArray( 3 , dim_proj, mxSINGLE_CLASS, mxREAL);
    proj = ((float *) mxGetData(plhs[0]));
    
    
    if ( geom_type == PROJ_RD ) {
        sof = cbct_proj_rd(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj );
    }else if (geom_type == PROJ_DD ) {
        sof = cbct_proj_dd(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj  );
    } else if (geom_type == PROJ_TF ) {
        sof = cbct_proj_tf(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj  );
    }
    
    return sof;
}


/*
* static bool fbct_proj_all( ... )
* forward projection call function
*/
static bool cbct_proj1_all
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
	Cmx mx_couchz,
	Cmx mx_detector_type,
	Cmx mx_image,
	Cmx mx_map,
	GEOM_TYPE geom_type){

		bool sof = false;

		int nx = mxGetPr_cint32(mx_nxyz)[0];
		int ny = mxGetPr_cint32(mx_nxyz)[1];
		int nz = mxGetPr_cint32(mx_nxyz)[2];
		float dx = mxGetPr_cfloat(mx_dxyz)[0];
		float dy = mxGetPr_cfloat(mx_dxyz)[1];
		float dz = mxGetPr_cfloat(mx_dxyz)[2];
		float offset_x = mxGetPr_cfloat(mx_cxyz)[0];
		float offset_y = mxGetPr_cfloat(mx_cxyz)[1];
		float offset_z = mxGetPr_cfloat(mx_cxyz)[2];

		int nu = mxGetPr_cint32(mx_nuv)[0];
		int nv = mxGetPr_cint32(mx_nuv)[1];
		float du = mxGetPr_cfloat(mx_duv)[0];
		float dv = mxGetPr_cfloat(mx_duv)[1];
		float offset_u = mxGetPr_cfloat(mx_cuv)[0];
		float offset_v = mxGetPr_cfloat(mx_cuv)[1];

		float sad = mxGetPr_cfloat(mx_sad)[0];
		float add = mxGetPr_cfloat(mx_add)[0];
		float sdd = mxGetPr_cfloat(mx_sdd)[0];
		int noviews = mxGetPr_cint32(mx_noviews)[0];
		int detector_type = mxGetPr_cint32(mx_detector_type)[0];

		double *betas;
		float *couchz;
		float *image;
		float *proj;
		bool *map;        
        const unsigned long* dim_image = mxGetDimensions( mx_image );
        const unsigned long* dim_map = mxGetDimensions( mx_map );
		mwSize dim_proj[3];

		/* get parameters for projection angle */
		if (noviews != mxGetM(mx_betas) * mxGetN(mx_betas) || !mxIsRealDouble(mx_betas)){
			mexErrMsgTxt("Error: mx_betas must have noviews X 1 double array.\n");
			return false;
		}

		betas = (double *) mxGetData(mx_betas);

		/* get parameters for couch position */
		if (noviews != mxGetM(mx_couchz) * mxGetN(mx_couchz) || !mxIsSingle(mx_couchz)){
			mexErrMsgTxt("Error: mx_couchz must have noviews X 1 double array.\n");
			return false;
		}

		couchz = (float *) mxGetData(mx_couchz);


		/* get image data */
		if ( nz != dim_image[2] || nx != dim_image[0] || ny != dim_image[1] ||  !mxIsSingle(mx_image) ){
			mexErrMsgTxt("Error: mx_image must have [nx ny nz] float array.\n");
			return false;
		}

		image = mxGetPr_cfloat(mx_image);


		/* get image map */
		if ( nx != dim_map[0] || ny != dim_map[1] ||  !mxIsLogical(mx_map) ){
			mexErrMsgTxt("Error: mx_map must have [nx ny] bool array.\n");
			return false;
		}

		map = mxGetPr_cbool(mx_map);


		/* create memory space for proj */
		dim_proj[0] = nv;
		dim_proj[1] = nu;
		dim_proj[2] = noviews;

		plhs[0] = mxCreateNumericArray( 3 , dim_proj, mxSINGLE_CLASS, mxREAL);
		proj = ((float *) mxGetData(plhs[0]));


		if (geom_type == PROJ_DD ) {
			sof = cbct_proj1_dd(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, map  );
		} else if (geom_type == PROJ_TF ) {
			sof = cbct_proj1_tf(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, map  );
		}

		return sof;
}



/*
 * static bool fbct_back_all( ... )
 * back projection call function
 */
static bool cbct_back_all
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
        Cmx mx_couchz,
        Cmx mx_detector_type,
        Cmx mx_proj,
        GEOM_TYPE geom_type){
    
    bool sof = false;
    
    int nx = mxGetPr_cint32(mx_nxyz)[0];
    int ny = mxGetPr_cint32(mx_nxyz)[1];
    int nz = mxGetPr_cint32(mx_nxyz)[2];
    float dx = mxGetPr_cfloat(mx_dxyz)[0];
    float dy = mxGetPr_cfloat(mx_dxyz)[1];
    float dz = mxGetPr_cfloat(mx_dxyz)[2];
    float offset_x = mxGetPr_cfloat(mx_cxyz)[0];
    float offset_y = mxGetPr_cfloat(mx_cxyz)[1];
    float offset_z = mxGetPr_cfloat(mx_cxyz)[2];
    
    int nu = mxGetPr_cint32(mx_nuv)[0];
    int nv = mxGetPr_cint32(mx_nuv)[1];
    float du = mxGetPr_cfloat(mx_duv)[0];
    float dv = mxGetPr_cfloat(mx_duv)[1];
    float offset_u = mxGetPr_cfloat(mx_cuv)[0];
    float offset_v = mxGetPr_cfloat(mx_cuv)[1];
    
    float sad = mxGetPr_cfloat(mx_sad)[0];
    float add = mxGetPr_cfloat(mx_add)[0];
    float sdd = mxGetPr_cfloat(mx_sdd)[0];
    int noviews = mxGetPr_cint32(mx_noviews)[0];
    int detector_type = mxGetPr_cint32(mx_detector_type)[0];
    
    double *betas;
    float *couchz;
    float *image;
    float *proj;
    const unsigned long *dim_proj = mxGetDimensions( mx_proj );
    mwSize dim_image[3];
    
    /* get parameters for geometry */
    if ( noviews != mxGetM(mx_betas) * mxGetN(mx_betas) || !mxIsRealDouble(mx_betas)){
        mexErrMsgTxt("Error: mx_betas must have noviews X 1 float array.\n");
        return false;
    }
    
    betas = (double *) mxGetData(mx_betas);
    
    /* get parameters for couch position */
    if (noviews != mxGetM(mx_couchz) * mxGetN(mx_couchz) || !mxIsSingle(mx_couchz)){
        mexErrMsgTxt("Error: mx_couchz must have noviews X 1 double array.\n");
        return false;
    }
    
    couchz = (float *) mxGetData(mx_couchz);
    
    
    /* get image data */
    if (nv != dim_proj[0] ||nu != dim_proj[1] || noviews != dim_proj[2] ||  !mxIsSingle(mx_proj) ){
        mexErrMsgTxt("Error: mx_proj must have [nv nu noviews] float array.\n");
        return false;
    }
    
    proj = mxGetPr_cfloat(mx_proj);
    
    /* create memory space for proj */
    dim_image[0] = (unsigned long) nx;
    dim_image[1] = (unsigned long) ny;
    dim_image[2] = (unsigned long) nz;
    
    plhs[0] = mxCreateNumericArray( 3 , dim_image, mxSINGLE_CLASS, mxREAL);
    image = ((float *) mxGetData(plhs[0]));
    
    if ( geom_type == BACK_PD ) {
        sof = cbct_back_pd(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, true );
    } else if ( geom_type == BACK_HPD ) {
        sof = cbct_back_pd(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, false );
    } else if (geom_type == BACK_DD ) {
        sof = cbct_back_dd(   nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj );
    } else if (geom_type == BACK_TF ) {
        sof = cbct_back_tf(   nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj );
    }
    
    return sof;
}




/*
* static bool fbct_back_all( ... )
* back projection call function
*/
static bool cbct_back1_all
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
	Cmx mx_couchz,
	Cmx mx_detector_type,
	Cmx mx_proj,
	Cmx mx_map,
	GEOM_TYPE geom_type){

		bool sof = false;

		int nx = mxGetPr_cint32(mx_nxyz)[0];
		int ny = mxGetPr_cint32(mx_nxyz)[1];
		int nz = mxGetPr_cint32(mx_nxyz)[2];
		float dx = mxGetPr_cfloat(mx_dxyz)[0];
		float dy = mxGetPr_cfloat(mx_dxyz)[1];
		float dz = mxGetPr_cfloat(mx_dxyz)[2];
		float offset_x = mxGetPr_cfloat(mx_cxyz)[0];
		float offset_y = mxGetPr_cfloat(mx_cxyz)[1];
		float offset_z = mxGetPr_cfloat(mx_cxyz)[2];

		int nu = mxGetPr_cint32(mx_nuv)[0];
		int nv = mxGetPr_cint32(mx_nuv)[1];
		float du = mxGetPr_cfloat(mx_duv)[0];
		float dv = mxGetPr_cfloat(mx_duv)[1];
		float offset_u = mxGetPr_cfloat(mx_cuv)[0];
		float offset_v = mxGetPr_cfloat(mx_cuv)[1];

		float sad = mxGetPr_cfloat(mx_sad)[0];
		float add = mxGetPr_cfloat(mx_add)[0];
		float sdd = mxGetPr_cfloat(mx_sdd)[0];
		int noviews = mxGetPr_cint32(mx_noviews)[0];
		int detector_type = mxGetPr_cint32(mx_detector_type)[0];

		double *betas;
		float *couchz;
		float *image;
		float *proj;
		bool *map;
        const unsigned long *dim_proj = mxGetDimensions( mx_proj );
        const unsigned long *dim_map = mxGetDimensions( mx_map );
        
		mwSize dim_image[3];

		/* get parameters for geometry */
		if ( noviews != mxGetM(mx_betas) * mxGetN(mx_betas) || !mxIsRealDouble(mx_betas)){
			mexErrMsgTxt("Error: mx_betas must have noviews X 1 float array.\n");
			return false;
		}

		betas = (double *) mxGetData(mx_betas);

		/* get parameters for couch position */
		if (noviews != mxGetM(mx_couchz) * mxGetN(mx_couchz) || !mxIsSingle(mx_couchz)){
			mexErrMsgTxt("Error: mx_couchz must have noviews X 1 double array.\n");
			return false;
		}

		couchz = (float *) mxGetData(mx_couchz);


		/* get image map */
		dim_map = mxGetDimensions( mx_map );
		if ( nx != dim_map[0] || ny != dim_map[1] ||  !mxIsLogical(mx_map) ){
			mexErrMsgTxt("Error: mx_map must have [nx ny] bool array.\n");
			return false;
		}

		map = mxGetPr_cbool(mx_map);


		/* get image data */
		dim_proj = mxGetDimensions( mx_proj );
		if (nv != dim_proj[0] ||nu != dim_proj[1] || noviews != dim_proj[2] ||  !mxIsSingle(mx_proj) ){
			mexErrMsgTxt("Error: mx_proj must have [nv nu noviews] float array.\n");
			return false;
		}

		proj = mxGetPr_cfloat(mx_proj);

		/* create memory space for proj */
		dim_image[0] = nx;
		dim_image[1] = ny;
		dim_image[2] = nz;

		plhs[0] = mxCreateNumericArray( 3 , dim_image, mxSINGLE_CLASS, mxREAL);
		image = ((float *) mxGetData(plhs[0]));

		if ( geom_type == BACK_PD ) {
			sof = cbct_back1_pd(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, map, true );
		} else if ( geom_type == BACK_HPD ) {
			sof = cbct_back1_pd(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, map, false );
		} else if (geom_type == BACK_DD ) {
			sof = cbct_back1_dd(   nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, map );
		} else if (geom_type == BACK_TF ) {
			sof = cbct_back1_tf(   nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, map );
		} 

		return sof;
}


/*
 * bool fbct_geom_mex( ... )
 * decide the type of forward or back projection
 */
bool cbct_geom_mex (int nlhs, mxArray *plhs[], int nrhs, const  mxArray *prhs[])
{   
    bool sof = false;
    char *arg;
    GEOM_TYPE geom_type;
    
    if ( nlhs <= 1 && !mxIsChar(prhs[0])) {
        mexErrMsgTxt("The 1st argument is invalid! \n");
        cbct_geom_mex_help();
        return false;
    }
    
    arg = mxu_getstring(prhs[0]);
    

    
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
        
        if( nrhs == 15 ){
            sof = cbct_proj_all(plhs, prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], \
                    prhs[6], prhs[7], prhs[8], prhs[9], prhs[10], \
                    prhs[11], prhs[12], prhs[13], prhs[14], geom_type );
        }else{
            sof = cbct_proj1_all(plhs, prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], \
                    prhs[6], prhs[7], prhs[8], prhs[9], prhs[10], \
                    prhs[11], prhs[12], prhs[13], prhs[14], prhs[15], geom_type );
        }
    }
    
    if ( strncmp( arg, "back", 4 ) == 0 ){
        
        if ( !strcmp( arg, "back,pd" ) )
            geom_type = BACK_PD;
        else if ( !strcmp( arg, "back,dd" ) )
            geom_type = BACK_DD;
        else if ( !strcmp( arg, "back,tf" ) )
            geom_type = BACK_TF;
        else if ( !strcmp( arg, "back,hpd" ) )
            geom_type = BACK_HPD;
        else{
            mexErrMsgTxt("Error: unknown geom type. \n");
            return false;
        }
        
        if( nrhs == 15 ){
            sof = cbct_back_all(plhs, prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], \
                    prhs[6], prhs[7], prhs[8], prhs[9], prhs[10], \
                    prhs[11], prhs[12], prhs[13], prhs[14], geom_type );
        }else{
            sof = cbct_back1_all(plhs, prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], \
                    prhs[6], prhs[7], prhs[8], prhs[9], prhs[10], \
                    prhs[11], prhs[12], prhs[13], prhs[14], prhs[15], geom_type );
        }
    }
    
    
    mxu_freestring(arg);
    return sof;
    
}

/*
 * mex entry function
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if (!nlhs && !nrhs) {
        cbct_geom_mex_help();
        return;
    }
    if (!cbct_geom_mex(nlhs, plhs, nrhs, prhs))
        mexErrMsgTxt("Error: cbct_geom_mex() failed. \n");
    return;
    
}



