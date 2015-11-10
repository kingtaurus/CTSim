#include "stdio.h"
#include "string.h"
#include "mex.h"
#include "math.h"
#include "ct_def.h"
#include "book.h"


#include "cbct_back_pd.c"

#define fdk_gpu TRUE

static void cbct_geom_mex_help(){
    enum_gpu();
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
    const int *dim_proj = mxGetDimensions( mx_proj );
    const int *dim_map = mxGetDimensions( mx_map );
    
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
    
#if fdk_gpu
    if ( geom_type == BACK_PD ) {
        sof = cbct_back_pd_gpu(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, map, true );
    } else if ( geom_type == BACK_HPD ) {
        sof = cbct_back_pd_gpu(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, map, false );
    }
#else
    if ( geom_type == BACK_PD ) {
        sof = cbct_back_pd_gpu(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, map, true );
    } else if ( geom_type == BACK_HPD ) {
        sof = cbct_back_pd_gpu(  nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, couchz, detector_type, image, proj, map, false );
    }
#endif

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
         mexErrMsgTxt("Not for foward projection. \n");
    }
    
    if ( strncmp( arg, "back", 4 ) == 0 ){
        
        if ( !strcmp( arg, "back,pd" ) )
            geom_type = BACK_PD;
        else if ( !strcmp( arg, "back,hpd" ) )
            geom_type = BACK_HPD;
        else{
            mexErrMsgTxt("Error: unknown geom type. \n");
            return false;
        }
        
        sof = cbct_back_all(plhs, prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], \
                    prhs[6], prhs[7], prhs[8], prhs[9], prhs[10], \
                    prhs[11], prhs[12], prhs[13], prhs[14], prhs[15], geom_type );
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

