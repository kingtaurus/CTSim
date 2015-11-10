#include "stdio.h"
#include "string.h"
#include "mex.h"
#include "math.h"
#include "ct_def.h"


#include "cbct_back.c"


static void cbct_geom_mex_help(){
    mexPrintf("\n\
            coming soon!\
            \n");
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
        Cmx mx_proj){
    
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
    
    double *betas;
    float *image;
    float *proj;
    int *dim_proj;
    mwSize dim_image[3];
    
    /* get parameters for geometry */
    if ( noviews != mxGetM(mx_betas) * mxGetN(mx_betas) || !mxIsRealDouble(mx_betas)){
        mexErrMsgTxt("Error: mx_betas must have noviews X 1 float array.\n");
        return false;
    }
    
    betas = (double *) mxGetData(mx_betas);
    
    
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
    
    
    sof = cbct_back_dd(   nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, offset_u, offset_v, sad, add, sdd, noviews, betas, image, proj );
    
    return sof;
}




/*
 * bool fbct_geom_mex( ... )
 * decide the type of forward or back projection
 */
bool cbct_geom_mex (int nlhs, mxArray *plhs[], int nrhs,   mxArray *prhs[])
{
    

    bool sof = false;
    char *arg;
    
    if ( nlhs <= 1 && !mxIsChar(prhs[0])) {
        mexErrMsgTxt("The 1st argument is invalid! \n");
        cbct_geom_mex_help();
        return false;
    }
    
    arg = mxu_getstring(prhs[0]);
    
    if (nrhs != 13 || nlhs != 1){
        mexErrMsgTxt("Error: incorrect number of arguments. \n");
        cbct_geom_mex_help();
        return false;
    }
    
    sof = cbct_back_all(plhs, prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], \
            prhs[6], prhs[7], prhs[8], prhs[9], prhs[10], \
            prhs[11], prhs[12] );

#if Debug1
mexPrintf("cbct_geom_mex() successed! \n");
#endif

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



