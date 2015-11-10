#include "stdio.h"
#include "string.h"
#include "mex.h"
#include "math.h"
#include "ct_def.h"
#include "cbct_back_pmat.c"
#include "cbct_proj_pmat.c"




static void cbct_geom_mex_help(){
    mexPrintf("\n\
            coming soon!\
            \n");
}

static bool cbct_proj_all(
        mxArray *plhs[],
        Cmx mx_nxyz,
        Cmx mx_dxyz,
        Cmx mx_cxyz,
        Cmx mx_nuv,
        Cmx mx_duv,
        Cmx mx_noviews,
        Cmx mx_pmat,
        Cmx mx_cameras,
        Cmx mx_sdds,
        Cmx mx_cuvs,
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
    float du = mxGetPr_cfloat(mx_duv)[1];
    float dv = mxGetPr_cfloat(mx_duv)[2];
    
    int noviews = mxGetPr_cint32(mx_noviews)[0];
    
    float * pmat;
    float * cameras;
    float * sdds;
    float * offset_uv;
    float * image;
    float * proj;
    bool  * map;
    
    const unsigned long * dim_image = mxGetDimensions( mx_image );
    const unsigned long * dim_pmat = mxGetDimensions( mx_pmat );
    
    mwSize dim_proj[3];
    
    /* get image data */
    if ( nz != dim_image[2] || nx != dim_image[0] || ny != dim_image[1] ||  !mxIsSingle(mx_image) ){
        mexErrMsgTxt("Error: mx_image must have [nx ny nz] float array.\n");
        return false;
    }
    image = (float *) mxGetData(mx_image);
    
    /* get projection matrices */
    if ( 3 != dim_pmat[0] || 4 != dim_pmat[1] || noviews != dim_pmat[2] ||  !mxIsSingle(mx_pmat) ){
        mexErrMsgTxt("Error: mx_pmat must have [3 4 noViews ] float array.\n");
        return false;
    }
    pmat = (float *) mxGetData(mx_pmat);
    
    /* get cameras positions */
    if ( 3 * noviews != mxGetM(mx_cameras) * mxGetN(mx_cameras) || !mxIsSingle(mx_cameras)){
        mexErrMsgTxt("Error: mx_cameras must have noviews X 2 float array.\n");
        return false;
    }
    cameras = (float *) mxGetData(mx_cameras);
    
    /* get detector offsets */
    if ( 2 * noviews != mxGetM(mx_cuvs) * mxGetN(mx_cuvs) || !mxIsSingle(mx_cuvs)){
        mexErrMsgTxt("Error: mx_cuvs must have noviews X 2 float array.\n");
        return false;
    }
    offset_uv = (float *) mxGetData(mx_cuvs);
    
    
    /* get sdds */
    if ( noviews != mxGetM(mx_sdds) * mxGetN(mx_sdds) || !mxIsSingle(mx_sdds)){
        mexErrMsgTxt("Error: mx_sdds must have noviews X 1 float array.\n");
        return false;
    }
    sdds = (float *) mxGetData(mx_sdds);
    
    
    /* get image map */
    if ( nx * ny != mxGetM( mx_map ) * mxGetN( mx_map ) ||  !mxIsLogical(mx_map) ){
        mexErrMsgTxt("Error: mx_map must have [nx ny] bool array.\n");
        return false;
    }
    map = (bool *) mxGetData(mx_map);
    
    
    /* create memory space for proj */
    dim_proj[0] = (unsigned long) nv;
    dim_proj[1] = (unsigned long) nu;
    dim_proj[2] = (unsigned long) noviews;
    
    plhs[0] = mxCreateNumericArray( 3 , dim_proj, mxSINGLE_CLASS, mxREAL);
    proj = ((float *) mxGetData(plhs[0]));
    
    if ( geom_type == PROJ_RD ) {
        mexErrMsgTxt("Error: ray driven method is not supported. Please consider using distance driven method. \n");
        sof = false;
    }else if (geom_type == PROJ_DD ) {
        sof = cbct_proj_dd( nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, noviews, pmat, cameras, sdds, offset_uv, image, proj, map, true);
    } else if (geom_type == PROJ_TF ) {
        mexWarnMsgTxt("Warning: using distance driven method instead. \n");
        sof = cbct_proj_dd( nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, noviews, pmat, cameras, sdds, offset_uv, image, proj, map, true);
    }
    
    return sof;
}


static bool cbct_back_all(
        mxArray *plhs[],
        Cmx mx_nxyz,
        Cmx mx_dxyz,
        Cmx mx_cxyz,
        Cmx mx_nuv,
        Cmx mx_duv,
        Cmx mx_noviews,
        Cmx mx_pmat,
        Cmx mx_cameras,
        Cmx mx_sdds,
        Cmx mx_cuvs,
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
    float du = mxGetPr_cfloat(mx_duv)[1];
    float dv = mxGetPr_cfloat(mx_duv)[2];
    
    int noviews = mxGetPr_cint32(mx_noviews)[0];
    
    float * pmat;
    float * cameras;
    float * offset_uv;
    float * sdds;
    float * image;
    float * proj;
    bool * map;
    
    const unsigned long *dim_proj = mxGetDimensions( mx_proj );
    const unsigned long* dim_pmat = mxGetDimensions( mx_pmat );
    
    mwSize dim_image[3];
    
    /* get projection data */
    if (nv != dim_proj[0] ||nu != dim_proj[1] || noviews != dim_proj[2] ||  !mxIsSingle(mx_proj) ){
        mexErrMsgTxt("Error: mx_proj must have [nv nu noviews] float array.\n");
        return false;
    }
    proj = (float *) mxGetData(mx_proj);
    
    /* get projection matrices */
    if ( 3 != dim_pmat[0] || 4 != dim_pmat[1] || noviews != dim_pmat[2] ||  !mxIsSingle(mx_pmat) ){
        mexErrMsgTxt("Error: mx_pmat must have [3 4 noViews ] float array.\n");
        return false;
    }
    pmat = (float *) mxGetData(mx_pmat);
    
    /* get cameras positions */
    if ( 3 * noviews != mxGetM(mx_cameras) * mxGetN(mx_cameras) || !mxIsSingle(mx_cameras)){
        mexErrMsgTxt("Error: mx_cameras must have noviews X 2 float array.\n");
        return false;
    }
    cameras = (float *) mxGetData(mx_cameras);
    
    /* get detector offsets */
    if ( 2 * noviews != mxGetM(mx_cuvs) * mxGetN(mx_cuvs) || !mxIsSingle(mx_cuvs)){
        mexErrMsgTxt("Error: mx_cuvs must have noviews X 2 float array.\n");
        return false;
    }
    offset_uv = (float *) mxGetData(mx_cuvs);
    
    
    /* get sdds */
    if ( noviews != mxGetM(mx_sdds) * mxGetN(mx_sdds) || !mxIsSingle(mx_sdds)){
        mexErrMsgTxt("Error: mx_sdds must have noviews X 1 float array.\n");
        return false;
    }
    sdds = (float *) mxGetData(mx_sdds);
    
    if ( nx * ny != mxGetM( mx_map ) * mxGetN( mx_map ) ||  !mxIsLogical(mx_map) ){
        mexErrMsgTxt("Error: mx_map must have [nx ny] bool array.\n");
        return false;
    }
    map = (bool *) mxGetData(mx_map);
    
    /* create memory space for image */
    dim_image[0] = (unsigned long) nx;
    dim_image[1] = (unsigned long) ny;
    dim_image[2] = (unsigned long) nz;
    
    plhs[0] = mxCreateNumericArray( 3 , dim_image, mxSINGLE_CLASS, mxREAL);
    image = ((float *) mxGetData(plhs[0]));
    
    if ( geom_type == BACK_PD ) {
        sof = cbct_back_pd( nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, noviews, pmat, cameras, sdds, offset_uv, image, proj, map, true);
    } else if ( geom_type == BACK_HPD ) {
        sof = cbct_back_pd( nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, noviews, pmat, cameras, sdds, offset_uv, image, proj, map, false);
    } else if (geom_type == BACK_DD ) {
        sof = cbct_back_dd( nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, noviews, pmat, cameras, sdds, offset_uv, image, proj, map, false);
    } else if (geom_type == BACK_TF ) {
        mexWarnMsgTxt("Warning: using distance driven method instead. \n");
        sof = cbct_back_dd( nx, ny, nz, dx, dy, dz, offset_x, offset_y, offset_z, nu, nv, du, dv, noviews, pmat, cameras, sdds, offset_uv, image, proj, map, false);
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
        
        sof = cbct_proj_all(plhs, prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], \
                prhs[6], prhs[7], prhs[8], prhs[9], prhs[10], prhs[11], prhs[12], geom_type );
    }
    
    if ( !strncmp( arg, "back", 4 ) ){
        
        if ( !strcmp( arg, "back,pd" ) )
            geom_type = BACK_PD;
        else if ( !strcmp( arg, "back,hpd" ) )
            geom_type = BACK_HPD;
        else if ( !strcmp( arg, "back,dd" ) )
            geom_type = BACK_DD;
        else if ( !strcmp( arg, "back,tf" ) )
            geom_type = BACK_TF;
        else{
            mexErrMsgTxt("Error: unknown geom type. \n");
            return false;
        }
        
        sof = cbct_back_all(plhs, prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], \
                prhs[6], prhs[7], prhs[8], prhs[9], prhs[10], prhs[11], prhs[12], geom_type );
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

