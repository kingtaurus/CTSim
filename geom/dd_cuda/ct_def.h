#include "math.h"

#ifndef _CT_DEF_H_
#define _CT_DEF_H_

#define Debug1 0
#define Debug3 0

#define sqr(x) (x * x)

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

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

typedef  const mxArray *Cmx;
typedef enum { PROJ_RD, PROJ_DD, PROJ_TF, BACK_PD, BACK_DD, BACK_TF, BACK_VAR, BACK_HPD }  GEOM_TYPE ;

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

static int iDivUp(int a, int b) {
	return (a % b != 0) ? (a / b + 1) : (a / b);
}

#endif

