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


#define mxGetPr_cint32(mx) ((  int *) mxGetPr(mx))
#define mxGetPr_cfloat(mx) ((  float *) mxGetPr(mx))
#define mxGetPr_cdouble(mx) ((cdouble *) mxGetPr(mx))

#define mxGetInt(mx) (*((  int *) mxGetData(mx)))
#define mxGetDouble(mx) (*((cdouble *) mxGetData(mx)))
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


typedef   mxArray *Cmx;
typedef enum { PROJ_RD, PROJ_DD, PROJ_TF, BACK_PD, BACK_DD, BACK_TF, BACK_VAR }  GEOM_TYPE ;



#endif

