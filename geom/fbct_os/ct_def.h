#include "math.h"


#ifndef _CT_DEF_H_
#define _CT_DEF_H_


#define sqr(x) (x * x)

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

typedef enum { PROJ_RD, PROJ_DD, PROJ_TF, BACK_PD, BACK_DD, BACK_TF, BACK_VAR }  GEOM_TYPE ;


#endif

