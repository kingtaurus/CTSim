#include "ct_def.h"

#ifndef _FBCT_BACK_H_
#define _FBCT_BACK_H_

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
	float *proj );

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
	float *proj  );


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
	float *proj  );

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
	float *proj );


#endif

