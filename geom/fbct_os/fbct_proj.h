#include "ct_def.h"

#ifndef _FBCT_PROJ_H_
#define _FBCT_PROJ_H_

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
	float *proj );



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
	float *proj );


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
	float *proj );


#endif

