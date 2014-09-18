#ifndef PARAMETERS
#define PARAMETERS

#include "stdafx.h"

//  size of simulation "m x m"
const float SIZE_X = 1.0, SIZE_Y = 1.0, SIZE_Z = 1.0;

//  support radius
const float h = 0.0457;

//  linked cell grid size
const int GRID_X = (int) (ceil((SIZE_X / h) + 0.1) + 0.5) + 1,
	GRID_Y = (int) (ceil((SIZE_Y / h) + 0.1) + 0.5) + 1,
	GRID_Z = (int) (ceil((SIZE_Z / h) + 0.1) + 0.5) + 1,
	GRIDSIZE = GRID_X * GRID_Y * GRID_Z;

//  gravity constant
const float g = -9.8;

//  time step
const float dt = 0.001;

//  rest density
const float p0 = 1000.0;

//  particle mass
const float mass = 0.02; // opt - N2

//  viscosity 
const float viscosity_const = 3.5;

//  surface tension 
const float st = 0.073;

//  threshold
const float l = 7.065;

//  gass stiffness
const float k = 3.5;

//  restitution
const float cr = 0.0;

// specification of initial particle grid size
const int init_slen_y = 13, init_slen_x = 13, init_slen_z = 13;
// make sure top values are put as floats
const float init_spacing_x = .5 / init_slen_x, 
	init_spacing_y = .5 / init_slen_y,
	init_spacing_z = .5 / init_slen_z;

const float c0 = 5.0;

const float c1 = 0.0;
const float c2 = 0.0;

const float gamma = 7.0;

const float bulk_modulus = p0 * c0 * c0 / gamma;

const float centre_x = 0.5, centre_y = 0.5, centre_z = 0.5;
const float radius = 0.3;

const float corr = 0.01 * h * h;

const float ggradius = 0.4 * h;

const float gglen = ggradius / sqrt(2.0);

const int ggx = (int) (ceil((SIZE_X / gglen) + 0.1) + 0.5), 
	ggy = (int) (ceil((SIZE_Y / gglen) + 0.1) + 0.5), 
	ggz = (int) (ceil((SIZE_Z / gglen) + 0.1) + 0.5), 
	ggsize = ggy * ggx * ggz;

const float thickness = 3.0;

const float gge = 1.085;

const float PI = 4.0*atan(1.0);
const float h6_term = 45 / PI / pow(h, 6);
const float h9_term = 945 / 32 / pow(h, 9);

#endif