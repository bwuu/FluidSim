#ifndef PARAMETERS
#define PARAMETERS

#include "stdafx.h"

//  size of simulation "m x m"
const double SIZE_X = 1.0, SIZE_Y = 1.0, SIZE_Z = 1.0;

//  support radius
const double h = 0.0457;

//  linked cell grid size
const int GRID_X = (int) (ceil((SIZE_X / h) + 0.1) + 0.5) + 1,
	GRID_Y = (int) (ceil((SIZE_Y / h) + 0.1) + 0.5) + 1,
	GRID_Z = (int) (ceil((SIZE_Z / h) + 0.1) + 0.5) + 1,
	GRIDSIZE = GRID_X * GRID_Y * GRID_Z;

//  gravity constant
const double g = -9.8;

//  time step
const double dt = 0.001;

//  rest density
const double p0 = 1000.0;

//  particle mass
const double mass = 0.02; // opt - N2

//  viscosity 
const double viscosity_const = 3.5;

//  surface tension 
const double st = 0.073;

//  threshold
const double l = 7.065;

//  gass stiffness
const double k = 3.5;

//  restitution
const double cr = 0.0;

// specification of initial particle grid size
const int init_slen_y = 13, init_slen_x = 13, init_slen_z = 13;
// make sure top values are put as doubles
const double init_spacing_x = .5 / init_slen_x, 
	init_spacing_y = .5 / init_slen_y,
	init_spacing_z = .5 / init_slen_z;

const double c0 = 5.0;

const double c1 = 0.0;
const double c2 = 0.0;

const double gamma = 7.0;

const double bulk_modulus = p0 * c0 * c0 / gamma;

const double centre_x = 0.5, centre_y = 0.5, centre_z = 0.5;
const double radius = 0.3;

const double corr = 0.01 * h * h;

const double ggradius = 0.4 * h;

const double gglen = ggradius / sqrt(2.0);

const int ggx = (int) (ceil((SIZE_X / gglen) + 0.1) + 0.5), 
	ggy = (int) (ceil((SIZE_Y / gglen) + 0.1) + 0.5), 
	ggz = (int) (ceil((SIZE_Z / gglen) + 0.1) + 0.5), 
	ggsize = ggy * ggx * ggz;

const double thickness = 3.0;

const double gge = 1.085;

const double PI = 4.0*atan(1.0);
const double h6_term = 45 / PI / pow(h, 6);
const double h9_term = 945 / 32 / pow(h, 9);

#endif