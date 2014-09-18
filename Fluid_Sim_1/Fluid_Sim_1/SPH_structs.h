#ifndef SPH_STRUCTS
#define SPH_STRUCTS

#include "stdafx.h"
#include "Vec3.h"

struct Particle {
	float x[3];
	float v[3];
	float f_pressure[3];
	float f_viscosity[3];
	float f_gravity[3];
	float f_surface_tension[3];
	float inward_normal[3];
	float norm;
	float color_laplace;
	float pressure;
	float mass_density;
	bool   ghost_flag;

	Particle(float a, float b, float c) 
	{
		for (int i = 0; i < 3; i++) 
		{
			x[i] = 0;
			v[i] = 0;
			f_pressure[i] = 0;
			f_viscosity[i] = 0;
			f_gravity[i] = 0;
			f_surface_tension[i] = 0;
			inward_normal[i] = 0;
		}
		norm = 0;
		color_laplace = 0;
		pressure = 0;
		mass_density = 0;
		ghost_flag = false;
	}

	Particle() 
	{
		for (int i = 0; i < 3; i++) 
		{
			x[i] = 0;
			v[i] = 0;
			f_pressure[i] = 0;
			f_viscosity[i] = 0;
			f_gravity[i] = 0;
			f_surface_tension[i] = 0;
			inward_normal[i] = 0;
		}
		norm = 0;
		color_laplace = 0;
		pressure = 0;
		mass_density = 0;
		ghost_flag = false;
	}
};


struct ParticleRef {
	std::vector<Particle>* vecref;
	int index;

	Particle* getref() 
	{ 
		return &((*vecref)[index]); 
	}

};

struct Cell {
	// ***VECREFS NOT POINTERS, IS A SPECIAL CLASS
	std::list<Particle*> ref_list;
};

struct Container {
	float x[3];
	float v[3];

	Container() {
		for (int i = 0; i < 3; i++)
		{
			x[i] = 0;
			v[i] = 0;
		}
	}

	void move(float a, float b, float c) {
		x[0] += a;
		x[1] += b;
		x[2] += c;
	}
};
	
#endif