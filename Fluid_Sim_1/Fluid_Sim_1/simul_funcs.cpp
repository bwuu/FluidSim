#include "stdafx.h"
#include "Parameters.h"
#include "SPH_Structs.h"

// check that gridsize and particle count are used properly 1

//  Particles are passed using their pointers, vectors are passed
//  by reference (&)

//  ************************************************************
//  Linked Cell Updater
//  ************************************************************

//  nc 1
//  updates all cells of c_vec
void update_cellgrid(std::vector<Cell>& c_vec) 
{
	//  iterate over each cell i in c_vec
	for (int i = 0; i < GRIDSIZE; i++)
	{
		//c_vec[i].surface_flag = false;
		std::list<Particle*>::iterator it = c_vec[i].ref_list.begin();

		//  iterate over particle*s of cell i
		while (it != c_vec[i].ref_list.end()) 
		{
			int correct_pos = ( (int) (floor((*it)->x[0] / h) + 0.5) )
				+ ( (int) (floor((*it)->x[1] / h) + 0.5) ) * GRID_X
				+ ( (int) (floor((*it)->x[2] / h) + 0.5) ) * GRID_Y * GRID_X;
			assert(correct_pos >= 0);
			assert(correct_pos < GRIDSIZE);

			//  erasing operation already handles it++
			if ( i != correct_pos ) 
			{
				c_vec[correct_pos].ref_list.push_back(*it);
				it = c_vec[i].ref_list.erase(it);
			}
			else it++;
		}
	}
	return;
}

//  ************************************************************
//  Mass-Density Updater
//  ************************************************************
//  calculate and update mutual mass-density(md) contributions
//  between particle *p and all particles in acell
void iterate_md_over_cell(Particle* p, Cell& acell)  
{
	std::list<Particle*>::iterator it = acell.ref_list.begin();

	// iterate over cell
	while (it != acell.ref_list.end()) 
	{
		const double sqr_dist = pow(p->x[0] - (*it)->x[0], 2) 
			+ pow(p->x[1] - (*it)->x[1], 2)
			+ pow(p->x[2] - (*it)->x[2], 2);

		//  break if dist > smoothing radius
		if (sqr_dist > h*h)
		{
			it++;
			continue;
		}

		//  can get huge
		const double kernel_res = 315 / (64 * PI * pow(h, 9)) * (pow((h*h - sqr_dist), 3));
		
		assert(kernel_res <= DBL_MAX && kernel_res >= -DBL_MAX);
		p->mass_density += kernel_res * mass;
		(*it)->mass_density += kernel_res * mass;

		it++;
	}
	return;
}


void reset_massd_and_forces(std::vector<Particle>& p_vec)
{
	for (int i = 0; i < p_vec.size(); i++) 
	{
		p_vec[i].mass_density = 0;

		// set gravity force
		p_vec[i].f_gravity[2] = g;

		//  reset other forces etc
		p_vec[i].f_pressure[0] = 0;
		p_vec[i].f_pressure[1] = 0;
		p_vec[i].f_pressure[2] = 0;
		p_vec[i].f_surface_tension[0] = 0;
		p_vec[i].f_surface_tension[1] = 0;
		p_vec[i].f_surface_tension[2] = 0;
		p_vec[i].f_viscosity[0] = 0;
		p_vec[i].f_viscosity[1] = 0;
		p_vec[i].f_viscosity[2] = 0;
		p_vec[i].inward_normal[0] = 0;
		p_vec[i].inward_normal[1] = 0;
		p_vec[i].inward_normal[2] = 0;
		p_vec[i].color_laplace = 0;
	}
}

//  update mass densities of all particles
void update_mass_densities(std::vector<Cell>& c_vec) 
{
	//  iterate over each cell i
	for (int i = 0; i < GRIDSIZE; i++)
	{
		//  iterate over cell contents using it1
		for (std::list<Particle*>::iterator it1 = c_vec[i].ref_list.begin(); 
				it1 != c_vec[i].ref_list.end(); it1++) 
		{
			//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! is this step correct/needed?
			//  do each particle with itself once  // self calc? -all
			(*it1)->mass_density += mass * 315 / (64 * PI * pow(h, 3));

			//  iterate over lexicographically greater cell contents using it2
			std::list<Particle*>::iterator it2 = it1;

			it2++;

			while (it2 != c_vec[i].ref_list.end()) 
			{ 
				double sqr_dist = pow((*it1)->x[0] - (*it2)->x[0], 2) 
					+ pow((*it1)->x[1] - (*it2)->x[1], 2)
					+ pow((*it1)->x[2] - (*it2)->x[2], 2);
				if (sqr_dist > h*h)
				{
					it2++;
					continue;
				}
				//  can get huge
				double kernel_res = 315 / (64 * PI * pow(h, 9)) * (pow((h*h - sqr_dist), 3)); //opt-minor -2

				(*it1)->mass_density += kernel_res * mass;
				(*it2)->mass_density += kernel_res * mass;
				it2++;
			}

			//  check all of these again !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			//  iterate over adjacent cells, up-left, up, up-right, right
			bool y_down = (i % (GRID_X * GRID_Y)) / GRID_X != 0;
			bool y_up = (i % (GRID_X * GRID_Y)) < GRID_X * (GRID_Y - 1);
			bool x_down = (i % GRID_X) != 0;
			bool x_up = (i % GRID_X) < (GRID_X - 1);
			bool z_up = i < GRID_X * GRID_Y * (GRID_Z - 1);
			if (y_up && x_down) 
			{
				iterate_md_over_cell(*it1, c_vec[i + GRID_X - 1]);
				if (z_up) iterate_md_over_cell(*it1, c_vec[i + GRID_X * GRID_Y + GRID_X - 1]);
			}
			if (y_up) 
			{
				iterate_md_over_cell(*it1, c_vec[i + GRID_X]);
				if (z_up) iterate_md_over_cell(*it1, c_vec[i + GRID_X * (GRID_Y + 1)]);
			}
			if (y_up && x_up) 
			{
				iterate_md_over_cell(*it1, c_vec[i + GRID_X + 1]);
				if (z_up) iterate_md_over_cell(*it1, c_vec[i + GRID_X * (GRID_Y + 1) + 1]);
			}
			if (x_up)
			{
				iterate_md_over_cell(*it1, c_vec[i + 1]);
				if (z_up) iterate_md_over_cell(*it1, c_vec[i + GRID_X * GRID_Y + 1]);
			}
			if (x_up && y_down && z_up) iterate_md_over_cell(*it1, c_vec[i + 1 + GRID_X * (GRID_Y - 1)]);
			if (y_down && z_up) iterate_md_over_cell(*it1, c_vec[i + GRID_X * (GRID_Y - 1)]);
			if (x_down && y_down && z_up) iterate_md_over_cell(*it1, c_vec[i - 1 + GRID_X * (GRID_Y - 1)]);
			if (x_down && z_up) iterate_md_over_cell(*it1, c_vec[i - 1 + GRID_X * GRID_Y]);

			if ((*it1)->mass_density < 1)
				(*it1)->mass_density = p0;

			assert((*it1)->mass_density <= DBL_MAX && (*it1)->mass_density >= -DBL_MAX);
		}
	}
	return;
}

//  ************************************************************
//  Pressure updater
//  ************************************************************
//  nc 1
//  updates pressures of all particles
void update_pressure_field(std::vector<Particle>& p_vec) 
{
	for (int i = 0; i < p_vec.size(); i++) 
	{
		// set pressure
		p_vec[i].pressure = (p_vec[i].mass_density >= p0) ?
			(p0 + bulk_modulus * (pow(p_vec[i].mass_density / p0, gamma) - 1)) : 0;

		//p_vec[i].pressure = k * (p_vec[i].mass_density - p0);
		assert(p_vec[i].pressure <= DBL_MAX && p_vec[i].pressure >= -DBL_MAX);
		assert(p_vec[i].pressure >= 0);
	}
	return;
}

//  ************************************************************
//  Force Updater - Pressure, Gravity, Viscosity, Surface Tension
//  ************************************************************

void iterate_forces_over_cell(Particle* p, Cell& acell)
{
	std::list<Particle*>::iterator it = acell.ref_list.begin();
	while (it != acell.ref_list.end()) 
	{ 
		const double x_disp = p->x[0] - (*it)->x[0];
		const double y_disp = p->x[1] - (*it)->x[1];
		const double z_disp = p->x[2] - (*it)->x[2];
		const double sqr_dist = x_disp * x_disp + y_disp * y_disp + z_disp * z_disp;

		if (sqr_dist > h*h)
		{
			it++;
			continue;
		}

		const double dist = sqrt(sqr_dist);
		const double inv_density_p = 1 / p->mass_density; /// wtf
		const double inv_density_it = 1 / (*it)->mass_density;
		double symmetric_calc_res;
		double x_force, y_force, z_force;
		
		//  pressure calc - force is actually acceleration lol 
		symmetric_calc_res = (p->pressure * pow(inv_density_p, 2) + (*it)->pressure * pow(inv_density_it, 2))
			* mass * h6_term / dist * pow(h - dist, 2);
		x_force = symmetric_calc_res * x_disp;
		y_force = symmetric_calc_res * y_disp;
		z_force = symmetric_calc_res * z_disp;
		p->f_pressure[0] += x_force;// / (*it)->mass_density;
		p->f_pressure[1] += y_force;// / (*it)->mass_density;
		p->f_pressure[2] += z_force;// / (*it)->mass_density;
		(*it)->f_pressure[0] -= x_force;// / p->mass_density;
		(*it)->f_pressure[1] -= y_force;// / p->mass_density;
		(*it)->f_pressure[2] -= z_force;// / p->mass_density;

		// viscosity calc
		symmetric_calc_res = viscosity_const * mass * h6_term * (h - dist);
		x_force = symmetric_calc_res * ((*it)->v[0] - p->v[0]);
		y_force = symmetric_calc_res * ((*it)->v[1] - p->v[1]);
		z_force = symmetric_calc_res * ((*it)->v[2] - p->v[2]);
		p->f_viscosity[0] += x_force * inv_density_it;
		p->f_viscosity[1] += y_force * inv_density_it;
		p->f_viscosity[2] += z_force * inv_density_it;
		(*it)->f_viscosity[0] -= x_force * inv_density_p;
		(*it)->f_viscosity[1] -= y_force * inv_density_p;
		(*it)->f_viscosity[2] -= z_force * inv_density_p;

		//  st calc  - not actually forces lol
		symmetric_calc_res = mass * -h9_term * (h*h - sqr_dist);
		x_force = symmetric_calc_res * (h*h-sqr_dist) * x_disp; //opt - also the dist components
		y_force = symmetric_calc_res * (h*h-sqr_dist) * y_disp;
		z_force = symmetric_calc_res * (h*h-sqr_dist) * z_disp;
		p->inward_normal[0] += x_force * inv_density_it;
		p->inward_normal[1] += y_force * inv_density_it;
		p->inward_normal[2] += z_force * inv_density_it;
		(*it)->inward_normal[0] -= x_force * inv_density_p;
		(*it)->inward_normal[1] -= y_force * inv_density_p;
		(*it)->inward_normal[2] -= z_force * inv_density_p;
		//  color field laplacian
		x_force = symmetric_calc_res * (3 * h * h - 7 * sqr_dist); 
		p->color_laplace += x_force * inv_density_it;
		(*it)->color_laplace += x_force * inv_density_p;

		it++;
	}
	return;
}

//  MAIN UPDATER
void update_forces(std::vector<Cell>& c_vec) 
{
	//  iterate over each cell i
	for (int i = 0; i < GRIDSIZE; i++) 
	{
		//  iterate over cell contents using it1
		for (std::list<Particle*>::iterator it1 = c_vec[i].ref_list.begin(); 
				it1 != c_vec[i].ref_list.end(); it1++) 
		{
			(*it1)->f_gravity[2] = g;

			//  iterate over lexicographically greater cell contents using it2
			std::list<Particle*>::iterator it2 = it1;
			it2++;

			//  iteration with it2
			while (it2 != c_vec[i].ref_list.end()) 
			{ 
				//  used by several subsequent calculations
				double sqr_dist = pow((*it1)->x[0] - (*it2)->x[0], 2) 
					+ pow((*it1)->x[1] - (*it2)->x[1], 2)
					+ pow((*it1)->x[2] - (*it2)->x[2], 2);

				if (sqr_dist > h*h)
				{
					it2++;
					continue;
				}

				double dist = sqrt(sqr_dist);
				double symmetric_calc_res;
				double x_force, y_force, z_force;
				double contrib;

				//  mutual pressure force calculation - negatives accounted for and cancelled
				/*
				double artdot = ((*it1)->x[0] - (*it2)->x[0]) * ((*it1)->v[0] - (*it2)->v[0])
					+ ((*it1)->x[1] - (*it2)->x[1]) * ((*it1)->v[1] - (*it2)->v[1])
					+ ((*it1)->x[2] - (*it2)->x[2]) * ((*it1)->v[2] - (*it2)->v[2]);
				double artific_kernel_res = 
					(artdot > 0) ? 0 : (artdot / (sqr_dist + corr));
					
				assert(artific_kernel_res <= DBL_MAX && artific_kernel_res >= -DBL_MAX);
				*/
				symmetric_calc_res = 
					h6_term * (h-dist) * (h-dist) * mass / dist
					* ((*it1)->pressure / pow((*it1)->mass_density, 2) 
					+ (*it2)->pressure / pow((*it2)->mass_density, 2));
					/*
					- c1 * c0 * 2 * h / ((*it1)->mass_density + (*it2)->mass_density) * artific_kernel_res
					+ c2 * 2 * h * h / ((*it1)->mass_density + (*it2)->mass_density) * artific_kernel_res * artific_kernel_res ); 
					*/
				/*
				(dist > h) ? 0 : (h6_term * (h-dist) * (h-dist) * mass / (dist + corr)
				* ( ((*it1)->pressure + (*it2)->pressure) / 2 ));
				*/
				/*
				- c1 * c0 * 2 * h / ((*it1)->mass_density + (*it2)->mass_density) * artific_kernel_res
				+ c2 * 2 * h / ((*it1)->mass_density + (*it2)->mass_density) * artific_kernel_res *
				artific_kernel_res ) );\*/
				

				assert(symmetric_calc_res <= DBL_MAX && symmetric_calc_res >= -DBL_MAX);
				
				x_force = symmetric_calc_res * ((*it1)->x[0] - (*it2)->x[0]); //opt
				y_force = symmetric_calc_res * ((*it1)->x[1] - (*it2)->x[1]);
				z_force = symmetric_calc_res * ((*it1)->x[2] - (*it2)->x[2]);
				assert(x_force <= DBL_MAX && x_force >= -DBL_MAX);
				assert(y_force <= DBL_MAX && y_force >= -DBL_MAX);
				assert(z_force <= DBL_MAX && z_force >= -DBL_MAX);
				(*it1)->f_pressure[0] += x_force;// / (*it2)->mass_density;
				(*it1)->f_pressure[1] += y_force;// / (*it2)->mass_density;
				(*it1)->f_pressure[2] += z_force;// / (*it2)->mass_density;
				(*it2)->f_pressure[0] -= x_force;// / (*it1)->mass_density;
				(*it2)->f_pressure[1] -= y_force;// / (*it1)->mass_density;
				(*it2)->f_pressure[2] -= z_force;// / (*it1)->mass_density;

				//  mutual viscosity force calculation
				symmetric_calc_res = 
					viscosity_const * h6_term * (h-dist) * mass;
				x_force = symmetric_calc_res * ((*it2)->v[0] - (*it1)->v[0]); 
				y_force = symmetric_calc_res * ((*it2)->v[1] - (*it1)->v[1]);
				z_force = symmetric_calc_res * ((*it2)->v[2] - (*it1)->v[2]);
				(*it1)->f_viscosity[0] += x_force / (*it2)->mass_density;
				(*it1)->f_viscosity[1] += y_force / (*it2)->mass_density;
				(*it1)->f_viscosity[2] += z_force / (*it2)->mass_density;
				(*it2)->f_viscosity[0] -= x_force / (*it1)->mass_density;
				(*it2)->f_viscosity[1] -= y_force / (*it1)->mass_density;
				(*it2)->f_viscosity[2] -= z_force / (*it1)->mass_density;


				//  surface tension properties calculat//  might need self calculationion
				symmetric_calc_res = 
					-h6_term * 21.0 / 32.0 / pow(h,3) * (h*h-sqr_dist) * mass;
				// inward normals
				x_force = symmetric_calc_res * (h*h-sqr_dist) * ((*it1)->x[0] - (*it2)->x[0]); //opt - also the dist components
				y_force = symmetric_calc_res * (h*h-sqr_dist) * ((*it1)->x[1] - (*it2)->x[1]);
				z_force = symmetric_calc_res * (h*h-sqr_dist) * ((*it1)->x[2] - (*it2)->x[2]);
				(*it1)->inward_normal[0] += x_force / (*it2)->mass_density;
				(*it1)->inward_normal[1] += y_force / (*it2)->mass_density;
				(*it1)->inward_normal[2] += z_force / (*it2)->mass_density;
				(*it2)->inward_normal[0] -= x_force / (*it1)->mass_density;
				(*it2)->inward_normal[1] -= y_force / (*it1)->mass_density;
				(*it2)->inward_normal[2] -= z_force / (*it1)->mass_density;
				//  color field laplacian
				contrib = symmetric_calc_res * (3 * h * h - 7 * sqr_dist); //opt - also the dist components
				(*it1)->color_laplace += contrib / (*it2)->mass_density;
				(*it2)->color_laplace += contrib / (*it1)->mass_density;

				it2++;

			}

			for (int c = 0; c < 3; c++) 
			{
				assert((*it1)->f_pressure[c] <= DBL_MAX && (*it1)->f_pressure[c] >= -DBL_MAX);
				assert((*it1)->f_viscosity[c] <= DBL_MAX && (*it1)->f_viscosity[c] >= -DBL_MAX);
			}

			bool y_down = (i % (GRID_X * GRID_Y)) / GRID_X != 0;
			bool y_up = (i % (GRID_X * GRID_Y)) < GRID_X * (GRID_Y - 1);
			bool x_down = (i % GRID_X) != 0;
			bool x_up = (i % GRID_X) < (GRID_X - 1);
			bool z_up = i < GRID_X * GRID_Y * (GRID_Z - 1);
			if (y_up && x_down) 
			{
				iterate_forces_over_cell(*it1, c_vec[i + GRID_X - 1]);
				if (z_up) 
				{
					iterate_forces_over_cell(*it1, c_vec[i + GRID_X * GRID_Y + GRID_X - 1]);
				}
			}
			if (y_up) 
			{
				iterate_forces_over_cell(*it1, c_vec[i + GRID_X]);
				if (z_up)
				{
					iterate_forces_over_cell(*it1, c_vec[i + GRID_X * (GRID_Y + 1)]);
				}
			}
			if (y_up && x_up) 
			{
				iterate_forces_over_cell(*it1, c_vec[i + GRID_X + 1]);
				if (z_up) 
				{
					iterate_forces_over_cell(*it1, c_vec[i + GRID_X * (GRID_Y + 1) + 1]);
				}
			}
			if (x_up)
			{
				iterate_forces_over_cell(*it1, c_vec[i + 1]);
				if (z_up) 
				{
					iterate_forces_over_cell(*it1, c_vec[i + GRID_X * GRID_Y + 1]);
				}
			}
			if (x_up && y_down && z_up) 
			{	
				iterate_forces_over_cell(*it1, c_vec[i + 1 + GRID_X * (GRID_Y - 1)]);
			}
			if (y_down && z_up) 
			{
				iterate_forces_over_cell(*it1, c_vec[i + GRID_X * (GRID_Y - 1)]);
			}
			if (x_down && y_down && z_up) 
			{
				iterate_forces_over_cell(*it1, c_vec[i - 1 + GRID_X * (GRID_Y - 1)]);
			}
			if (x_down && z_up) 
			{
				iterate_forces_over_cell(*it1, c_vec[i - 1 + GRID_X * GRID_Y]);
			}

			(*it1)->norm = sqrt(pow((*it1)->inward_normal[0],2) 
				+ pow((*it1)->inward_normal[1],2)
				+ pow((*it1)->inward_normal[2],2));

			if ((*it1)->norm > l) {
				(*it1)->f_surface_tension[0] = -st * (*it1)->color_laplace * (*it1)->inward_normal[0] / (*it1)->norm;
				(*it1)->f_surface_tension[1] = -st * (*it1)->color_laplace * (*it1)->inward_normal[1] / (*it1)->norm;
				(*it1)->f_surface_tension[2] = -st * (*it1)->color_laplace * (*it1)->inward_normal[2] / (*it1)->norm;
			}
			for (int c = 0; c < 3; c++) 
			{
				assert((*it1)->f_pressure[c] <= DBL_MAX && (*it1)->f_pressure[c] >= -DBL_MAX);
				assert((*it1)->f_viscosity[c] <= DBL_MAX && (*it1)->f_viscosity[c] >= -DBL_MAX);
			}
		}
		//opt

	}
	return;
}

//  ************************************************************
//  Velocity Updater
//  ************************************************************
// this DOES NOT IMPLEMENT LEAP FROG
void update_velocities(std::vector<Particle>& p_vec, Container& cont)
{
	for (int i = 0; i < p_vec.size(); i++) 
	{
		p_vec[i].v[0] += dt * (p_vec[i].f_gravity[0]
			+ p_vec[i].f_pressure[0] + (p_vec[i].f_viscosity[0] + p_vec[i].f_surface_tension[0]) / p_vec[i].mass_density
			- cont.v[0]);
		p_vec[i].v[1] += dt * (p_vec[i].f_gravity[1]
			+ p_vec[i].f_pressure[1] + (p_vec[i].f_viscosity[1] + p_vec[i].f_surface_tension[1]) / p_vec[i].mass_density
			- cont.v[1]);
		p_vec[i].v[2] += dt * (p_vec[i].f_gravity[2]
			+ p_vec[i].f_pressure[2] + (p_vec[i].f_viscosity[2] + p_vec[i].f_surface_tension[2]) / p_vec[i].mass_density
			- cont.v[2]);
		/*
		if (i == 5) std::cout << "pressures" << p_vec[i].f_pressure[0] << p_vec[i].f_pressure[1];
		if (i == 5) std::cout << "gravity" << p_vec[i].f_gravity[0] << p_vec[i].f_gravity[1];
		if (i == 5) std::cout << "visc" << p_vec[i].f_viscosity[0] << p_vec[i].f_viscosity[1];
		if (i == 5) std::cout << "st" << p_vec[i].f_surface_tension[0] << p_vec[i].f_surface_tension[1];
		if (i == 5) std::cout << "velos" << p_vec[i].v[0] << p_vec[i].v[1];
		if (i == 5) std::cout << "field" << p_vec[i].pressure;
		*/
	}
	return;
}

//  ************************************************************
//  Position Updater
//  ************************************************************

void update_positions(std::vector<Particle>& p_vec) {
	for (int i = 0; i < p_vec.size(); i++) {
		p_vec[i].x[0] += p_vec[i].v[0] * dt;
		p_vec[i].x[1] += p_vec[i].v[1] * dt;
		p_vec[i].x[2] += p_vec[i].v[2] * dt;
	}
}

double sgn(double val) {
    return (0 < val) - (val < 0);
}

void dumb_collisions(std::vector<Particle>& p_vec, Container& cont) 
{
	for (int i = 0; i < p_vec.size(); i++) 
	{
		double xlocal = p_vec[i].x[0] - centre_x;
		double ylocal = p_vec[i].x[1] - centre_y;
		double zlocal = p_vec[i].x[2] - centre_z;
		double ax = abs(xlocal), 
			ay = abs(ylocal),
			az = abs(zlocal);
		while (std::max(std::max(ax, ay), az) - radius > 0.0) 
		{ 
			double cpx = centre_x + 0.999 * std::min(radius, std::max(-radius, xlocal));
			double cpy = centre_y + 0.999 * std::min(radius, std::max(-radius, ylocal));
			double cpz = centre_z + 0.999 * std::min(radius, std::max(-radius, zlocal));
			double normalx = 0.0,
				normaly = 0.0,
				normalz = 0.0;
			if (ax > radius) normalx = sgn(cpx - centre_x);
			if (ay > radius) normaly = sgn(cpy - centre_y);
			if (az > radius) normalz = sgn(cpz - centre_z);
			double normnorm = sqrt(normalx * normalx + normaly * normaly + normalz * normalz);
			normalx = normalx / normnorm;
			normaly = normaly / normnorm;
			normalz = normalz / normnorm;

			assert(abs(p_vec[i].x[0] / p_vec[i].v[0]) > 0.0 ||
				!(std::cerr << p_vec[i].x[0] << ":" << p_vec[i].v[0] << ":" << abs(p_vec[i].x[0] / p_vec[i].v[0])));
			
			double dot = p_vec[i].v[0] * normalx 
				+ p_vec[i].v[1] * normaly
				+ p_vec[i].v[2] * normalz;
			double dcoeff = 1 + cr;
			p_vec[i].v[0] -= dcoeff * dot * normalx;
			p_vec[i].v[1] -= dcoeff * dot * normaly;
			p_vec[i].v[2] -= dcoeff * dot * normalz;
			
			p_vec[i].x[0] = cpx + 0.5 * dt * p_vec[i].v[0];
			p_vec[i].x[1] = cpy + 0.5 * dt * p_vec[i].v[1];
			p_vec[i].x[2] = cpz + 0.5 * dt * p_vec[i].v[2];
			
			xlocal = p_vec[i].x[0] - centre_x;
			ylocal = p_vec[i].x[1] - centre_y;
			zlocal = p_vec[i].x[2] - centre_z;
			ax = abs(xlocal);
			ay = abs(ylocal);
			az = abs(zlocal);

		}
		/*
		double indicate = pow(p_vec[i].x[0] - centre_x, 2) 
			+ pow(p_vec[i].x[1] - centre_y, 2)
			+ pow(p_vec[i].x[2] - centre_z, 2) - radius * radius;
		if (indicate > 0.0) 
		{
			double dist = sqrt(pow(p_vec[i].x[0] - centre_x, 2) 
				+ pow(p_vec[i].x[1] - centre_y, 2)
				+ pow(p_vec[i].x[2] - centre_z, 2));
			double cp[3] = { centre_x + 0.998 * radius * (p_vec[i].x[0] - centre_x) / dist, 
				centre_y + 0.999 * radius * (p_vec[i].x[1] - centre_y) / dist,
				centre_z + 0.999 * radius * (p_vec[i].x[2] - centre_z) / dist };
			double depth = abs(dist - radius); 
			double normal[3] = { (cp[0] - centre_x) / radius, 
				(cp[1] - centre_y) / radius,
				(cp[2] - centre_z) / radius };
			p_vec[i].x[0] = cp[0];
			p_vec[i].x[1] = cp[1];
			p_vec[i].x[2] = cp[2];

			double dot = p_vec[i].v[0] * normal[0] 
				+ p_vec[i].v[1] * normal[1]
				+ p_vec[i].v[2] * normal[2];
			double dcoeff = 1 + cr;//* depth / dt / sqrt(pow(p_vec[i].v[0], 2) + pow(p_vec[i].v[1], 2));
			p_vec[i].v[0] = p_vec[i].v[0] - dcoeff * dot * normal[0];
			p_vec[i].v[1] = p_vec[i].v[1] - dcoeff * dot * normal[1];
			p_vec[i].v[2] = p_vec[i].v[2] - dcoeff * dot * normal[2];
		}
		*/
		
	}
}
