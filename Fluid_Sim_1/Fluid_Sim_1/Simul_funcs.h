#ifndef SIMUL_FUNCS
#define SIMUL_FUNCS

struct Cell;
struct Particle;
struct Container;

#include "stdafx.h"

void update_cellgrid(std::vector<Cell>& c_vec);
void reset_massd_and_forces(std::vector<Particle>& p_vec);
void update_mass_densities(std::vector<Cell>& c_vec);
void update_pressure_field(std::vector<Particle>& p_vec);
void update_forces(std::vector<Cell>& c_vec);
void update_velocities(std::vector<Particle>& p_vec, Container& cont);
void update_positions(std::vector<Particle>& p_vec);
void dumb_collisions(std::vector<Particle>& p_vec, Container& cont);

#endif