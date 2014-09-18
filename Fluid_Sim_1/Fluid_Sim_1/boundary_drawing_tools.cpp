#include "Parameters.h"
#include "stdafx.h"
#include "SPH_structs.h"

// lean lower left
float ggposx(int i) 
{  
	return (i % ggx) * gglen;  
}

float ggposy(int i) 
{
	return ( (i % (ggx * ggy)) / ggx) * gglen;  
}

float ggposz(int i) 
{
	return (i / (ggx * ggy)) * gglen;  
}

void gg_rand_ingrid(int i, float ghostpos[])
{
	ghostpos[0] = ggposx(i) + ( gglen * rand() ) / RAND_MAX;
	ghostpos[1] = ggposy(i) + ( gglen * rand() ) / RAND_MAX;
	ghostpos[2] = ggposz(i) + ( gglen * rand() ) / RAND_MAX;
}

void gg_rand_ow_ext(float ghostpos[])
{
	float randvec[2] = { (rand() % 10001 - 5000) / 5000.0, (rand() % 10001 - 5000) / 5000.0 };
	float root = sqrt(pow(randvec[0], 2) + pow(randvec[1], 2) + pow(randvec[2], 2));
	randvec[0] /= root;
	randvec[1] /= root;

	float ax = abs(ghostpos[0] - centre_x);
	float ay = abs(ghostpos[1] - centre_y);
	float az = abs(ghostpos[2] - centre_z);
	if (ax >= ay && ax >= az)
	{
		ghostpos[1] += gglen * gge * randvec[0];
		ghostpos[2] += gglen * gge * randvec[1];
	}
	else if (ay >= ax && ay >= az)
	{
		ghostpos[0] += gglen * gge * randvec[0];
		ghostpos[2] += gglen * gge * randvec[1];
	}
	else if (az >= ay && az >= ax)
	{
		ghostpos[0] += gglen * gge * randvec[0];
		ghostpos[1] += gglen * gge * randvec[1];
	}
	else assert(false);// remove
}

void gg_rand_ow_int(float ghostpos[])
{
	float randvec[3] = { (rand() % 10001 - 5000) / 5000.0, (rand() % 10001 - 5000) / 5000.0 };
	float root = sqrt(pow(randvec[0], 2) + pow(randvec[1], 2) + pow(randvec[2], 2));
	randvec[0] /= root;
	randvec[1] /= root;
	randvec[2] /= root;
	assert(randvec[0] < 1.1);
	assert(randvec[1] < 1.1);
	assert(randvec[2] < 1.1);

	float randnum = (rand() % 1001) / 1000.0;
	ghostpos[0] += randvec[0] * (1 + randnum) * ggradius;
	ghostpos[1] += randvec[1] * (1 + randnum) * ggradius;
	ghostpos[2] += randvec[2] * (1 + randnum) * ggradius;
}

int getmax_argpos(float a, float b, float c)
{
	if (a > b && a > c) return 1;
	if (b > a && b > c) return 2;
	if (c > a && c > b) return 3;
	//   if they equal-write a better catch case later
	return 1;
}


void gg_project(float ghostpos[])
{
	float xlocal = ghostpos[0] - centre_x;
	float ylocal = ghostpos[1] - centre_y;
	float zlocal = ghostpos[2] - centre_z;
	float ax = abs(xlocal); 
	float ay = abs(ylocal);
	float az = abs(zlocal);
	//  ????  does it work when there equal -all
	bool x2inner = ax < radius + gglen * thickness / 2;
	bool y2inner = ay < radius + gglen * thickness / 2;
	bool z2inner = az < radius + gglen * thickness / 2;
	if (x2inner && y2inner && z2inner)
	{
		switch(getmax_argpos(ax, ay, az)) 
		{
		case 1:
			if (xlocal > 0) ghostpos[0] = centre_x + radius + gglen * thickness;
			else ghostpos[0] = centre_x - radius - gglen * thickness;
			break;
		case 2:
			if (ylocal > 0) ghostpos[1] = centre_y + radius + gglen * thickness;
			else ghostpos[1] = centre_y - radius - gglen * thickness;
			break;
		case 3:
			if (zlocal > 0) ghostpos[2] = centre_z + radius + gglen * thickness;
			else ghostpos[2] = centre_z - radius - gglen * thickness;
			break;
		default: 
			assert(false || !(std::cerr << "getmax_argpos returned invalid value"));
		}
		return;
	}
	if (!x2inner && y2inner && z2inner)
	{
		if (xlocal >= 0) ghostpos[0] = centre_x + radius + gglen * thickness;
		else ghostpos[0] = centre_x - radius - gglen * thickness;
		return;
	}
	if (x2inner && !y2inner && z2inner)
	{
		if (ylocal >= 0) ghostpos[1] = centre_y + radius + gglen * thickness;
		else ghostpos[1] = centre_y - radius - gglen * thickness;
		return;
	}
	if (x2inner && y2inner && !z2inner)
	{
		if (zlocal >= 0) ghostpos[2] = centre_z + radius + gglen * thickness;
		else ghostpos[2] = centre_z - radius - gglen * thickness;
		return;
	}
	if (!x2inner && !y2inner && z2inner)
	{
		if (xlocal >= 0) ghostpos[0] = centre_x + radius + gglen * thickness;
		else ghostpos[0] = centre_x - radius - gglen * thickness;
		if (ylocal >= 0) ghostpos[1] = centre_y + radius + gglen * thickness;
		else ghostpos[1] = centre_y - radius - gglen * thickness;
		return;
	}
	if (!x2inner && y2inner && !z2inner)
	{
		if (xlocal >= 0) ghostpos[0] = centre_x + radius + gglen * thickness;
		else ghostpos[0] = centre_x - radius - gglen * thickness;
		if (zlocal >= 0) ghostpos[2] = centre_z + radius + gglen * thickness;
		else ghostpos[2] = centre_z - radius - gglen * thickness;
		return;
	}
	if (x2inner && !y2inner && !z2inner)
	{
		if (ylocal >= 0) ghostpos[1] = centre_y + radius + gglen * thickness;
		else ghostpos[1] = centre_y - radius - gglen * thickness;
		if (zlocal >= 0) ghostpos[2] = centre_z + radius + gglen * thickness;
		else ghostpos[2] = centre_z - radius - gglen * thickness;
		return;
	}
	else
	{
		if (xlocal >= 0) ghostpos[0] = centre_x + radius + gglen * thickness;
		else ghostpos[0] = centre_x - radius - gglen * thickness;
		if (ylocal >= 0) ghostpos[1] = centre_y + radius + gglen * thickness;
		else ghostpos[1] = centre_y - radius - gglen * thickness;
		if (zlocal >= 0) ghostpos[2] = centre_z + radius + gglen * thickness;
		else ghostpos[2] = centre_z - radius - gglen * thickness;
		return;
	}
	assert(false);
}

void sample_ext(std::vector<Particle>& static_vec) 
{
	std::vector<int> active_grid_vec;
	//  get active grid
	int count = 0;
	for (int i = 0; i < ggsize; i++)
	{
		bool x_down_strip = (centre_x - radius > ggposx(i) && centre_x - radius < ggposx(i) + gglen)
			|| ( (centre_x - radius - thickness * gglen > ggposx(i)) && (centre_x - radius - thickness * gglen < ggposx(i) + gglen) );
		bool x_up_strip = (centre_x + radius > ggposx(i) && centre_x + radius < ggposx(i) + gglen)
			|| ( (centre_x + radius + thickness * gglen > ggposx(i)) && (centre_x + radius + thickness * gglen < ggposx(i) + gglen) );

		bool y_up_strip = (centre_y + radius > ggposy(i) && centre_y + radius < ggposy(i) + gglen)
			|| ( (centre_y + radius + thickness * gglen > ggposy(i)) && (centre_y + radius + thickness * gglen < ggposy(i) + gglen) );
		bool y_down_strip = (centre_y - radius > ggposy(i) && centre_y - radius < ggposy(i) + gglen)
			|| ( (centre_y - radius - thickness * gglen > ggposy(i)) && (centre_y - radius - thickness * gglen < ggposy(i) + gglen) );

		bool z_up_strip = (centre_z + radius > ggposz(i) && centre_z + radius < ggposz(i) + gglen)
			|| ( (centre_z + radius + thickness * gglen > ggposz(i)) && (centre_z + radius + thickness * gglen < ggposz(i) + gglen) );
		bool z_down_strip = (centre_z - radius > ggposz(i) && centre_z - radius < ggposz(i) + gglen)
			|| ( (centre_z - radius - thickness * gglen > ggposz(i)) && (centre_z - radius - thickness * gglen < ggposz(i) + gglen) );

		bool x_bounded = (ggposx(i) < centre_x + radius + thickness * gglen && ggposx(i) + gglen > centre_x - radius - thickness * gglen);
		bool y_bounded = (ggposy(i) < centre_y + radius + thickness * gglen && ggposy(i) + gglen > centre_y - radius - thickness * gglen);
		bool z_bounded = (ggposz(i) < centre_z + radius + thickness * gglen && ggposz(i) + gglen > centre_z - radius - thickness * gglen);

		if (   ( (x_down_strip || x_up_strip) && y_bounded && z_bounded)
			|| ( (y_down_strip || y_up_strip) && x_bounded && z_bounded)
			|| ( (z_down_strip || z_up_strip) && x_bounded && y_bounded) )
		{
			active_grid_vec.push_back(i);
			count++;
		}
	}
	std::cout << count;

	//  seed
	srand((unsigned)time(0));

	bool sample_success;
	for ( ; !active_grid_vec.empty(); active_grid_vec.pop_back() )
	{
		float ghost_pos[3];
		sample_success = false;
		for (int i = 0; i < 30; i++)
		{
			gg_rand_ingrid(active_grid_vec.back(), ghost_pos);
			gg_project(ghost_pos);
			std::vector<Particle>::iterator it = static_vec.begin();
			bool obstructed = false;
			while (it != static_vec.end())
			{
				if ( (pow(it->x[0] - ghost_pos[0], 2)
				    + pow(it->x[1] - ghost_pos[1], 2)  
					+ pow(it->x[2] - ghost_pos[2], 2)) < ggradius * ggradius)
				{ 
					obstructed = true;
					break;
				}
				it++;
			}
			if (obstructed) continue;
			else 
			{
				sample_success = true;
				break;
			}
		}

		while (sample_success)
		{
			Particle p (ghost_pos[0], ghost_pos[1], ghost_pos[2]);
			p.ghost_flag = true;
			static_vec.push_back(p);
			gg_rand_ow_ext(ghost_pos);
			gg_project(ghost_pos);

			std::vector<Particle>::iterator it = static_vec.begin();
			bool obstructed = false;
			while (it != static_vec.end())
			{
				if ( (pow(it->x[0] - ghost_pos[0], 2)
				    + pow(it->x[1] - ghost_pos[1], 2)
					+ pow(it->x[2] - ghost_pos[2], 2)) < ggradius * ggradius)
				{ 
					obstructed = true;
					break;
				}
				it++;
			}
			if (obstructed) sample_success = false;
		}
	}
}

//  im scared of move in memory vectors
void sample_int(std::vector<Particle>& static_vec)
{
	srand((unsigned)time(0));

	for (int i = 0; i < static_vec.size(); i++)
	{
		float base_ghostpos[3] = { static_vec[i].x[0], static_vec[i].x[1], static_vec[i].x[2] };
		float ghostpos[3];

		bool sample_success = false;
		for (int j = 0; j < 30; j++)
		{
			ghostpos[0] = base_ghostpos[0]; 
			ghostpos[1] = base_ghostpos[1];
			ghostpos[2] = base_ghostpos[2];
			gg_rand_ow_int(ghostpos);

			float xlocal = ghostpos[0] - centre_x;
			float ylocal = ghostpos[1] - centre_y;
			float zlocal = ghostpos[1] - centre_z;
			float ax = abs(xlocal);
			float ay = abs(ylocal);
			float az = abs(zlocal);
			
			bool x_bounded = ax < radius + thickness * gglen;
			bool y_bounded = ay < radius + thickness * gglen;
			bool z_bounded = az < radius + thickness * gglen;
			if ( !((ax > radius && ax < radius + gglen * thickness && y_bounded && z_bounded)
				|| (ay > radius && ay < radius + gglen * thickness && x_bounded && z_bounded)
				|| (az > radius && az < radius + gglen * thickness && x_bounded && y_bounded)) )
				continue;

			std::vector<Particle>::iterator it = static_vec.begin();
			bool obstructed = false;
			while (it != static_vec.end())
			{
				if ( (pow(it->x[0] - ghostpos[0], 2)
				    + pow(it->x[1] - ghostpos[1], 2)
				    + pow(it->x[2] - ghostpos[2], 2) ) < ggradius * ggradius)
				{ 
					obstructed = true;
					break;
				}
				it++;
			}
			if (obstructed) continue;
			else 
			{
				sample_success = true;
				break;
			}
		}

		if (sample_success)
		{
			Particle p (ghostpos[0], ghostpos[1], ghostpos[2]);
			p.ghost_flag = true;
			static_vec.push_back(p);
		}
	}
}
