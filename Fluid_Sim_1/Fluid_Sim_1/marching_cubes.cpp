#include "Parameters.h"
#include "SPH_structs.h"
#include "stdafx.h"
#include "mc_tables.h"

//  ---------------------------------------------------------------
//  parameters/constants for Marching Cubes rendering
//  ---------------------------------------------------------------

//  attempt to optimize calculations
const double mc_term = 315 / 64 / PI / pow(h, 9);

//  rendering resolution
const double mclen = 0.018;

//  threshold for mcubes. higher == easier for particles to be considered on the "surface"
const double thresh = .45;

//  various mc grid properties, determined by resolution
const int mcx = (int) (ceil((SIZE_X / mclen) + 0.1) + 0.5) + 1; 
const int mcy = (int) (ceil((SIZE_Y / mclen) + 0.1) + 0.5) + 1; 
const int mcz = (int) (ceil((SIZE_Z / mclen) + 0.1) + 0.5) + 1; 
const int mcsize = mcx * mcy * mcz;

//  ---------------------------------------------------------------
//  structs for mc
//  ---------------------------------------------------------------
//  for storing vertexes ready for mc algorithm or for
//  OpenGl rendering
struct mc_vertex {
	double pos[3];
	double val;
	double normal[3];

	mc_vertex ()
	{
		for (int i = 0; i < 3; i++)
		{
			pos[i] = 0;
			normal[i] = 0;
		}
		val = 0;
	}
};

//  stores 8 vertexes, ready to pass onto the Polygonise function
struct GRIDCELL {
   mc_vertex vert[8];
};



//  ---------------------------------------------------------------
//  some peripheral functions
//  ---------------------------------------------------------------

mc_vertex VertexInterp(double isolevel, mc_vertex& vert1, mc_vertex& vert2);

int Polygonise(GRIDCELL& grid,double isolevel, mc_vertex* triangles);

//  return actual value of position along a certain axes, given index into mc grid
double mc_posx(int i) {  return (i % mcx) * mclen;  }

double mc_posy(int i) {  return ( (i % (mcx * mcy)) / mcx) * mclen;  }

double mc_posz(int i) {  return (i / (mcx * mcy)) * mclen;  }


//  ---------------------------------------------------------------
//  main functions
//  ---------------------------------------------------------------

//  iterate color field calculations for a vertex over a cell
void mc_iterate_cf_over_cell(double* vertex, double pos[], Cell& acell)
{
	//  get iterator for cell
	std::list<Particle*>::iterator it = acell.ref_list.begin();

	//  iterate over cell
	while (it != acell.ref_list.end()) 
	{
		const double sqr_dist = pow(pos[0] - (*it)->x[0], 2) 
			+ pow(pos[1] - (*it)->x[1], 2)
			+ pow(pos[2] - (*it)->x[2], 2);

		(*vertex) += (sqr_dist > h*h) ? 0 : (mc_term 
			* pow(h*h-sqr_dist, 3) * mass / (*it)->mass_density);

		assert(*vertex <= DBL_MAX && *vertex >= -DBL_MAX);
		it++;
	}
	return;
}

//  render and display image using all particles in c_vec
void mc_render_image(std::vector<Cell>& c_vec)
{
	std::vector<mc_vertex> mc_field (mcsize);
	for (int i = 0; i < mcsize; i++)
	{
		mc_field[i].pos[0] = mc_posx(i);
		mc_field[i].pos[1] = mc_posy(i);
		mc_field[i].pos[2] = mc_posz(i);

		int grid_num = ( (int) (mc_field[i].pos[0] / h + 0.0000001) )
				+ ( (int) (mc_field[i].pos[1] / h + 0.0000001) ) * GRID_X
				+ ( (int) (mc_field[i].pos[2] / h + 0.0000001) ) * GRID_Y * GRID_X;

		assert(grid_num < GRIDSIZE);
		assert(grid_num >= 0);

		if (c_vec[grid_num].ref_list.size() > 0)
		{
			mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num]);

			bool y_down = (grid_num % (GRID_X * GRID_Y)) / GRID_X != 0;
			bool y_up = (grid_num % (GRID_X * GRID_Y)) / GRID_X >= 1;
			bool x_down = (grid_num % GRID_X) != 0;
			bool x_up = (grid_num % GRID_X) < (GRID_X - 1);
			bool z_down = (grid_num >= GRID_X * GRID_Y);
			bool z_up = grid_num < GRID_X * GRID_Y * (GRID_Z - 1);
			if (x_up) 
			{
				mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + 1]);
				if (y_up) 
				{
					mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X + 1]);
					if (z_up) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X * GRID_Y + GRID_X + 1]);
					if (z_down) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X * GRID_Y + GRID_X + 1]);
				}
				if (y_down)
				{
					mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X + 1]);
					if (z_up) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X * GRID_Y - GRID_X + 1]);
					if (z_down) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X * GRID_Y - GRID_X + 1]);
				}
				if (z_up) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X * GRID_Y + 1]);
				if (z_down) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X * GRID_Y + 1]);
			}
			if (x_down) 
			{
				mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - 1]);
				if (y_up) 
				{
					mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X - 1]);
					if (z_up) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X * GRID_Y + GRID_X - 1]);
					if (z_down) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X * GRID_Y + GRID_X - 1]);
				}
				if (y_down)
				{
					mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X - 1]);
					if (z_up) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X * GRID_Y - GRID_X - 1]);
					if (z_down) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X * GRID_Y - GRID_X - 1]);
				}
				if (z_up) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X * GRID_Y - 1]);
				if (z_down) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X * GRID_Y - 1]);
			}

			if (y_up) 
			{
				mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X]);
				if (z_up) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X * GRID_Y + GRID_X]);
				if (z_down) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X * GRID_Y + GRID_X]);
			}
			if (y_down)
			{
				mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X]);
				if (z_up) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X * GRID_Y - GRID_X]);
				if (z_down) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X * GRID_Y - GRID_X]);
			}
			if (z_up) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num + GRID_X * GRID_Y]);
			if (z_down) mc_iterate_cf_over_cell(&(mc_field[i].val), mc_field[i].pos, c_vec[grid_num - GRID_X * GRID_Y]);
		}
	}
	
	for (int i = 0; i < mcsize; i++)
	{
		if (i % mcx == mcx - 1 // x-top
			|| i % mcx == 0    // x-low
			|| i % (mcx * mcy) >= mcx * (mcy - 1) // y-top
			|| i % (mcx * mcy) < mcx // y-low
			|| i >= mcx * mcy * (mcz - 1)
			|| i / (mcx * mcy) == 0
			) continue;

		mc_field[i].normal[0] = (mc_field[i - 1].val - mc_field[i + 1].val) / mcx; 
		mc_field[i].normal[1] = (mc_field[i - mcx].val - mc_field[i + mcx].val) / mcy; 
		mc_field[i].normal[2] = (mc_field[i - mcx * mcy].val - mc_field[i + mcx * mcy].val) / mcz; 

		double norm = sqrt(pow(mc_field[i].normal[0], 2) 
			+ pow(mc_field[i].normal[1], 2)
			+ pow(mc_field[i].normal[2], 2));

		if (norm > 0.0000001)
		{
			mc_field[i].normal[0] /= norm;
			mc_field[i].normal[1] /= norm;
			mc_field[i].normal[2] /= norm;
		}
	}

	mc_vertex triangles[15];
	GRIDCELL gcell;

	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_shininess[] = { 50.0 };
	
	GLfloat mat_color1[] = { 0.5, 0.5, 0.5, 0.5 };

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_color1);

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	for (int i = 0; i < mcsize; i++)
	{
		if (i % mcx == mcx - 1 // x-top
			|| i % mcx == 0    // x-low
			|| i % (mcx * mcy) >= mcx * (mcy - 1) // y-top
			|| i % (mcx * mcy) < mcx // y-low
			|| i >= mcx * mcy * (mcz - 1)
			|| i / (mcx * mcy) == 0
			) continue;

		gcell.vert[3] = mc_field[i];

		gcell.vert[0] = mc_field[i + mcx];

		gcell.vert[1] = mc_field[i + mcx + 1];

		gcell.vert[2] = mc_field[i + 1];

		gcell.vert[7] = mc_field[i + mcx * mcy];

		gcell.vert[4] = mc_field[i + mcx * mcy + mcx];

		gcell.vert[5] = mc_field[i + mcx * mcy + mcx + 1];

		gcell.vert[6] = mc_field[i + mcx * mcy + 1];

		int ntri = Polygonise(gcell, thresh, triangles);

		glBegin( GL_TRIANGLES );
		for ( int j = 0; j < ntri; j += 3 )
		{
			glNormal3f( triangles[j].normal[0], triangles[j].normal[1], triangles[j].normal[2]);
			glVertex3f( triangles[j].pos[0], triangles[j].pos[1], triangles[j].pos[2]);

			glNormal3f( triangles[j + 1].normal[0], triangles[j + 1].normal[1], triangles[j + 1].normal[2]);
			glVertex3f( triangles[j + 1].pos[0], triangles[j + 1].pos[1], triangles[j + 1].pos[2]);

			glNormal3f( triangles[j + 2].normal[0], triangles[j + 2].normal[1], triangles[j + 2].normal[2]);
			glVertex3f( triangles[j + 2].pos[0], triangles[j + 2].pos[1], triangles[j + 2].pos[2]);
		}
		glEnd();
	}
	
	return;
}




/*
   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangles"
   will be loaded up with the vertices at most 5 triangular facets.
	0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/
int Polygonise(GRIDCELL& grid, double isolevel, mc_vertex* triangles)
{
   int ntriang;
   int cubeindex;
   mc_vertex vertlist[12];

   /*
      Determine the index into the edge table which
      tells us which vertices are inside of the surface
   */
cubeindex = 0;
if (grid.vert[0].val < isolevel) cubeindex |= 1;
if (grid.vert[1].val < isolevel) cubeindex |= 2;
if (grid.vert[2].val < isolevel) cubeindex |= 4;
if (grid.vert[3].val < isolevel) cubeindex |= 8;
if (grid.vert[4].val < isolevel) cubeindex |= 16;
if (grid.vert[5].val < isolevel) cubeindex |= 32;
if (grid.vert[6].val < isolevel) cubeindex |= 64;
if (grid.vert[7].val < isolevel) cubeindex |= 128;

/* Cube is entirely in/out of the surface */
if (edgeTable[cubeindex] == 0)
	return(0);

/* Find the vertices where the surface intersects the cube */
if (edgeTable[cubeindex] & 1)
	vertlist[0] =
	VertexInterp(isolevel,grid.vert[0],grid.vert[1]);
if (edgeTable[cubeindex] & 2)
	vertlist[1] =
	VertexInterp(isolevel,grid.vert[1],grid.vert[2]);
if (edgeTable[cubeindex] & 4)
	vertlist[2] =
	VertexInterp(isolevel,grid.vert[2],grid.vert[3]);
if (edgeTable[cubeindex] & 8)
	vertlist[3] =
	VertexInterp(isolevel,grid.vert[3],grid.vert[0]);
if (edgeTable[cubeindex] & 16)
	vertlist[4] =
	VertexInterp(isolevel,grid.vert[4],grid.vert[5]);
if (edgeTable[cubeindex] & 32)
	vertlist[5] =
	VertexInterp(isolevel,grid.vert[5],grid.vert[6]);
if (edgeTable[cubeindex] & 64)
	vertlist[6] =
	VertexInterp(isolevel,grid.vert[6],grid.vert[7]);
if (edgeTable[cubeindex] & 128)
	vertlist[7] =
	VertexInterp(isolevel,grid.vert[7],grid.vert[4]);
if (edgeTable[cubeindex] & 256)
	vertlist[8] =
	VertexInterp(isolevel,grid.vert[0],grid.vert[4]);
if (edgeTable[cubeindex] & 512)
	vertlist[9] =
	VertexInterp(isolevel,grid.vert[1],grid.vert[5]);
if (edgeTable[cubeindex] & 1024)
	vertlist[10] =
	VertexInterp(isolevel,grid.vert[2],grid.vert[6]);
if (edgeTable[cubeindex] & 2048)
	vertlist[11] =
	VertexInterp(isolevel,grid.vert[3],grid.vert[7]);

/* Create the triangle */
ntriang = 0;
for (int i = 0; triTable[cubeindex][i] != -1; i += 3) {
	triangles[ntriang    ] = vertlist[triTable[cubeindex][i  ]];
	triangles[ntriang + 1] = vertlist[triTable[cubeindex][i+1]];
	triangles[ntriang + 2] = vertlist[triTable[cubeindex][i+2]];
	ntriang += 3;
}

return(ntriang);
}

/*
Linearly interpolate the position where an isosurface cuts
an edge between two vertices, each with their own scalar value
*/
mc_vertex VertexInterp(double isolevel, mc_vertex& vert1, mc_vertex& vert2) 
{
	double mu;
	mc_vertex ans;

	if (abs(isolevel - vert1.val) < 0.00001)
		return(vert1);
	if (abs(isolevel - vert2.val) < 0.00001)
		return(vert2);
	if (abs(vert1.val - vert2.val) < 0.00000001)
		return(vert1);
	mu = (isolevel - vert1.val) / (vert2.val - vert1.val);
	ans.pos[0] = vert1.pos[0] + mu * (vert2.pos[0] - vert1.pos[0]);
	ans.pos[1] = vert1.pos[1] + mu * (vert2.pos[1] - vert1.pos[1]);
	ans.pos[2] = vert1.pos[2] + mu * (vert2.pos[2] - vert1.pos[2]);

	ans.normal[0] = vert1.normal[0] + mu * (vert2.normal[0] - vert1.normal[0]);
	ans.normal[1] = vert1.normal[1] + mu * (vert2.normal[1] - vert1.normal[1]);
	ans.normal[2] = vert1.normal[2] + mu * (vert2.normal[2] - vert1.normal[2]);

	const double norm = sqrt(pow(ans.normal[0],2) 
		+ pow(ans.normal[1],2) 
		+ pow(ans.normal[2],2));
	
	if (norm > 0.0000001)
	{
	ans.normal[0] /= norm;
	ans.normal[1] /= norm;
	ans.normal[2] /= norm;
	}
	else {
		ans.normal[0] = 1.0;
	}
	return(ans);
}

