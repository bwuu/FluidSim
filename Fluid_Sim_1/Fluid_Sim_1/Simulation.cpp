
#include "stdafx.h"

#include "parameters.h"
#include "SPH_Structs.h"
#include "Simul_funcs.h"
#include "marching_cubes.h"
#include "boundary_drawing_tools.h"
#include "GL/glut.h"

void Shut_Down(int return_code);
void Reshape(GLFWwindow* window, int w, int h);
void save_ghosts(std::vector<Particle>& static_vec);

	
std::vector<Particle> static_vec;
std::vector<Particle> p_vec(init_slen_x * init_slen_y * init_slen_z);
std::vector<Cell> c_vec(GRIDSIZE);
GLFWwindow* window;

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

//  init glfw
void init( ) 
{ 
	if(!glfwInit())
        exit(1);

	window = glfwCreateWindow (640, 480, "Simulation", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);

	glClearColor (1.0, 1.0, 1.0, 1.0);
	//glShadeModel (GL_SMOOTH);

    /*
	GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };

	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	GLfloat mat_ambient[] = { 0.7, 0.7, 0.7, 1.0};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, mat_ambient);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
    */

	glfwSetWindowSizeCallback(window, Reshape);
    glPointSize( 3.0 );
}    

//  forward simulation one step - does not display
void update(std::vector<Particle>& p_vec, std::vector<Particle>& static_vec, std::vector<Cell>& c_vec, Container& cont) 
{
	update_cellgrid(c_vec);

	reset_massd_and_forces(p_vec);
	reset_massd_and_forces(static_vec);
	update_mass_densities(c_vec);

	update_pressure_field(p_vec);
	update_pressure_field(static_vec);

	update_forces(c_vec);

	update_velocities(p_vec, cont);
	update_positions(p_vec);

	dumb_collisions(p_vec, cont);
}

//  display current frame of simulation
void display(std::vector<Particle>& p_vec, std::vector<Cell>& c_vec, Container cont)
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	
	gluLookAt(0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

	//glTranslated(-0.0, -0.0, 0.8);

	//glTranslated(cont.x[0], cont.x[2], -cont.x[1]);

	
	glBegin( GL_POINTS );
	glColor3f( 0.4f, 0.3f, 0.5f );
	for ( int i = 0; i < p_vec.size(); i++ )
	{
			glVertex3f( p_vec[i].x[0] - .5, p_vec[i].x[2] - .5, -p_vec[i].x[1]);
	}
	glEnd();

	glBegin(GL_LINES);
	glColor3d(0.0, 0.0, 0.0);
		glVertex3d(0.2, 0.2, -0.2);
		glVertex3d(0.2, 0.8, -0.2);

		glVertex3d(0.2, 0.2, -0.8);
		glVertex3d(0.2, 0.8, -0.8);

		glVertex3d(0.8, 0.2, -0.2);
		glVertex3d(0.8, 0.8, -0.2);

		glVertex3d(0.8, 0.2, -0.8);
		glVertex3d(0.8, 0.8, -0.8);

		glVertex3d(0.2, 0.2, -0.2);
		glVertex3d(0.2, 0.2, -0.8);

		glVertex3d(0.2, 0.2, -0.2);
		glVertex3d(0.8, 0.2, -0.2);

		glVertex3d(0.8, 0.2, -0.8);
		glVertex3d(0.2, 0.2, -0.8);

		glVertex3d(0.8, 0.2, -0.8);
		glVertex3d(0.8, 0.2, -0.2);
	glEnd();

	//glRotated(-90.0, 1.0, 0.0, 0.0);
	//mc_render_image(c_vec);
    std::cout << "swapping buffers";
    glfwSwapBuffers(window);
    std::cout << "swapped";
}



int main()
{
	//  OpenGL state variables
	int simulating = GL_TRUE;
	
	//  set boundaries
	Container bound_box;
	//  uncomment both to enable ghost particles - addt 30 sec loading time
	//sample_ext(static_vec);
	//sample_int(static_vec);

	//  create simulation particles and linked cell grid

	//  initialize particles
	for (int i = 0; i < init_slen_x; i++) {
		for (int j = 0; j < init_slen_y; j++) {
			for (int k = 0; k < init_slen_z; k++) 
			{
				int a = i + j * init_slen_x + k * init_slen_x * init_slen_z;
				p_vec[a].x[0] = .21 + j * init_spacing_x; 
				p_vec[a].x[1] = .21 + i * init_spacing_y;
				p_vec[a].x[2] = .21 + k * init_spacing_z;
			}
		}
	}
	std::cout << " YAY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
	std::cout << p_vec.max_size();

	//  initialize particle positions in linked cell grid
	for (int i = 0; i < p_vec.size(); i++) 
		c_vec[1].ref_list.push_back(&p_vec[i]);

	//  ghost particle grid insertion
	for (int i = 0; i < static_vec.size(); i++) 
		c_vec[1].ref_list.push_back(&static_vec[i]);

	init ();
	float old_time = glfwGetTime();

	while (1) {
		if (simulating == GL_TRUE) 
			update(p_vec, static_vec, c_vec, bound_box);
        /*
		for (int i = 0; i < 3; i++)
			bound_box.x[i] += bound_box.v[i] * dt;
            */

		display(p_vec, c_vec, bound_box);
        /*
		if (glfwGetKey(GLFW_KEY_ESC) == GLFW_PRESS)
			break;
		if (glfwGetKey(GLFW_KEY_UP) == GLFW_PRESS)
			bound_box.v[1] = 4.0;
		else bound_box.v[1] = 0.0;
		if (glfwGetKey(GLFW_KEY_DOWN) == GLFW_PRESS)
			bound_box.v[1] = -4.0;
		if (glfwGetKey(GLFW_KEY_LEFT) == GLFW_PRESS)
			bound_box.v[0] = -4.0;
		else if (glfwGetKey(GLFW_KEY_RIGHT) == GLFW_PRESS)
			bound_box.v[0] = 4.0;
		else bound_box.v[0] = 0.0;

		if (glfwGetMouseButton(GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
			if (simulating == GL_FALSE) simulating = GL_TRUE;
			//  else simulating = GL_FALSE;
		}
        */
	}
	Shut_Down(0);
	return 0;
}

//  shut down glfw
void Shut_Down(int return_code)
{
  glfwTerminate();
  exit(return_code);
}

//  reshape func
void Reshape(GLFWwindow* window, int w, int h)
{
   glViewport (0, 0, (GLsizei) w, (GLsizei) h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glFrustum (-1.0, 1.0, -1.0, 1.0, 3.5, 20.0);
   glMatrixMode(GL_MODELVIEW);
}

void save_ghosts(std::vector<Particle>& static_vec)
{
	std::ofstream myfile ("boundary_particles.txt", std::ios::trunc);
	if (!myfile.is_open()) 
	{
		std::cout << "error: failed to open boundary_particles.txt";
		return;
	}

}

