
#include "loadShader.h"

#include "parameters.h"
#include "SPH_Structs.h"
#include "Simul_funcs.h"
//#include "marching_cubes.h"
//#include "boundary_drawing_tools.h"
#include "GL/glut.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

void Shut_Down(int return_code);
void Reshape(GLFWwindow* window, int w, int h);
void save_ghosts(std::vector<Particle>& static_vec);

std::vector<Particle> static_vec;
std::vector<Particle> p_vec(init_slen_x * init_slen_y * init_slen_z);
std::vector<Cell> c_vec(GRIDSIZE);

GLFWwindow* window;
GLuint programID;
GLuint VertexArrayID;
GLuint vertexbuffer;


static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
	{
        std::cout << "exiting...\n";
        glfwSetWindowShouldClose(window, GL_TRUE);
	}
}

//  init glfw
void init( ) 
{ 
	if(!glfwInit())
        exit(1);

	window = glfwCreateWindow (640, 480, "Simulation", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);

	glClearColor (0.5, 0.8, 0.8, 1.0);        // Initialize GLEW

    glewExperimental = true; // Needed for core profile
    if (glewInit() != GLEW_OK) {
            fprintf(stderr, "Failed to initialize GLEW\n");
            return;
    }

    glEnable(GL_POINT_SPRITE);
    glEnable(GL_PROGRAM_POINT_SIZE);
    //glEnable(GL_DEPTH_TEST);

    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);
	glfwSetWindowSizeCallback(window, Reshape);
    glGenBuffers(1, &vertexbuffer);

    programID = LoadShaders("../Fluid_Sim_1/PointSprites.vertexshader", "../Fluid_Sim_1/PointSprites.fragmentshader");
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

float* posdata;
//  display current frame of simulation
void display(std::vector<Particle>& p_vec, std::vector<Cell>& c_vec, Container cont, int index)
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glUseProgram(programID);

    //  vertex shader
	glm::mat4 modelview = glm::mat4(1.0);
	modelview = glm::rotate<float>(modelview, 30.0, glm::vec3(1, 0, 0));
	modelview = glm::rotate<float>(modelview, 30.0, glm::vec3(0, 1, 0));
    modelview = glm::translate<float>(modelview, glm::vec3(0.5, 0, -0.3));
    GLint modelview_uniform_loc = glGetUniformLocation(programID, "modelview");
	glUniformMatrix4fv(modelview_uniform_loc, 1, GL_FALSE, &modelview[0][0]);

	glm::mat4 projection = glm::mat4(1.0);
	projection = glm::perspective(60.0f, 4.0f/3.0f, 0.1f, 100.0f); 
    GLint projection_uniform_loc = glGetUniformLocation(programID, "projection");
	glUniformMatrix4fv(projection_uniform_loc, 1, GL_FALSE, &projection[0][0]);

    GLint screenSize_uniform_loc = glGetUniformLocation(programID, "screenSize");
	glUniform2f(screenSize_uniform_loc, 640, 480);

    GLint spriteSize_uniform_loc = glGetUniformLocation(programID, "spriteSize");
	glUniform1f(spriteSize_uniform_loc, 0.035f);

    //  fragment shader
    GLint color_uniform_loc = glGetUniformLocation(programID, "Color");
    GLint lightdir_uniform_loc = glGetUniformLocation(programID, "lightDir");
    glUniform3f(color_uniform_loc, 0.9, 0.9, 0.9); 
    glUniform3f(lightdir_uniform_loc, 0.6, -0.6, 1); 

	//gluLookAt(0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	//glTranslated(cont.x[0], cont.x[2], -cont.x[1]);

	for (int i = 0; i < p_vec.size() * 4; i+=4) 
	{
        posdata[i]   = p_vec[i / 4].x[0] - .5;
        posdata[i+1] = p_vec[i / 4].x[2] - .5;
		posdata[i+2] = -p_vec[i / 4].x[1];
	}
    //posdata[p_vec.size()*4]   = p_vec[index].x[0] - .5;
    //posdata[p_vec.size()*4+1] = p_vec[index].x[1] - .5;
    //posdata[p_vec.size()*4+2] = -p_vec[index].x[2];

    glBufferData(GL_ARRAY_BUFFER, p_vec.size() * sizeof(float), posdata, GL_DYNAMIC_DRAW);

	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glVertexAttribPointer(
		0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
		4,                  // size
		GL_FLOAT,           // type
		GL_FALSE,           // normalized?
		0,                  // stride
		(void*)0            // array buffer offset
		);

	// Draw the points !
	glDrawArrays(GL_POINTS, 0, p_vec.size());
    //glUniform3f(color_uniform_loc, 0.0, 0.0, 0.0); 
	//glDrawArrays(GL_POINTS, p_vec.size()*4, 1);

	glDisableVertexAttribArray(0);

	// Swap buffers
	glfwSwapBuffers(window);
	glfwPollEvents();
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
				p_vec[a].x[2] = .31 + k * init_spacing_z;
			}
		}
	}
	std::cout << " Initialized correctly \n";
	std::cout << p_vec.max_size();

	//  initialize particle positions in linked cell grid
	for (int i = 0; i < p_vec.size(); i++) 
		c_vec[1].ref_list.push_back(&p_vec[i]);

	//  ghost particle grid insertion
	for (int i = 0; i < static_vec.size(); i++) 
		c_vec[1].ref_list.push_back(&static_vec[i]);

	posdata = new float[p_vec.size() * 4 + 4];
    for (int i = 3; i < p_vec.size() * 4 + 4; i+=4)
        posdata[i] = 1.0;
	init ();
	float old_time = glfwGetTime();
    
    int counter = 0;

	while (1) {
        counter++;
        if (counter % 100 == 0) std::cout << counter << "\n";

        /*
        float maxx = 0;
        int index;
        for (int i = 0; i < p_vec.size(); i++)
		{
            if (p_vec[i].x[2] > maxx) 
            {
                maxx = p_vec[i].x[2];
                index = i;
			}
		}
        std::cout << "I am max. My height is " << maxx << "==" << p_vec[index].x[2] << "\n";
        std::cout << "my index is " << index << "\n";
        std::cout << "presssure:\n";
        std::cout << p_vec[index].f_pressure[0] << "\n";
        std::cout << p_vec[index].f_pressure[1] << "\n";
        std::cout << p_vec[index].f_pressure[2] << "\n";
        std::cout << "viscosity:\n";
        std::cout << p_vec[index].f_viscosity[0] << "\n";
        std::cout << p_vec[index].f_viscosity[1] << "\n";
        std::cout << p_vec[index].f_viscosity[2] << "\n";
        std::cout << "velocity:\n";
        std::cout << p_vec[index].v[0] << "\n";
        std::cout << p_vec[index].v[1] << "\n";
        std::cout << p_vec[index].v[2] << "\n";
        std::cout << "\n";
        */

		if (simulating == GL_TRUE) 
			update(p_vec, static_vec, c_vec, bound_box);

		display(p_vec, c_vec, bound_box, 1);
	}
	Shut_Down(0);
    delete[] posdata;
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

