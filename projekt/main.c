// New version by Ingemar 2010
// Removed all dependencies of the Wild Magic (wml) library.
// Replaced it with VectorUtils2 (in source)
// Replaced old shader module with the simpler "ShaderUtils" unit.

// 2013: Adapted to VectorUtils3 and MicroGlut.

// gcc skinning.c ../common/*.c -lGL -o skinning -I../common

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "common/objects.h"
#ifdef __APPLE__
// Mac
	#include <OpenGL/gl3.h>
	#include "MicroGlut.h"
//uses framework Cocoa
#else
	#ifdef WIN32
// MS
		#include <stdio.h>
		#include <GL/glew.h>
		#include <GL/glut.h>
	#else
// Linux
		#include <GL/gl.h>
		#include "../common/MicroGlut.h" // #include <GL/glut.h>
	#endif
#endif

#include "GL_utilities.h"
#include "VectorUtils3.h"
#include "loadobj.h"

// Ref till shader
GLuint g_shader;

bone_s * bone;
cow_s * cow;

float cam_angle = 0;

mat4 modelViewMatrix, projectionMatrix;

void DisplayWindow()
{
	glClearColor(0.4, 0.4, 0.2, 1);
	glClear(GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT);


	draw_cow(&cow, g_shader);
	draw_bone(&bone, g_shader);
	glutSwapBuffers();
	//printf("runrunrun\n");
}


void OnTimer(int value)
{

	glutTimerFunc(20, &OnTimer, value);


	mat4 proj_matrix = frustum(-1, 1, -1, 1, 2, 1750.0);
	mat4 cam_matrix = lookAt(100*cos(cam_angle), 0,  100*sin(cam_angle), 0, 0, 0, 0.0, 1.0, 0.0);
	//mat4 cam_matrix = lookAt(200+200*cos(cam_angle), 0,  200, 0, 8, 0, 0.0, 1.0, 0.0);

	//g_shader = loadShaders("shader.vert" , "shader.frag");
	glUseProgram(g_shader);

	glUniformMatrix4fv(glGetUniformLocation(g_shader, "proj_matrix"), 1, GL_TRUE, proj_matrix.m);
	glUniformMatrix4fv(glGetUniformLocation(g_shader, "cam_matrix"), 1, GL_TRUE, cam_matrix.m);
	//glutDisplayFunc(DisplayWindow);

	//draw_cow(&cow, g_shader);
	//draw_bone(&bone, g_shader);


	GLfloat t = (GLfloat)glutGet(GLUT_ELAPSED_TIME)/1000.0;
	glUniform1f(glGetUniformLocation(g_shader, "time"), t);

	if(keyIsDown('a'))
	  cam_angle += .2;
	if(keyIsDown('e'))
	  cam_angle -= .2;


	glutPostRedisplay();
}

void keyboardFunc( unsigned char key, int x, int y)
{
// Add any keyboard control you want here
	if(key == 27)	//Esc
		exit(-1);
}

//void init_models()
//{
//  cow_model = LoadModelPlus("res/ko_rough.obj");
//}



/////////////////////////////////////////
//		M A I N
//
int main(int argc, char **argv)
{



	glutInit(&argc, argv);
	initKeymapManager();

	glutInitWindowSize(512, 512);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitContextVersion(3, 2); // Might not be needed in Linux
	glutCreateWindow("Farm Escape");

	glutDisplayFunc(DisplayWindow);
	//glutKeyboardFunc( keyboardFunc ); 
	//glutReshapeFunc(reshape);

	// Set up depth buffer
	glEnable(GL_DEPTH_TEST);

	// initiering
#ifdef WIN32
	glewInit();
#endif

	//mat4 proj_matrix = frustum(-1, 1, -1, 1, 2, 750.0);
	//mat4 cam_matrix = lookAt(0, 0, -20, 0, 0, 0, 0.0, 1.0, 0.0);


	//glUniformMatrix4fv(glGetUniformLocation(g_shader, "proj_matrix"), 1, GL_TRUE, proj_matrix.m);
	//glUniformMatrix4fv(glGetUniformLocation(g_shader, "cam_matrix"), 1, GL_TRUE, cam_matrix.m);
	//glutDisplayFunc(DisplayWindow);

	//init_models();
	create_cow(&cow);
	create_bone(&bone, SetVector(0,0,0));

	g_shader = loadShaders("shader.vert" , "shader.frag");
	glUseProgram(g_shader);

	glutTimerFunc(20, &OnTimer, 0);

	glutMainLoop();
	exit(0);
}
