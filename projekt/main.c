// New version by Ingemar 2010
// Removed all dependencies of the Wild Magic (wml) library.
// Replaced it with VectorUtils2 (in source)
// Replaced old shader module with the simpler "ShaderUtils" unit.

// 2013: Adapted to VectorUtils3 and MicroGlut.

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

bone_s bone;
joint_s foot_joint[4];
joint_s knee_joint[4];
joint_s legbase_joint[4];
joint_s tail_joint[4];
joint_s head_joint[3];

cow_s cow;

float cam_angle = 0;
float cam_dist = 6;

mat4 modelViewMatrix, projectionMatrix;

void DisplayWindow()
{
	glClearColor(0.4, 0.4, 0.2, 1);
	glClear(GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT);

	glUniform1i(glGetUniformLocation(g_shader, "draw_cow"), 1);
	draw_cow(&cow, g_shader);
	glUniform1i(glGetUniformLocation(g_shader, "draw_cow"), 0);
	int i;
	for(i=0; i<4;i++)
	{
	  draw_joint(&foot_joint[i], g_shader);
	  draw_joint(&knee_joint[i], g_shader);
	  draw_joint(&legbase_joint[i], g_shader);
	  draw_joint(&tail_joint[i], g_shader);
	  if(i < 3) 
	    draw_joint(&head_joint[i], g_shader);

	}

	//draw_bone(&bone, g_shader);
	glutSwapBuffers();
	//printf("runrunrun\n");
}


void OnTimer(int value)
{


	glutTimerFunc(20, &OnTimer, value);

	mat4 proj_matrix = frustum(-1, 1, -1, 1, 1, 750.0);
	mat4 cam_matrix = lookAt(cam_dist*cos(cam_angle), 1,  cam_dist*sin(cam_angle), 0, 0, 0, 0.0, 1.0, 0.0);
	//mat4 cam_matrix = lookAt(200+200*cos(cam_angle), 0,  200, 0, 8, 0, 0.0, 1.0, 0.0);

	//g_shader = loadShaders("shader.vert" , "shader.frag");
	//glUseProgram(g_shader);

	glUniformMatrix4fv(glGetUniformLocation(g_shader, "proj_matrix"), 1, GL_TRUE, proj_matrix.m);
	glUniformMatrix4fv(glGetUniformLocation(g_shader, "cam_matrix"), 1, GL_TRUE, cam_matrix.m);
	//glutDisplayFunc(DisplayWindow);

	//draw_cow(&cow, g_shader);
	//draw_bone(&bone, g_shader);


	GLfloat t = (GLfloat)glutGet(GLUT_ELAPSED_TIME)/1000.0;
	glUniform1f(glGetUniformLocation(g_shader, "time"), t);


	mat4 tmp, testjoint, testjoint2;
	testjoint = IdentityMatrix();
	legbase_joint[3].T = T(.1, 2, 0);
	legbase_joint[3].M = T(.1, 2, 0);

	knee_joint[3].T = T(.4, .8, 0);
	knee_joint[3].M = T(.4, .8, 0);

	knee_joint[3].Minv = InvertMat4(T(.4, .8, 0));
	legbase_joint[3].Minv = InvertMat4(T(.1, 2, 0));

	legbase_joint[3].R = ArbRotate(SetVector(0,0,1), sin(5*t)/2);
	legbase_joint[3].Mp = Mult(legbase_joint[3].T, legbase_joint[3].R);

	testjoint2 = Mult(legbase_joint[3].Mp, legbase_joint[3].Minv);
	knee_joint[3].R = ArbRotate(SetVector(0,0,1), cos(4*t)/2);
	knee_joint[3].M = Mult(knee_joint[3].T, legbase_joint[3].Minv);

	knee_joint[3].Mp = Mult(knee_joint[3].T, knee_joint[3].R);
	knee_joint[3].Mp = Mult(legbase_joint[3].Minv, knee_joint[3].Mp);
	knee_joint[3].Minv = InvertMat4(knee_joint[3].M);
	//testjoint = testjoint2;
	mat4 inverts = Mult(knee_joint[3].Minv, legbase_joint[3].Minv);
	mat4 primes = Mult(legbase_joint[3].Mp, knee_joint[3].Mp);
	testjoint = Mult(primes, inverts);
	//testjoint = testjoint2;
	//testjoint = Mult(Mult(legbase_joint[3].Mp, knee_joint[3].Mp),
	//Mult(knee_joint[3].Minv, legbase_joint[3].Minv));

 	/*
	mat4 M1, M2, M1_p, M2_p;

	testjoint = IdentityMatrix();
	M1_p = Mult( T(.1,2.0,0), ArbRotate(SetVector(0,0,1), sin(5*t)/2));
	M1 = T(.1,2.0,0);
	testjoint2 = Mult(M1_p, InvertMat4(M1));


	M2 = Mult(T(.4,.8,0), InvertMat4(M1));
	M2_p = Mult(T(.4,.8,0), ArbRotate(SetVector(0,0,1), cos(4*t)/2));
	M2_p = Mult(InvertMat4(M1), M2_p);

	mat4 inverts = Mult(InvertMat4(M2), InvertMat4(M1));
	mat4 primes = Mult(M1_p, M2_p);
	testjoint = Mult(primes, inverts);
*/

	glUniformMatrix4fv(glGetUniformLocation(g_shader, "testjoint"), 1, GL_TRUE, testjoint.m);
	glUniformMatrix4fv(glGetUniformLocation(g_shader, "testjoint2"), 1, GL_TRUE, testjoint2.m);


	if(keyIsDown('e'))
	  cam_angle += .05;
	if(keyIsDown('i'))
	  cam_angle -= .05;

	if(keyIsDown('u'))
	  cam_dist += .05;

	if(keyIsDown('p'))
	  cam_dist -= .05;


	glutPostRedisplay();
}

void keyboardFunc( unsigned char key, int x, int y)
{
// Add any keyboard control you want here
	if(key == 27)	//Esc
		exit(-1);
}


/////////////////////////////////////////
//		M A I N
//
int main(int argc, char **argv)
{

	glutInit(&argc, argv);
	initKeymapManager();

	//glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(512, 512);
	glutInitContextVersion(3, 2); // Might not be needed in Linux
	glutCreateWindow("Farm Escape");
	glutDisplayFunc(DisplayWindow);

	create_joint(&legbase_joint[0], SetVector(-2.3, 1.8, -.4));
	create_joint(&legbase_joint[1], SetVector(-2.3, 1.8, .4));
	create_joint(&legbase_joint[2], SetVector(.1, 2, -.4));
	create_joint(&legbase_joint[3], SetVector(.1, 2.0, 0));

	create_joint(&knee_joint[0], SetVector(-2.6, .9, -.4));
	create_joint(&knee_joint[1], SetVector(-2.6, .9, .4));
	create_joint(&knee_joint[2], SetVector(.4, .9, -.4));
	create_joint(&knee_joint[3], SetVector(.4, .8, 0));

	create_joint(&foot_joint[0], SetVector(.34, 0, -.4));
	create_joint(&foot_joint[1], SetVector(.34, 0, .4));
	create_joint(&foot_joint[2], SetVector(-2.75, 0, -.4));
	create_joint(&foot_joint[3], SetVector(-2.75, 0, .4));


	//SET PARENTS!
	foot_joint[0].parent = &knee_joint[0];
	foot_joint[1].parent = &knee_joint[1];
	foot_joint[2].parent = &knee_joint[2];
	foot_joint[3].parent = &knee_joint[3];

	knee_joint[0].parent = &legbase_joint[0];
	knee_joint[1].parent = &legbase_joint[1];
	knee_joint[2].parent = &legbase_joint[2];
	knee_joint[3].parent = &legbase_joint[3];


	//TAIL JOINTS
	create_joint(&tail_joint[0], SetVector(.7, 3.75, 0));
	create_joint(&tail_joint[1], SetVector(2.2, 3.8, 0));
	create_joint(&tail_joint[2], SetVector(3.0, 3.75, 0));
	create_joint(&tail_joint[3], SetVector(3.9, 3.65, 0));

	//HEAD JOINTS
	create_joint(&head_joint[0], SetVector(-2.9, 3.2, 0));
	create_joint(&head_joint[1], SetVector(-3.85, 4, 0));
	create_joint(&head_joint[2], SetVector(-4.7, 3, 0));


	create_cow(&cow);


	g_shader = loadShaders("shader.vert" , "shader.frag");
	glUseProgram(g_shader);

	//glutKeyboardFunc( keyboardFunc ); 
	//glutReshapeFunc(reshape);

	// Set up depth buffer
	glEnable(GL_DEPTH_TEST);

	//glEnable (GL_BLEND);
	//glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);



	// initiering
//#ifdef WIN32
//	glewInit();
//#endif

	glutTimerFunc(20, &OnTimer, 0);

	glutMainLoop();
	exit(0);
}
