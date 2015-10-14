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

	joint_s * j = &tail_joint[0];
	joint_s * jc = j->child;
	joint_s * jcc = jc->child;
	//joint_s * jp = j->parent;
	//joint_s * jpp = jp->parent;
/*
	jp->R = ArbRotate(SetVector(0,0,1), cos(4*t)*1.2);
	jp->M = Mult(jp->T, jpp->Minv);

	jp->Mp = Mult(jp->T, jp->R);
	jp->Mp = Mult(jpp->Minv, jp->Mp);
	jp->Minv = InvertMat4(jp->M);
*/
	float Ms[8][4*4];
	float currpos[8*3] = {{0}};
	float bonepos[8*3] = {{0}};
	//GLfloat * Mtmp = Ms;
	int i=0, ii=0;
	j->R = ArbRotate(SetVector(0,0,1), cos(4*t/(i+1))/1.5);
	//j->R = ArbRotate(SetVector(0,0,1), 0);
	jc->R = ArbRotate(SetVector(0,0,1), sin(4*t/(3+1))/2);
	jcc->R = ArbRotate(SetVector(0,0,1), cos(8*t/(2+1))/4);
	mat4 Mpacc, Minvacc, tmptrans, tmppp, invtrans;
	Mpacc = IdentityMatrix();
	Minvacc = IdentityMatrix();
	tmppp = IdentityMatrix();

	//Mpacc = InvertMat4(Mult(j->T, j->R));
	//Minvacc = InvertMat4(j->Minv);
	while(j->child != NULL)
	{
	  //jp = j->child;
	  //j->R = ArbRotate(SetVector(0,0,1), cos(4*t/(i+1))/4);
	  //j->M = Mult(j->T, jp->Minv);
	  vec3 tmp_bonepos;
	  tmptrans = j->T;
	  tmppp = Mult(tmppp, tmptrans);
	  invtrans = InvertMat4(tmptrans);
	  tmppp = Mult(tmppp, Mult(j->R, invtrans));

	  tmp_bonepos = ScalarMult(
		VectorAdd(j->pos, jc->pos), .5); //middle of bone
	  //tmp_bonepos = MultVec3(tmppp, tmp_bonepos);
	  //j->Mp = Mult(j->T, j->R);
	  //jp->Mp = Mult(jp->T, jp->R);
	  //j->Mp = Mult(jp->Minv, j->Mp); //super special?

	  //jp->R = ArbRotate(SetVector(0,0,1), sin(t*(i+1)/6)/4);
	  //jp->Mp = Mult(jp->T, jp->R);

	  //j->Mp = Mult(jp->Mp, j->Mp);
	  //j->Minv = Mult(j->Minv, jp->Minv);
	  //j->Mp = Mult(Mult(jp->Mp, j->Mp), Mult(j->Minv, jp->Minv));
	  //Mpacc = Mult(Mpacc, j->Mp);
	  //Minvacc = Mult(j->Minv, Minvacc);
	  Mpacc = Mult(Mpacc, j->Mp);
	  Minvacc = Mult(Minvacc, j->Minv);
	  //j->Mp = Mult(Mult(jp->Mp, jp->Minv), j->Mp);
	  //j->Mp = Mult(jp->Mp, jp->Minv);
	  //j->Mp = Mult(jp->Mp, jp->Minv);
	  //jp->Mtot = Mult(jp->Mp, jp->Minv);
	  mat4 tmpp = Mult(Mpacc, Minvacc);
	  //mat4 tmpp = Mult(j->Mp, j->Minv);

	  for(ii=0;ii<16;ii++)
	    Ms[i][ii] = (tmppp).m[ii];
	  //for(ii=0;ii<16;ii++)
	  //  Ms[1][ii] = (jp->Mp).m[ii];
	  //for(ii=0;ii<16;ii++)
	  //  Ms[2][ii] = (jpp->Mtot).m[ii];

	  //glUniformMatrix4fv(glGetUniformLocation(g_shader, "testjoint2"), 1, GL_TRUE, testjoint2.m);


	  currpos[i*3] = 10*j->pos.x;
	  currpos[i*3+1] = 10*j->pos.y;
	  currpos[i*3+2] = 10*j->pos.z;

	  bonepos[i*3] = 10*tmp_bonepos.x;
	  bonepos[i*3+1] = 10*tmp_bonepos.y;
	  bonepos[i*3+2] = 10*tmp_bonepos.z;

	  //printf("%i\n", i);
	  i++;

	  j = j->child;
	}
	  glUniformMatrix4fv(glGetUniformLocation(g_shader, "testjoint"), 8, GL_TRUE, Ms);
	  glUniform3fv(glGetUniformLocation(g_shader, "currpos"), 8, currpos);
	  glUniform3fv(glGetUniformLocation(g_shader, "bonepos"), 8, bonepos);

	//test with tail!
/*
	//start at the end of the tail
	j = &tail_joint[3];

	while(j->parent != NULL)
	{
	  jp = j->parent;
	  j->R = ArbRotate(SetVector(0,0,1), cos(4*t)*1.2);
	  j->M = Mult(j->T, jp->Minv);
	  j->Mp = Mult(j->T, j->R);
	  j->Mp = Mult(jp->Minv, j->Mp);
	  j->Minv = InvertMat4(j->M);

	  jp->R = ArbRotate(SetVector(0,0,1), sin(t/4));
	  jp->Mp = Mult(jp->T, jp->R);

	  j->Mtot = Mult(Mult(jp->Mp, j->Mp), Mult(j->Minv, jp->Minv));
	  jp->Mtot = Mult(jp->Mp, jp->Minv);

	  glUniformMatrix4fv(glGetUniformLocation(g_shader, "tailjoint"), 1, GL_TRUE, j->Mtot.m);
	  glUniformMatrix4fv(glGetUniformLocation(g_shader, "tailjoint2"), 1, GL_TRUE, jp->Mtot.m);

	  glUniform3f(glGetUniformLocation(g_shader, "currtailpos"), 10*legbase_joint[0].pos.x, 10*legbase_joint[0].pos.y, 10*legbase_joint[0].pos.z);
	  glUniform3f(glGetUniformLocation(g_shader, "currtailpos2"), 10*knee_joint[0].pos.x, 10*knee_joint[0].pos.y, 10*knee_joint[0].pos.z);

	  j = j->parent;
	}

*/


/*
	jpp->R = ArbRotate(SetVector(0,0,1), sin(t/4));
	//legbase_joint[0].R = ArbRotate(SetVector(0,0,1), 0);
	jpp->Mp = Mult(jpp->T, jpp->R);

	testjoint2 = Mult(jpp->Mp, jpp->Minv);

	mat4 inverts = Mult(jp->Minv, jpp->Minv);
	mat4 primes = Mult(jpp->Mp, jp->Mp);
	testjoint = Mult(primes, inverts);
*/

	//testjoint = j->Mtot;
	//testjoint2 = jp->Mtot;

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

	create_joint(&legbase_joint[0], SetVector(-2.3, 1.8, .4));
	create_joint(&legbase_joint[1], SetVector(-2.3, 1.8, -.4));
	create_joint(&legbase_joint[2], SetVector(.1, 2, -.4));
	create_joint(&legbase_joint[3], SetVector(.1, 2.0, 0));

	create_joint(&knee_joint[0], SetVector(-2.6, .9, .4));
	create_joint(&knee_joint[1], SetVector(-2.6, .9, -.4));
	create_joint(&knee_joint[2], SetVector(.4, .9, -.4));
	create_joint(&knee_joint[3], SetVector(.4, .8, 0));

	create_joint(&foot_joint[0], SetVector(-2.75, 0, .4));
	create_joint(&foot_joint[1], SetVector(-2.75, 0, -.4));
	create_joint(&foot_joint[2], SetVector(.34, 0, -.4));
	create_joint(&foot_joint[3], SetVector(.34, 0, .4));


	//SET PARENTS!
	foot_joint[0].parent = &knee_joint[0];
	foot_joint[1].parent = &knee_joint[1];
	foot_joint[2].parent = &knee_joint[2];
	foot_joint[3].parent = &knee_joint[3];

	knee_joint[0].parent = &legbase_joint[0];
	knee_joint[1].parent = &legbase_joint[1];
	knee_joint[2].parent = &legbase_joint[2];
	knee_joint[3].parent = &legbase_joint[3];

	
	legbase_joint[0].parent = NULL;
	legbase_joint[1].parent = NULL;
	legbase_joint[2].parent = NULL;
	legbase_joint[3].parent = NULL;


	head_joint[0].child = &legbase_joint[0];
	legbase_joint[0].child = &knee_joint[0];
	knee_joint[0].child = &foot_joint[0];
	foot_joint[0].child = NULL;



	//TAIL JOINTS
	create_joint(&tail_joint[0], SetVector(.7, 3.75, 0));
	create_joint(&tail_joint[1], SetVector(2.2, 3.8, 0));
	create_joint(&tail_joint[2], SetVector(3.0, 3.75, 0));
	create_joint(&tail_joint[3], SetVector(3.9, 3.65, 0));

	tail_joint[3].parent = &tail_joint[2];
	tail_joint[2].parent = &tail_joint[1];
	tail_joint[1].parent = &tail_joint[0];
	tail_joint[0].parent = NULL;


	tail_joint[0].child = &tail_joint[1];
	tail_joint[1].child = &tail_joint[2];
	tail_joint[2].child = &tail_joint[3];
	tail_joint[3].child = NULL;


	//HEAD JOINTS
	create_joint(&head_joint[0], SetVector(-2.9, 3.2, 0));
	create_joint(&head_joint[1], SetVector(-3.85, 4, 0));
	create_joint(&head_joint[2], SetVector(-4.7, 3, 0));

	
	legbase_joint[0].parent = &head_joint[0];
	legbase_joint[1].parent = &head_joint[0];
	legbase_joint[2].parent = &head_joint[0];
	legbase_joint[3].parent = &head_joint[0];
	
	head_joint[0].parent = NULL;
	//head_joint[0].parent = &head_joint[1];


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
