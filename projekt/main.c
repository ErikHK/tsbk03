// New version by Ingemar 2010
// Removed all dependencies of the Wild Magic (wml) library.
// Replaced it with VectorUtils2 (in source)
// Replaced old shader module with the simpler "ShaderUtils" unit.

// 2013: Adapted to VectorUtils3 and MicroGlut.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "common/objects.h"
#include "SFML/Graphics.h"
#include <sys/time.h>
#include "zpr.h"
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

GLfloat delta_t, old_t, t, current_time;

static double startTime = 0;

bone_s bone;
joint_s foot_joint[4];
joint_s knee_joint[4];
joint_s thigh_joint[4];
joint_s legbase_joint[4];
joint_s tail_joint[4];
joint_s head_joint[3];
joint_s body_joint[3];

cow_s cow;
floor_s f;
float cam_angle = 0;
float cam_dist = 8;

float mouse_x, mouse_y, old_mouse_x, old_mouse_y;
float m_angle;

mat4 modelViewMatrix, projectionMatrix;

void DisplayWindow()
{
	glClearColor(0.4, 0.4, 0.2, 1);
	glClear(GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT);

	glUniform1i(glGetUniformLocation(g_shader, "draw_cow"), 1);
	draw_cow(&cow, g_shader);
	glUniform1i(glGetUniformLocation(g_shader, "draw_cow"), 0);
	draw_floor(&f, g_shader);

/*
	int i;
	for(i=0; i<4;i++)
	{
	  draw_joint(&foot_joint[i], g_shader);
	  draw_joint(&knee_joint[i], g_shader);
	  draw_joint(&thigh_joint[i], g_shader);
	  draw_joint(&legbase_joint[i], g_shader);
	  draw_joint(&tail_joint[i], g_shader);
	  if(i < 3) 
	  {
	    draw_joint(&head_joint[i], g_shader);
	    draw_joint(&body_joint[i], g_shader);
	  }

	}
*/
	//draw_bone(&bone, g_shader);
	glutSwapBuffers();
	//printf("runrunrun\n");
}

void calc_bone_transform(joint_s * j, int acc)
{
  joint_s * jc;
  joint_s * rootj = j;
  mat4 tmp, tmptrans, invtrans;
  //tmp = IdentityMatrix();
  if(acc)
    tmp = j->parent->tmp;
  else
    tmp = IdentityMatrix();

  GLfloat Ms[8][16];
  int i=0,ii=0, k=0;
  float currpos[8*3] = {0};
  float bonepos[8*3] = {0};

  while(j->child[0] != NULL)
  {
    for(k=1; k<MAX_CHILDREN; k++)
    {
      if(j->child[k] != NULL)
        calc_bone_transform(j->child[k], 1);
    }
    jc = j->child[0];
    vec3 tmp_bonepos;
    tmptrans = j->T;
    tmp = Mult(tmp, tmptrans);
    invtrans = InvertMat4(tmptrans);
    tmp = Mult(tmp, Mult(j->R, invtrans));
    j->tmp = tmp;
    j->isnull = 0;

    //middle of bone
    tmp_bonepos = ScalarMult(
		VectorAdd(j->pos, jc->pos), .5);

    for(ii=0;ii<16;ii++)
      Ms[i][ii] = (tmp).m[ii];

    currpos[i*3] = 10*j->pos.x;
    currpos[i*3+1] = 10*j->pos.y;
    currpos[i*3+2] = 10*j->pos.z;

    bonepos[i*3] = 10*tmp_bonepos.x;
    bonepos[i*3+1] = 10*tmp_bonepos.y;
    bonepos[i*3+2] = 10*tmp_bonepos.z;

    j = j->child[0];
    i++;
  }

  glUniformMatrix4fv(glGetUniformLocation(g_shader, rootj->Mvar), 8, GL_TRUE, Ms[0]);
  glUniform3fv(glGetUniformLocation(g_shader, rootj->posvar), 8, currpos);
  glUniform3fv(glGetUniformLocation(g_shader, rootj->boneposvar), 8, bonepos);



}


void OnTimer(int value)
{
	glutTimerFunc(20, &OnTimer, value);

	old_t = t;
	t = (GLfloat)glutGet(GLUT_ELAPSED_TIME)/1000.0;
	delta_t = (t - old_t);
	//printf("%f\n", delta_t);

	turn_cow(&cow, -m_angle);
	move_cow(&cow, m_angle);
	update_cow(&cow, delta_t);
	update_floor(&f, &cow);

	//printf("%f, %f\n", delta_t, old_t);
	glUniform1f(glGetUniformLocation(g_shader, "time"), t);


	mat4 proj_matrix = frustum(-1, 1, -1, 1, 1, 750.0);
	mat4 cam_matrix = lookAt(cam_dist*cos(m_angle), 9,  cam_dist*sin(m_angle), 0, 8.5, 0, 0.0, 1.0, 0.0);
	//mat4 cam_matrix = lookAt(200+200*cos(cam_angle), 0,  200, 0, 8, 0, 0.0, 1.0, 0.0);

	//g_shader = loadShaders("shader.vert" , "shader.frag");
	//glUseProgram(g_shader);

	glUniformMatrix4fv(glGetUniformLocation(g_shader, "proj_matrix"), 1, GL_TRUE, proj_matrix.m);
	glUniformMatrix4fv(glGetUniformLocation(g_shader, "cam_matrix"), 1, GL_TRUE, cam_matrix.m);
	//glutDisplayFunc(DisplayWindow);

	//draw_cow(&cow, g_shader);
	//draw_bone(&bone, g_shader);

	



	mat4 tmp, testjoint, testjoint2;
	testjoint = IdentityMatrix();

	joint_s * j = &head_joint[0];
	joint_s * jc = j->child[0];
	joint_s * jcc = jc->child[0];
	//joint_s * jp = j->parent;
	//joint_s * jpp = jp->parent;

	float Ms[8][4*4];
	float legMs[8][4*4];
	float currpos[8*3] = {0};
	float bonepos[8*3] = {0};
	float legcurrpos[8*3] = {0};
	float legbonepos[8*3] = {0};
	//GLfloat * Mtmp = Ms;
	int i=0, ii=0;
	//jc->R = ArbRotate(SetVector(0,0,1), cos(4*t/(i+1))/2.5);
	//j->R = ArbRotate(SetVector(0,0,1), 0);
	j->R = ArbRotate(SetVector(0,0,1), sin(4*t/(3+1))/1.2);
	//jcc->R = ArbRotate(SetVector(0,0,1), cos(8*t/(2+1))/4);
	mat4 Mpacc, Minvacc, tmptrans, tmppp, invtrans;
	tmppp = IdentityMatrix();


	j = &thigh_joint[0];
	jc = j->child[0];
	j->R = ArbRotate(SetVector(0,0,1), sin(7*t/(i+1))/1.5);
	jc->R = ArbRotate(SetVector(0,0,1), cos(7*t/(i+1))/2.5);

	j = &thigh_joint[1];
	jc = j->child[0];
	j->R = ArbRotate(SetVector(0,0,1), -cos(7*t/(i+1))/2.5);
	jc->R = ArbRotate(SetVector(0,0,1), sin(7*t/(i+1))/2.5);


	j = &thigh_joint[2];
	jc = j->child[0];
	j->R = ArbRotate(SetVector(0,0,1), sin(7*t/(i+1))/1.5);
	jc->R = ArbRotate(SetVector(0,0,1), cos(7*t/(i+1))/2.5);

	j = &thigh_joint[3];
	jc = j->child[0];
	j->R = ArbRotate(SetVector(0,0,1), -cos(7*t/(i+1))/2.5);
	jc->R = ArbRotate(SetVector(0,0,1), sin(7*t/(i+1))/2.5);


	j = &legbase_joint[0];
	//j->R = ArbRotate(SetVector(1,0,0), -cos(7*t/(i+1)));
	//j->R = ArbRotate(SetVector(1,0,0), -M_PI/6);

	j = &thigh_joint[0];
	//j->R = Mult(ArbRotate(SetVector(1,0,0), cos(7*t/(i+1))), j->R);
	//j->R = Mult(ArbRotate(SetVector(1,0,0), M_PI/6), j->R);

	i = 0;
	tmppp = IdentityMatrix();
	//legbase_joint[0].R = ArbRotate(SetVector(0,0,1), cos(3*t)/2);

	//legs
	calc_bone_transform(&legbase_joint[0], 0);
	calc_bone_transform(&legbase_joint[1], 0);
	calc_bone_transform(&legbase_joint[2], 0);
	calc_bone_transform(&legbase_joint[3], 0);

	j = &tail_joint[0];
	j->R = ArbRotate(SetVector(0,0,1), -M_PI/2.5);
	j = &tail_joint[1];
	j->R = ArbRotate(SetVector(0,0,1), -M_PI/8.5);

	//j = &tail_joint[1];
	//j->R = ArbRotate(SetVector(0,0,1), 0);
	j = &tail_joint[2];
	//j->R = ArbRotate(SetVector(0,0,1), 0);
	j = &tail_joint[3];
	j->R = ArbRotate(SetVector(0,0,1), 0);

	//body
	j = &body_joint[1];
	//j->R = ArbRotate(SetVector(0,0,1), cos(t*7)/9.5);
	j->R = ArbRotate(SetVector(0,1,0), m_angle - 0);
	j->R = Mult(j->R, ArbRotate(SetVector(0,0,1), cos(t*7)/9));

	calc_bone_transform(&body_joint[0],0);

	//tail
	//j = &tail_joint[1];
	//j->R = ArbRotate(SetVector(0,0,1), -M_PI/2.5);

	//j = &tail_joint[1];
	//j->R = ArbRotate(SetVector(0,0,1), 0);

	//j = &tail_joint[2];
	//j->R = ArbRotate(SetVector(0,0,1), 0);

	//j = &tail_joint[3];
	//j->R = ArbRotate(SetVector(0,0,1), 0);

	//calc_bone_transform(&tail_joint[0],0);

	//head
	j = &head_joint[0];
	j->R = ArbRotate(SetVector(0,0,1), sin(7*t)/9.5);

	calc_bone_transform(&head_joint[0],0);


	//glUniformMatrix4fv(glGetUniformLocation(g_shader, "testjoint"), 8, GL_TRUE, Ms);
	//glUniform3fv(glGetUniformLocation(g_shader, "currpos"), 8, currpos);
	//glUniform3fv(glGetUniformLocation(g_shader, "bonepos"), 8, bonepos);

	//glUniformMatrix4fv(glGetUniformLocation(g_shader, "legtestjoint"), 8, GL_TRUE, legMs);
	//glUniform3fv(glGetUniformLocation(g_shader, "legcurrpos"), 8, legcurrpos);
	//glUniform3fv(glGetUniformLocation(g_shader, "legbonepos"), 8, legbonepos);


	

/*
	if(keyIsDown('e'))
	  cam_angle += .05;
	if(keyIsDown('i'))
	  cam_angle -= .05;

	if(keyIsDown('u'))
	  cam_dist += .05;

	if(keyIsDown('p'))
	  cam_dist -= .05;
*/



	glutPostRedisplay();
}

void keyboardFunc( unsigned char key, int x, int y)
{
// Add any keyboard control you want here
	if(key == 27)	//Esc
		exit(-1);
}

void mouse(int x, int y)
{

  old_mouse_x = mouse_x;
  old_mouse_y = mouse_y;
  if(x > 700 || x < 100)
    glutWarpPointer(400,y);
  if(y > 500 || y < 100)
    glutWarpPointer(x,300);


  mouse_x = x;
  mouse_y = y;

  if(m_angle > -.8)
  {
  if(mouse_x - old_mouse_x > 0)
  {
    m_angle -= .05;
  }
  }

  if(m_angle < .8)
  {
  if(mouse_x - old_mouse_x < 0)
  {
    m_angle += .05;
  }
  }
  


}


/////////////////////////////////////////
//		M A I N
//
int main(int argc, char **argv)
{

	mouse_x = 0;
	old_mouse_x = 0;
	mouse_y = 0;
	old_mouse_y = 0;

	m_angle = 0;

	//glutSetCursor(GLUT_CURSOR_NONE);


	glutInit(&argc, argv);
	initKeymapManager();
	glutPassiveMotionFunc(mouse);

	//glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(800, 600);
	glutInitContextVersion(3, 2); // Might not be needed in Linux
	glutCreateWindow("Farm Escape");
	glutDisplayFunc(DisplayWindow);

	create_floor(&f);

	create_joint(&legbase_joint[0], SetVector(-2.2, 3.8, .7), 
	"legjoint0", "legcurrpos0", "legbonepos0", 0);
	create_joint(&legbase_joint[1], SetVector(-2.2, 3.8, -.7),
	"legjoint1", "legcurrpos1", "legbonepos1", 0);
	create_joint(&legbase_joint[2], SetVector(.2, 3.9, -.7),
	"legjoint2", "legcurrpos2", "legbonepos2", 0);
	create_joint(&legbase_joint[3], SetVector(.2, 3.9, .7),
	"legjoint3", "legcurrpos3", "legbonepos3", 0);

	create_joint(&thigh_joint[0], SetVector(-2.3, 1.8, .4),
	NULL, NULL, NULL, 0);
	create_joint(&thigh_joint[1], SetVector(-2.3, 1.8, -.4),
	NULL, NULL, NULL, 0);
	create_joint(&thigh_joint[2], SetVector(.1, 2, -.4),
	NULL, NULL, NULL, 0);
	create_joint(&thigh_joint[3], SetVector(.1, 2.0, .4),
	NULL, NULL, NULL, 0);

	create_joint(&knee_joint[0], SetVector(-2.6, .9, .4),
	NULL, NULL, NULL, 0);
	create_joint(&knee_joint[1], SetVector(-2.6, .9, -.4),
	NULL, NULL, NULL, 0);
	create_joint(&knee_joint[2], SetVector(.4, .9, -.4),
	NULL, NULL, NULL, 0);
	create_joint(&knee_joint[3], SetVector(.4, .9, .4),
	NULL, NULL, NULL, 0);
	create_joint(&foot_joint[0], SetVector(-2.75, 0, .4),
	NULL, NULL, NULL, 0);
	create_joint(&foot_joint[1], SetVector(-2.75, 0, -.4),
	NULL, NULL, NULL, 0);
	create_joint(&foot_joint[2], SetVector(.34, 0, -.4),
	NULL, NULL, NULL, 0);
	create_joint(&foot_joint[3], SetVector(.34, 0, .4),
	NULL, NULL, NULL, 0);


	//BODY JOINTS
	create_joint(&body_joint[0], SetVector(-1.8, 3.6, 0),
	"bodyjoint", "bodycurrpos", "bodybonepos", 0);
	create_joint(&body_joint[1], SetVector(-1, 3.3, 0),
	NULL, NULL, NULL, 0);
	create_joint(&body_joint[2], SetVector(.14, 3.9, 0),
	NULL, NULL, NULL, 0);


	//TAIL JOINTS
	create_joint(&tail_joint[0], SetVector(.7+.3, 3.75, 0),
	"tailjoint", "tailcurrpos", "tailbonepos", 0);
	create_joint(&tail_joint[1], SetVector(2.2, 3.8, 0),
	NULL, NULL, NULL, 0);
	create_joint(&tail_joint[2], SetVector(3.0, 3.75, 0),
	NULL, NULL, NULL, 0);
	create_joint(&tail_joint[3], SetVector(3.9, 3.65, 0),
	NULL, NULL, NULL, 0);

	//HEAD JOINTS
	create_joint(&head_joint[0], SetVector(-2.9, 3.2, 0),
	"headjoint", "headcurrpos", "headbonepos", 0);
	create_joint(&head_joint[1], SetVector(-3.85, 4, 0),
	NULL, NULL, NULL, 0);
	create_joint(&head_joint[2], SetVector(-4.7, 3, 0),
	NULL, NULL, NULL, 0);


	//SET CHILDREN!
	body_joint[0].child[0] = &body_joint[1];
	body_joint[1].child[0] = &body_joint[2];
	body_joint[2].child[0] = &tail_joint[0];
	body_joint[2].child[1] = &legbase_joint[2];
	body_joint[2].child[2] = &legbase_joint[3];


	tail_joint[0].child[0] = &tail_joint[1];
	tail_joint[1].child[0] = &tail_joint[2];
	tail_joint[2].child[0] = &tail_joint[3];
	tail_joint[3].child[0] = NULL;


	head_joint[0].child[0] = &head_joint[1];
	head_joint[1].child[0] = &head_joint[2];
	head_joint[2].child[0] = NULL;

	//head_joint[0].child = &legbase_joint[0];
	int i;
        for(i=0;i<4;i++)
	{
	  legbase_joint[i].child[0] = &thigh_joint[i];
	  thigh_joint[i].child[0] = &knee_joint[i];
	  knee_joint[i].child[0] = &foot_joint[i];
	  foot_joint[i].child[0] = NULL;
	}

	//SET PARENTS
	legbase_joint[2].parent = &body_joint[2];
	legbase_joint[3].parent = &body_joint[2];


	create_cow(&cow);

	g_shader = loadShaders("shader.vert" , "shader.frag");
	glUseProgram(g_shader);

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
