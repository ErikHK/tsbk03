// New version by Ingemar 2010
// Removed all dependencies of the Wild Magic (wml) library.
// Replaced it with VectorUtils2 (in source)
// Replaced old shader module with the simpler "ShaderUtils" unit.

// 2013: Adapted to VectorUtils3 and MicroGlut.

//#include <fftw3.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "common/objects.h"
//#include "SFML/Graphics.h"
#include <sys/time.h>
#include <time.h>
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

typedef struct COMPLEX_t
{
  double real;
  double imag;

} COMPLEX;

int FFT(int dir,int m,double *x,double *y);
int Powerof2(int n,int *m,int *twopm);

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
joint_s left_ear_joint[3];
joint_s right_ear_joint[3];

cow_s cow;
floor_s f;
ball_s ball;
wall_s wall;
ragdoll_s ragdoll;
farmer_s farmer;
plank_s p;
fence_s ff;

float cam_angle = 0;
float cam_dist = 8;

float mouse_x, mouse_y, old_mouse_x, old_mouse_y;
float m_angle;

mat4 modelViewMatrix, projectionMatrix;
mat4 cam_matrix;
Model * terr;
Model * skybox;
void draw_skybox(Model * skybox, GLuint program);


void filter(GLfloat * Xre, GLfloat * Xim, int size)
{
  int i,j;

  for(i=0;i < size;i++)
  {
    for(j=0;j < size;j++)
    {
      Xre[j + i*size] *= -(i+j)*.7/size;
      Xim[j + i*size] *= -(i+j)*.7/size;

    }

  }
}

float cosf(float x)
{
  return 1-x*x/2;
}

float sinf(float x)
{
  return x-x*x*x/6.;
}

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
      the dimensions are not powers of 2
*/
int FFT2D(COMPLEX **c,int nx,int ny,int dir)
{
   int i,j;
   int m,twopm;
   double *real,*imag;

   /* Transform the rows */
   real = (double *)malloc(nx * sizeof(double));
   imag = (double *)malloc(nx * sizeof(double));
   if (real == NULL || imag == NULL)
      return(0);
   if (!Powerof2(nx,&m,&twopm) || twopm != nx)
      return(0);
   for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
         real[i] = c[i][j].real;
         imag[i] = c[i][j].imag;
      }
      FFT(dir,m,real,imag);
      for (i=0;i<nx;i++) {
         c[i][j].real = real[i];
         c[i][j].imag = imag[i];
      }
   }
   free(real);
   free(imag);

   /* Transform the columns */
   real = (double *)malloc(ny * sizeof(double));
   imag = (double *)malloc(ny * sizeof(double));
   if (real == NULL || imag == NULL)
      return 0;
   if (!Powerof2(ny,&m,&twopm) || twopm != ny)
      return 0;
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         real[j] = c[i][j].real;
         imag[j] = c[i][j].imag;
      }
      FFT(dir,m,real,imag);
      for (j=0;j<ny;j++) {
         c[i][j].real = real[j];
         c[i][j].imag = imag[j];
      }
   }
   free(real);
   free(imag);

   return 1;
}


/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

     Formula: forward
                  N-1
                  ---
              1   \          - j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = forward transform
              N   /                                n=0..N-1
                  ---
                  k=0

      Formula: reverse
                  N-1
                  ---
                  \          j k 2 pi n / N
      X(n) =       >   x(k) e                    = forward transform
                  /                                n=0..N-1
                  ---
                  k=0
*/
int FFT(int dir,int m,double *x,double *y)
{
   long nn,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   nn = 1;
   for (i=0;i<m;i++)
      nn *= 2;

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<nn;i++) {
         x[i] /= (double)nn;
         y[i] /= (double)nn;
      }
   }

   return(1);
}

/*-------------------------------------------------------------------------
   Calculate the closest but lower power of two of a number
   twopm = 2**m <= n
   Return TRUE if 2**m == n
*/
int Powerof2(int n,int *m,int *twopm)
{
   if (n <= 1) {
      *m = 0;
      *twopm = 1;
      return(0);
   }

   *m = 1;
   *twopm = 2;
   do {
      (*m)++;
      (*twopm) *= 2;
   } while (2*(*twopm) <= n);

   if (*twopm != n)
      return(0);
   else
      return(1);
}



void fft(GLfloat * vertices, int size, GLfloat * Xre, GLfloat * Xim)
{
  int n,m,k,l,z=1;
  GLfloat xn=0;
  //GLfloat Xre[size][size];
  //GLfloat Xim[size][size];

  //Xre = malloc(sizeof(GLfloat)*size*size);
  //Xim = malloc(sizeof(GLfloat)*size*size);

  for(m=0;m < size; m++)
  {
  for(n=0;n < size; n++)
  {
    Xim[n + m*size] = 0;
    Xre[n + m*size] = 0;
    for(k=0;k < size; k++)
    {
    for(l=0;l < size; l++)
    {
      xn = vertices[(n + (m)*size)*3 + 1];
      Xre[l + size*k] += xn*cos(2*M_PI*(k*m+l*n)/size)/(size);
      Xim[l + size*k] -= xn*sin(2*M_PI*(k*m+l*n)/size)/(size);
    }
    }
  }
  }

}

void ifft(GLfloat * vertices, int size, GLfloat * Xre, GLfloat * Xim)
{
  int n,m,k,l,z=1;
  GLfloat xn=0, xn2=0;
  //GLfloat Xre[size][size];
  //GLfloat Xim[size][size];

  //Xre = malloc(sizeof(GLfloat)*size*size);
  //Xim = malloc(sizeof(GLfloat)*size*size);

  for(m=0;m < size; m++)
  {
  for(n=0;n < size; n++)
  {

    for(k=0;k < size; k++)
    {
    for(l=0;l < size; l++)
    {
      //xn = vertices[(n + (m)*size)*3 + 1];
      xn = Xre[n+size*m] + Xim[n + size*m];
      //xn = Xre[n + size*m];
      //xn2 = Xim[n + size*m];
      Xre[l + size*k] += xn*cos(2*M_PI*(k*m+l*n)/size)/size;
      Xim[l + size*k] += xn*sin(2*M_PI*(k*m+l*n)/size)/size;
      //vertices[(k + size*l)*3 + 1] = Xre[k + size*l] + Xim[k + size*l];

    }
    }

    vertices[(n + size*m)*3 + 1] = Xre[n + size*m] + Xim[n + size*m];

  }
  }

  
}


float random()
{
  return (float)rand() / (float)RAND_MAX;

}

void spatial_smooth(GLfloat * vertices, int size)
{
  int x,z;
  for(x=0;x < size; x++)
  {
    for(z = 0;z < size; z++)
    {
      if(x > 1 && z > 1 && x < size-1 && z < size-1)
        vertices[(x + z*size)*3 + 1] += vertices[(x+1 + z*size)*3 + 1] +
        vertices[(x + (z+1)*size)*3 + 1] + 
        vertices[(x+1 + (z-1)*size)*3 + 1] +
        vertices[(x + (z+1)*size)*3 + 1] +
        vertices[(x + (z-1)*size)*3 + 1] +
        vertices[(x-1 + (z-1)*size)*3 + 1] +
        vertices[(x-1 + (z+1)*size)*3 + 1] +
        vertices[(x+1 + (z+1)*size)*3 + 1];


        vertices[(x + (z)*size)*3 + 1] /= 9;
    }
  }
}


void calc_normal(GLfloat * vertexArray, int x, int z, int width, Point3D *normal)
{
        Point3D vec1, vec2;

        if(x > 0 && z > 0 && x < width && z < width)
        {
                vec1.x = vertexArray[(x-1 + z * width)*3 + 0] - 
                vertexArray[(x + z * width)*3 + 0];

                vec1.y = vertexArray[(x-1 + z * width)*3 + 1] - 
                vertexArray[(x + z * width)*3 + 1];

                vec1.z = vertexArray[(x-1 + z * width)*3 + 2] - 
                vertexArray[(x + z * width)*3 + 2];


                vec2.x = vertexArray[(x + (z+1) * width)*3 + 0] - 
                vertexArray[(x + z * width)*3 + 0];

                vec2.y = vertexArray[(x + (z+1) * width)*3 + 1] - 
                vertexArray[(x + z * width)*3 + 1];

                vec2.z = vertexArray[(x + (z+1) * width)*3 + 2] - 
                vertexArray[(x + z * width)*3 + 2];


                *normal = Normalize(CrossProduct(vec1, vec2));
      }

}



Model * generate_terrain(int size)
{
  int i, j, x, z;
  int tri_count = (size-1)*(size-1)*2;
  int count = size*size;
  Point3D tmp_normal;

  GLfloat *vertices = malloc(sizeof(GLfloat) * 3 * size * size);
  GLfloat *normals = malloc(sizeof(GLfloat) * 3 * size * size);
  GLuint *indices = malloc(sizeof(GLuint) * 3 * tri_count);

  COMPLEX * X = malloc(sizeof(GLfloat)*3*6*size*size);

  for(x=0;x < size; x++)
  {
    for(z = 0;z < size; z++)
    {
      vertices[(z + x*size)*3 + 0] = x*6;
      vertices[(z + x*size)*3 + 1] = 20*(random()-.5);
      vertices[(z + x*size)*3 + 2] = z*6;

      X[(z + x*size)*3 + 0].real = x*6;
      X[(z + x*size)*3 + 1].real = 20*(random()-.5);
      X[(z + x*size)*3 + 2].real = z*6;

      X[(z + x*size)*3 + 0].imag = 0;
      X[(z + x*size)*3 + 1].imag = 0;
      X[(z + x*size)*3 + 2].imag = 0;
    }
  }


  for(x=0;x < size; x++)
  {
    for(z = 0;z < size; z++)
    {
      calc_normal(vertices, x, z, size, &tmp_normal);

      normals[(z + x*size)*3 + 0] = tmp_normal.x;
      normals[(z + x*size)*3 + 1] = tmp_normal.y;
      normals[(z + x*size)*3 + 2] = tmp_normal.z;

    }

  }

  for(x=0;x < size-1; x++)
  {
    for(z = 0;z < size-1; z++)
    {

      indices[(x + z*(size-1))*6 + 0] = x + z*size;
      indices[(x + z*(size-1))*6 + 1] = x + (z+1)*size;
      indices[(x + z*(size-1))*6 + 2] = x + 1 + z*size;

      indices[(x + z*(size-1))*6 + 3] = x + 1 + z*size;
      indices[(x + z*(size-1))*6 + 4] = x + (z+1)*size;
      indices[(x + z*(size-1))*6 + 5] = x + 1 + (z+1)*size;
    }
  }

  //int j=0;
  //for(j=0;j < 30;j++)
  //  spatial_smooth(vertices, size);

  GLfloat Xre[size*size*6*6];
  GLfloat Xim[size*size*6*6];

  //COMPLEX X[size*size*6*6];

  //Xre = malloc(sizeof(GLfloat)*size*size);
  //Xim = malloc(sizeof(GLfloat)*size*size);

  //fft(vertices, size, Xre, Xim);

  //FFT2D(&X, size, size, 1);
  //filter(Xre, Xim, size);
  //ifft(vertices, size, Xre, Xim);

  //printf("%f\n", Xre[0][0]);

  Model * m = LoadDataToModel(
	vertices,
	normals,
	NULL,
	NULL,
	indices,
	size*size,
	(size-1)*(size-1)*2*3);

  return m;
}

void DisplayWindow()
{

	draw_skybox(skybox, g_shader);

	glUniformMatrix4fv(glGetUniformLocation(g_shader, "cam_matrix"), 1, GL_TRUE, cam_matrix.m);


	//glClearColor(0.4, 0.4, 0.8, 1);
	//glClear(GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT);

	glUniform1i(glGetUniformLocation(g_shader, "draw_cow"), 1);
	draw_cow(&cow, g_shader);
	glUniform1i(glGetUniformLocation(g_shader, "draw_cow"), 0);

	glUniform1i(glGetUniformLocation(g_shader, "draw_floor"), 1);
	draw_floor(&f, g_shader);
	glUniform1i(glGetUniformLocation(g_shader, "draw_floor"), 0);
	/*
	glUniform1i(glGetUniformLocation(g_shader, "draw_ball"), 1);
	draw_ball(&ball, g_shader);
	glUniform1i(glGetUniformLocation(g_shader, "draw_ball"), 0);
	draw_wall(&wall, g_shader);
	*/
	//draw_ragdoll(&ragdoll, g_shader);


	//glUniform1i(glGetUniformLocation(g_shader, "draw_farmer"), 1);
	//glEnable (GL_BLEND);
	//glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glDisable(GL_DEPTH_TEST);

	draw_farmer(&farmer, g_shader);
	//glEnable(GL_DEPTH_TEST);
	//glDisable(GL_BLEND);
	//glUniform1i(glGetUniformLocation(g_shader, "draw_farmer"), 0);
	glUniform1i(glGetUniformLocation(g_shader, "draw_plank"), 1);
	//draw_plank(&p, g_shader);
	draw_fence(&ff, g_shader);
	glUniform1i(glGetUniformLocation(g_shader, "draw_plank"), 0);


	//draw_debug_sphere(&ball, wall.bb.pos, g_shader);

	//draw_debug_sphere(&ball, cow.bb.pos, g_shader);
	draw_debug_sphere(&ball, cow.head_pos, g_shader);
	
	int i;
	for(i=0;i < 8;i++)
	{
	  draw_debug_sphere(&ball, cow.bb.vertices[i], g_shader);
	  //draw_debug_sphere(&ball, wall.bb.vertices[i], g_shader);
	}
	
	//draw_debug_sphere(&ball, VectorAdd(cow.bb.pos, cow.bb.size), g_shader);


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

void calc_bone_transform_cow(joint_s * j, int acc)
{
  joint_s * jc;
  joint_s * rootj = j;
  mat4 tmp, tmptrans, invtrans;
  //tmp = IdentityMatrix();
  if(acc)
    tmp = j->parent->tmp;
    //tmp = j->parent->Mtot;
  else
    //tmp = IdentityMatrix();
    tmp = Mult(T(cow.pos.x, cow.pos.y, cow.pos.z), Ry(cow.angle));

  GLfloat Ms[8][16];
  int i=0,ii=0, k=0;
  float currpos[8*3] = {0};
  float bonepos[8*3] = {0};

  while(j->child[0] != NULL)
  {
    //calc_bone_transform(j->child[0], 1);
    for(k=1; k < MAX_CHILDREN; k++)
    {
      if(j->child[k] != NULL)
        calc_bone_transform_cow(j->child[k], 1);
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

  if(rootj->posvar != NULL)
  {
  //printf(rootj->posvar);
  //printf("\n");
  }

  //printf("%f %f %f\n", currpos[0], currpos[1], currpos[2]);
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

	GLfloat tmpppp[3] = {cow.pos.x, cow.pos.y, cow.pos.z};
	glUniform3fv(glGetUniformLocation(g_shader, "cow_pos"), 1, tmpppp);

	vec3 tmpz = SetVector(1,1,1);
	tmpz = MultVec3(Ry(cow.angle), tmpz);
	GLfloat tmppppp[3] = {tmpz.x, tmpz.y, tmpz.z};
	glUniform3fv(glGetUniformLocation(g_shader, "light_normal"), 1, tmppppp);


	turn_cow(&cow, -m_angle);
	update_cow(&cow, delta_t);
	move_cow(&cow, m_angle);
	update_floor(&f, &cow);
	//update_wall(&wall, &cow, delta_t);
	//update_ball(&ball, &cow, delta_t);
	update_farmer(&farmer, t);
	update_fence(&ff, &cow, delta_t);


	//check collision between cow and farmer
	if(((check_sphere_collision(cow.head_pos, 
	farmer.skeleton.joints[8].pos, 1, 1))
	|| (check_sphere_collision(cow.head_pos, 
	farmer.skeleton.joints[0].pos, 1.5, 1.5))) && farmer.animate)
	{
	  farmer.skeleton.joints[0].leader = 1;
	  farmer.skeleton.joints[0].speed = ScalarMult(cow.speed, 1.5);
	  cow.momentum = ScalarMult(cow.momentum, .6);
	  //cow.pos = VectorAdd(cow.pos, ScalarMult(Normalize(cow.momentum), -.8));
	  glUniform1i(glGetUniformLocation(g_shader, "collision"), 1);
	  //switch to ragdoll mode!
	  farmer.animate = 0;
          reset_farmer_matrices(&farmer);
	  //printf("JAPP\n");
	}
	else
	{
	  glUniform1i(glGetUniformLocation(g_shader, "collision"), 0);
	  //printf("NOPP\n");
	}


        if(!farmer.animate)
        {
	  //update_ragdoll(&ragdoll, delta_t);
	  //if(keyIsDown('k'))
            update_ragdoll(&farmer.skeleton, delta_t);
        }



	if(check_collision_2(&cow.bb, &wall.bb))
	{
	  //printf("yey\n");
	  glUniform1i(glGetUniformLocation(g_shader, "collision"), 1);
	  //cow.momentum = SetVector(-cow.momentum.x, cow.momentum.y, -cow.momentum.z);
	  //wall_collision(wall, cow);
	}
	//else
	//  glUniform1i(glGetUniformLocation(g_shader, "collision"), 0);



	//printf("%f, %f\n", delta_t, old_t);
	glUniform1f(glGetUniformLocation(g_shader, "time"), t);


	mat4 proj_matrix = frustum(-1, 1, -1, 1, 1, 750.0);
	cam_matrix = lookAt(cam_dist*cos(m_angle)+cow.pos.x, 9+cow.pos.y,  cam_dist*sin(m_angle)+cow.pos.z, cow.pos.x, 8.5+cow.pos.y, cow.pos.z, 0.0, 1.0, 0.0);
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


	float freq = 7;//Norm(cow.speed)/4;

	if(Norm(cow.momentum) < .5)
	  freq = 0;
	else if(Norm(cow.momentum) > 5)
	  freq = 20;
	//printf("%f\n", freq);

	j = &thigh_joint[0];
	jc = j->child[0];
	j->R = ArbRotate(SetVector(0,0,1), sin(freq*t)/1.5);
	jc->R = ArbRotate(SetVector(0,0,1), cos(freq*t)/2.5);

	j = &thigh_joint[1];
	jc = j->child[0];
	j->R = ArbRotate(SetVector(0,0,1), -cos(freq*t)/2.5);
	jc->R = ArbRotate(SetVector(0,0,1), sin(freq*t)/2.5);


	j = &thigh_joint[2];
	jc = j->child[0];
	j->R = ArbRotate(SetVector(0,0,1), sin(freq*t)/1.5);
	jc->R = ArbRotate(SetVector(0,0,1), cos(freq*t)/2.5);

	j = &thigh_joint[3];
	jc = j->child[0];
	j->R = ArbRotate(SetVector(0,0,1), -cos(freq*t)/2.5);
	jc->R = ArbRotate(SetVector(0,0,1), sin(freq*t)/2.5);


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
	calc_bone_transform_cow(&legbase_joint[0], 0);
	calc_bone_transform_cow(&legbase_joint[1], 0);
	calc_bone_transform_cow(&legbase_joint[2], 0);
	calc_bone_transform_cow(&legbase_joint[3], 0);

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
	j->R = ArbRotate(SetVector(0,1,0), m_angle + cow.angle);
	j->R = Mult(j->R, ArbRotate(SetVector(0,0,1), cos(freq*t)/9));

	calc_bone_transform_cow(&body_joint[0],0);

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

	calc_bone_transform_cow(&head_joint[0],0);

	glutPostRedisplay();
}

void keyboardFunc( unsigned char key, int x, int y)
{
// Add any keyboard control you want here
	if(key == 27)	//Esc
		exit(-1);
}


void draw_skybox(Model * skybox, GLuint program)
{

  glDisable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glUniform1i(glGetUniformLocation(program, "draw_skybox"), 1);

  //curr_trans = &trans1;
  //glUniformMatrix4fv(glGetUniformLocation(program, "myMatrix"), 1, GL_TRUE, myMatrix);
  //glUniformMatrix4fv(glGetUniformLocation(program, "camMatrix"), 1, GL_TRUE, look.m);

  GLfloat skyboxmat[16];
  memcpy(skyboxmat, cam_matrix.m, sizeof(skyboxmat));
  skyboxmat[3] = 0;
  skyboxmat[7] = 0;
  skyboxmat[11] = 0;
  skyboxmat[15] = 1;

  //glUniformMatrix4fv(glGetUniformLocation(program, "projMatrix"), 1, GL_TRUE, skyboxmat);
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, IdentityMatrix().m);
  glUniformMatrix4fv(glGetUniformLocation(program, "cam_matrix"), 1, GL_TRUE, skyboxmat);
  //glBindTexture(GL_TEXTURE_2D, tex1);
  DrawModel(skybox, program, "inPosition", "inNormal", "inTexCoord");

  glEnable(GL_DEPTH_TEST);
  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glUniform1i(glGetUniformLocation(program, "draw_skybox"), 0);
}


void mouse(int x, int y)
{

  old_mouse_x = mouse_x;
  old_mouse_y = mouse_y;
  if(x > 700 || x < 150)
    glutWarpPointer(400,y);
  if(y > 500 || y < 200)
    glutWarpPointer(x,300);


  mouse_x = x;
  mouse_y = y;

  //if(m_angle > -.8)
 // {
  if(mouse_x - old_mouse_x > 0)
  {
    m_angle -= .03;
  }
  //}

  //if(m_angle < .8)
  //{
  if(mouse_x - old_mouse_x < 0)
  {
    m_angle += .03;
  }
  //}
  


}


/////////////////////////////////////////
//		M A I N
//
int main(int argc, char **argv)
{
	srand(time(NULL));
	mouse_x = 0;
	old_mouse_x = 0;
	mouse_y = 0;
	old_mouse_y = 0;

	m_angle = -8*M_PI;


	glutInit(&argc, argv);


	initKeymapManager();
	glutPassiveMotionFunc(mouse);

	//glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GL_MULTISAMPLE);
	glutInitWindowSize(800, 600);
	glutInitContextVersion(3, 2); // Might not be needed in Linux
	glutCreateWindow("Farm Escape");
	glutDisplayFunc(DisplayWindow);

	create_floor(&f);
	create_fence(&ff, 14, SetVector(0,0,30));
	//f.model = generate_terrain(32);
	double * X = (double *)malloc(sizeof(double)*32);
	double * Y = (double *)malloc(sizeof(double)*32);
	int ii;
	for(ii = 0;ii < 32;ii++)
	{
	  X[ii] = 0;
	  Y[ii] = 0;
	}
	FFT(1, 5, X, Y);
	create_ball(&ball, SetVector(5,0,0));
	create_wall(&wall, SetVector(0,10,0), SetVector(2,5,2));

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

	//EAR JOINTS
	//create_joint(&left_ear_joint[0], SetVector(-2.9, 3.2, -.4),
	//"leftearjoint", "leftearcurrpos", "leftearbonepos", 0);
	create_joint(&left_ear_joint[0], SetVector(-2.9, 3.2, .6),
	NULL, NULL, NULL, 0);
	create_joint(&right_ear_joint[0], SetVector(-2.9, 3.2, -.6),
	NULL, NULL, NULL, 0);
	//create_joint(&left_ear_joint[1], SetVector(-2.9, 3.2, .6),
	//NULL, NULL, NULL, 0);
	//create_joint(&right_ear_joint[1], SetVector(-2.9, 3.2, -.6),
	//NULL, NULL, NULL, 0);

	//create_joint(&right_ear_joint[0], SetVector(-2.9, 3.2, -.6),
	//"rightearjoint", "rightearcurrpos", "rightearbonepos", 0);

	//create_joint(&right_ear_joint[1], SetVector(-2.9, 3.2, -.8),
	//"rightearjoint", "rightearcurrpos", "rightearbonepos", 0);


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
	//head_joint[2].child[0] = &right_ear_joint[1];
	head_joint[2].child[1] = &left_ear_joint[0];
	head_joint[2].child[2] = &right_ear_joint[0];

	//left_ear_joint[0].child[0] = &left_ear_joint[1];
	//right_ear_joint[0].child[0] = &right_ear_joint[1];

	//left_ear_joint[1].child[0] = NULL;
	//right_ear_joint[1].child[0] = NULL;

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

	//left_ear_joint[1].parent = &left_ear_joint[0];
	left_ear_joint[0].parent = &head_joint[2];

	//right_ear_joint[1].parent = &right_ear_joint[0];
	right_ear_joint[0].parent = &head_joint[2];


	create_cow(&cow);
	create_farmer(&farmer, SetVector(0,0,0));
	skybox = LoadModelPlus("./res/sphere.obj");

	//create_ragdoll(&ragdoll);

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

	//glutSetCursor(0);

	//terr = generate_terrain(512);


	glutTimerFunc(20, &OnTimer, 0);

	glutMainLoop();
	exit(0);
}
