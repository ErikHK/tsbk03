#include "loadobj.h"
#include "VectorUtils3.h"
#include "GL_utilities.h"
#include "LoadTGA.h"

#define MAX_CHILDREN		4

#define COW_FORCE		120
#define COW_MAX_MOMENTUM	60
#define FLOOR_FRICTION		60

typedef struct cow_t
{
  Model * main_body;
  Model * udder;
  mat4 matrix;

  vec3 pos;
  vec3 speed;
  vec3 acc;
  vec3 force;
  vec3 torque;
  vec3 momentum;
  float mass;
  float angle;
  float d_angle;
  int jumping;
  GLuint body_tex;

} cow_s;

typedef struct ball_t
{
  Model * model;
  mat4 matrix;
  vec3 pos;
  vec3 speed;
  vec3 acc;
  vec3 force;
  vec3 torque;
  vec3 momentum;
  float mass;
  GLuint tex;

} ball_s;


typedef struct floor_t
{
  Model * model;
  mat4 matrix;
  GLuint tex;
} floor_s;

typedef struct bone_t
{
  Model * body;
  vec3 orig_pos;
  vec3 pos;
  mat4 matrix;

} bone_s;

typedef struct joint_t
{
  Model * body;
  vec3 orig_pos;
  vec3 pos;
  mat4 body_matrix;
  mat4 T;
  mat4 R;
  mat4 M;
  mat4 Mp;
  mat4 Minv;
  mat4 Mtot;
  mat4 tmp;
  int id;
  struct joint_t *parent;
  struct joint_t *child[MAX_CHILDREN];
  char * Mvar;
  char * posvar;
  char * boneposvar;
  int isnull;

} joint_s;


void create_cow();
void draw_cow();
void turn_cow(cow_s * c, float angle);
void update_cow(cow_s * c, GLfloat dT);
void update_floor(floor_s * f, cow_s * c);
void move_cow(cow_s * c, float angle);
void create_bone();
void draw_bone();
void create_joint();
void draw_joint();

void create_ball(ball_s * b, vec3 pos);
void update_ball(ball_s * b, cow_s * c, GLfloat dT);
void draw_ball(ball_s * b, GLuint program);