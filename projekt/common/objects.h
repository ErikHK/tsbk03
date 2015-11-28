#include "loadobj.h"
#include "VectorUtils3.h"
#include "GL_utilities.h"
#include "LoadTGA.h"

#define MAX_CHILDREN		4

#define COW_FORCE		140
#define COW_MAX_MOMENTUM	60
#define FLOOR_FRICTION		80
#define COW_GRAVITY		-80



typedef struct bounding_box_t
{
  vec3 pos;
  vec3 size;
  vec3 vertices[8];
  vec3 center;
} bounding_box_s;



typedef struct cow_t
{
  Model * main_body;
  Model * udder;
  Model * debug_sphere;
  mat4 matrix;

  vec3 pos;
  vec3 head_pos;
  vec3 speed;
  vec3 acc;
  vec3 force;
  vec3 torque;
  vec3 momentum;
  mat4 R;
  float mass;
  float angle;
  vec3 omega;
  float d_angle;
  int jumping;
  GLuint tex;
  struct bounding_box_t bb;

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
  vec3 speed;
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
  float dist_to_parent;
  float dist_to_child[MAX_CHILDREN];
  vec3 force;
  vec3 calculated_force;
  char * Mvar;
  char * posvar;
  char * boneposvar;
  int isnull;

} joint_s;

typedef struct skeleton_t
{
  joint_s joints[20];
  float distances[20];
} skeleton_s;

typedef struct farmer_t
{
  Model * body;
  mat4 matrix;
  struct skeleton_t skeleton;

  vec3 pos;
  vec3 speed;
  vec3 acc;
  vec3 force;
  vec3 torque;
  vec3 momentum;
  mat4 R;
  float mass;
  vec3 omega;
  GLuint tex;

} farmer_s;


typedef struct ragdoll_t
{
  //Model * joint_model;
  joint_s joints[20];

} ragdoll_s;

typedef struct wall_t
{
  Model * wall_model;
  vec3 pos;
  vec3 size;
  mat4 T;
  mat4 R;
  mat4 matrix;
  mat4 orig_matrix;
  vec3 speed;
  vec3 acc;
  vec3 force;
  vec3 torque;
  vec3 momentum;
  vec3 angular_momentum;
  vec3 omega;
  float mass;

  struct bounding_box_t bb;

} wall_s;

typedef struct plank_t
{
  Model * body;
  Model * debug_sphere;
  mat4 R;
  mat4 T;
  vec3 pos;
  vec3 size;
  int type;

  float mass;

} plank_s;

typedef struct fence_t
{
  plank_s planks[400];
  vec3 pos;
  int width;

} fence_s;


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

void create_floor(floor_s * f);

void create_ball(ball_s * b, vec3 pos);
void update_ball(ball_s * b, cow_s * c, GLfloat dT);
void update_wall(wall_s * w, cow_s * c, GLfloat dT);
void draw_ball(ball_s * b, GLuint program);
void create_bb(bounding_box_s * bb, vec3 pos, vec3 size);
void create_wall();
void update_vertices(bounding_box_s * bb, vec3 pos, mat4 R);
int check_collision(bounding_box_s * b1, bounding_box_s * b2);
int check_collision_2(bounding_box_s * b1, bounding_box_s * b2);
int sign(float x);
void draw_debug_sphere(ball_s * b, vec3 pos, GLuint program);
void draw_floor(floor_s * f, GLuint program);
void create_ragdoll(ragdoll_s * r);
void create_ragdoll_joint(joint_s * j, vec3 pos);
void update_ragdoll(ragdoll_s * r, GLfloat dT);
void draw_ragdoll(ragdoll_s * r, GLuint program);

void create_farmer(farmer_s * f, vec3 pos);
void draw_farmer(farmer_s * f, GLuint program);
void update_farmer(farmer_s * f);

void create_plank(plank_s * p, vec3 pos, int type);
void draw_plank(plank_s * p, GLuint program);

void create_fence(fence_s * f, int width, vec3 pos);
void draw_fence(fence_s * f, GLuint program);
//void update_fence(fence_s * f, cow_s * c, GLfloat dT);
