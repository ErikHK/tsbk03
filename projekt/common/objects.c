#include "objects.h"
#include <stdlib.h>

int check_collision(bounding_box_s * b1, bounding_box_s * b2)
{
  if(abs(b1->pos.x - b2->pos.x) < b1->size.x + b2->size.x)
  {
    if(abs(b1->pos.y - b2->pos.y) < b1->size.y + b2->size.y)
    {
      if(abs(b1->pos.z - b2->pos.z) < b1->size.z + b2->size.z)
      {
        return 1;
      }
    }
  }
  return 0;
}

int sign(float x)
{
  if(x > 0)
    return 1;
  else if(x <= 0)
    return -1;
  return 0;
}

int check_collision_2(bounding_box_s * b1, bounding_box_s * b2)
{
  //check 6 different planes on one box, against 8 points on the other
  //6 normals to the planes
  vec3 n[6];
  //6 points on the planes
  vec3 p[6];
  vec3 tmp1, tmp2;
  //projection vectors for SAT
  vec3 sat[3];

  tmp1 = VectorSub(b1->vertices[4], b1->vertices[0]); //size.x
  tmp2 = VectorSub(b1->vertices[2], b1->vertices[0]); //size.y
  sat[0] = tmp1;
  sat[1] = tmp2;
  n[0] = CrossProduct(tmp2, tmp1); //-z
  p[0] = b1->vertices[0];

  tmp1 =  VectorSub(b1->vertices[1], b1->vertices[0]); //size.z
  sat[2] = tmp1;
  tmp2 =  VectorSub(b1->vertices[2], b1->vertices[0]); //size.y
  n[1] = CrossProduct(tmp1, tmp2);  //-x
  p[1] = b1->vertices[0];

  tmp1 = VectorSub(b1->vertices[4], b1->vertices[0]); //size.x
  tmp2 = VectorSub(b1->vertices[1], b1->vertices[0]); //size.z
  n[2] = CrossProduct(tmp1, tmp2); //-y
  p[2] = b1->vertices[0];


  n[3] = ScalarMult(n[0], -1); //z
  p[3] = b1->vertices[7];

  n[4] = ScalarMult(n[1], -1); //x
  p[4] = b1->vertices[7];

  n[5] = ScalarMult(n[2], -1); //y
  p[5] = b1->vertices[7];

  //distance between centers
  float dist = Norm(VectorSub(b2->center, b1->center));
  //printf("%f\n", dist);

  //project b2 to sat[]
  


  vec3 v1, v2;
  int i;
  for(i=0;i < 8;i++)
  {
    v1 = VectorSub(b2->vertices[i], p[0]);
    float res1 = DotProduct(v1, n[0]);

    v2 = VectorSub(b2->vertices[i], p[3]);
    float res2 = DotProduct(v2, n[3]);

    //if(sign(res1) == sign(res2) && b2->vertices[i].z > b1->pos.z && b2->vertices[i].z < b1->pos.z + b1->size.z)
    if(sign(res1) == sign(res2))
    {
      v1 = VectorSub(b2->vertices[i], p[1]);
      res1 = DotProduct(v1, n[1]);

      v2 = VectorSub(b2->vertices[i], p[4]);
      res2 = DotProduct(v2, n[4]);


      //if(sign(res1) == sign(res2) && b2->vertices[i].x > b1->pos.x && b2->vertices[i].x < b1->pos.x + b1->size.x)
      if(sign(res1) == sign(res2))
      {

      v1 = VectorSub(b2->vertices[i], p[2]);
      res1 = DotProduct(v1, n[2]);

      v2 = VectorSub(b2->vertices[i], p[5]);
      res2 = DotProduct(v2, n[5]);

        if(sign(res1) == sign(res2))
          return 1;
      }

    }

  }

  return 0;
}

void create_wall(wall_s * w, vec3 pos, vec3 size)
{
  w->pos = pos;
  w->size = size;
  w->acc = SetVector(0,0,0);
  w->force = SetVector(0,0,0);
  w->momentum = SetVector(0,0,0);
  w->angular_momentum = SetVector(0,0,0);
  w->torque = SetVector(0,0,0);
  w->R = T(0,0,0);
  w->mass = 2;
  w->omega = SetVector(0,0,0);

  w->wall_model = LoadModelPlus("./res/cube_2.obj");
  //w->orig_matrix = S(w->size.x, w->size.y, w->size.z);
  w->orig_matrix = T(0,0,0);
  //w->matrix = Mult(S(size.x, size.y, size.z), T(pos.x, pos.y, pos.z));
  w->matrix = T(pos.x, pos.y, pos.z);
  create_bb(&w->bb, w->pos, w->size);
}

void draw_wall(wall_s * w, GLuint program)
{
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, w->matrix.m);
  //glBindTexture(GL_TEXTURE_2D, f->tex);
  DrawModel(w->wall_model, program, "inPosition", "inNormal", "inTexCoord");
}

void draw_ragdoll(ragdoll_s * r, GLuint program)
{
  int i;
  for(i=0;i<12;i++)
  {
    glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, r->joints[i].T.m);
    DrawModel(r->joints[i].body, program, "inPosition", "inNormal", "inTexCoord");
  }
}



void draw_debug_sphere(ball_s * b, vec3 pos, GLuint program)
{
  mat4 mat = Mult(T(pos.x, pos.y, pos.z), S(.4,.4,.4));
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, mat.m);
  //glBindTexture(GL_TEXTURE_2D, f->tex);
  DrawModel(b->model, program, "inPosition", "inNormal", "inTexCoord");
}

void update_vertices(bounding_box_s * bb, vec3 pos, mat4 RR)
{
  mat4 R = Mult(T(pos.x, pos.y, pos.z) ,Mult(RR, T(-pos.x, -pos.y, -pos.z)));
  //000
  bb->vertices[0] = MultVec3(R, bb->pos);
  //001
  bb->vertices[1] = MultVec3(R, SetVector(bb->pos.x, bb->pos.y, bb->pos.z+bb->size.z));
  //010
  bb->vertices[2] = MultVec3(R, SetVector(bb->pos.x, bb->pos.y+bb->size.y, bb->pos.z));
  //011
  bb->vertices[3] = MultVec3(R, SetVector(bb->pos.x, bb->pos.y+bb->size.y, bb->pos.z+bb->size.z));
  //100
  bb->vertices[4] = MultVec3(R, SetVector(bb->pos.x+bb->size.x, bb->pos.y, bb->pos.z));
  //101
  bb->vertices[5] = MultVec3(R, SetVector(bb->pos.x+bb->size.x, bb->pos.y, bb->pos.z+bb->size.z));
  //110
  bb->vertices[6] = MultVec3(R, SetVector(bb->pos.x+bb->size.x, bb->pos.y+bb->size.y, bb->pos.z));
  //111
  bb->vertices[7] = MultVec3(R, SetVector(bb->pos.x+bb->size.x, bb->pos.y+bb->size.y, bb->pos.z+bb->size.z));

  bb->center = MultVec3(R, VectorAdd(bb->pos, ScalarMult(bb->size, .5)));


  //printf("%f %f %f \n", bb->vertices[7].x, bb->vertices[7].y, bb->vertices[7].z);
}

void create_bb(bounding_box_s * bb, vec3 pos, vec3 size)
{
  bb->pos = pos;
  bb->size = size;
  bb->center = VectorAdd(bb->pos, ScalarMult(bb->size, .5));

  //update_vertices(bb, pos, SetVector(0,0,0));
  update_vertices(bb, pos, IdentityMatrix());
}

void update_bb(bounding_box_s * bb, vec3 pos, vec3 angle)
{
  //mat4 R = Mult(ArbRotate(angle, Norm(angle)), InvertMat4(T(pos.x, pos.y, pos.z)));
  //mat4 R = ArbRotate(angle, Norm(angle));
  //bb->pos = MultVec3(R,pos);
  bb->pos = pos;
  //update_vertices(bb, angle);
}

void create_cow(cow_s * c)
{
  c->main_body = LoadModelPlus("./res/ko_fine.obj");
  c->matrix = S(.1,.1,.1); //makes it roughly 5 high
  c->pos = SetVector(5,0,5);
  c->speed = SetVector(0,0,0);
  c->acc = SetVector(0,0,0);
  c->force = SetVector(0,0,0);
  c->torque = SetVector(0,0,0);
  c->momentum = SetVector(0,0,0);
  c->mass = 2.0;
  c->angle = 0;
  c->d_angle = 0;
  c->jumping = 0;
  LoadTGATextureSimple("./res/texture.tga", &(c->tex));

  create_bb(&c->bb, SetVector(c->pos.x-2, c->pos.y, c->pos.z-1), SetVector(2,4,2));


}

void update_cow(cow_s * c, GLfloat dT)
{
  vec3 dP, dX;

  dP = ScalarMult(c->force, dT);

  c->momentum = VectorAdd(c->momentum, dP);
  c->speed = ScalarMult(c->momentum, 1.0/c->mass);

  dX = ScalarMult(c->speed, dT);

  c->pos = VectorAdd(c->pos, dX);

  if(c->pos.y <= 0)
  {
    c->pos.y=0;
    c->jumping = 0;
    c->momentum.y = 0;
    c->force.y=0;
  }

  //printf("%f\n", c->pos.y);

  //c->matrix = Mult(S(.1, .1, .1), T(c->pos.x, c->pos.y, c->pos.z));
  //update_bb(&c->bb, SetVector(c->pos.x-2, c->pos.y, c->pos.z-4));
  update_bb(&c->bb, SetVector(c->pos.x-2, c->pos.y, c->pos.z-1), SetVector(0,c->angle,0));
  //update_vertices(&c->bb, c->pos, SetVector(0,c->angle, 0));
  c->R = ArbRotate(SetVector(0, -c->angle, 0), -c->angle);

  update_vertices(&c->bb, c->pos, c->R);
  //update_bb(&c->bb, c->pos, SetVector(0,-c->angle,0));
}

void move_cow(cow_s * c, float angle)
{
  vec3 move_force = SetVector(0,0,0);
  static int jump_timer = 0;

  if(keyIsDown(0x20))
  {

    if(jump_timer < 20)
    {
    c->momentum.y = 20;
    c->jumping = 1;
    }

    jump_timer++;
  }else
    jump_timer = 0;

  if(!c->jumping)
  {
  if(keyIsDown('p'))
    move_force = VectorAdd(move_force, SetVector(-COW_FORCE*cos(angle),0,-COW_FORCE*sin(angle)));
  if(keyIsDown('u'))
    move_force = VectorAdd(move_force, SetVector(COW_FORCE*cos(angle),0,COW_FORCE*sin(angle)));
  if(keyIsDown('e'))
    move_force = VectorAdd(move_force, SetVector(-COW_FORCE*sin(angle),0,COW_FORCE*cos(angle)));
  if(keyIsDown('i'))
    move_force = VectorAdd(move_force, SetVector(COW_FORCE*sin(angle),0,-COW_FORCE*cos(angle)));
  }

  if(Norm(SetVector(c->momentum.x, 0, c->momentum.z)) > COW_MAX_MOMENTUM )
    move_force = SetVector(0,0,0);


  if(!c->jumping)
  {
    vec3 moment = SetVector(-c->momentum.x, 0, -c->momentum.z);
    if(Norm(moment) == 0)
      move_force = VectorAdd(move_force, moment);
    else
      move_force = VectorAdd(move_force, ScalarMult(Normalize(SetVector(-c->momentum.x,0,-c->momentum.z)), FLOOR_FRICTION));
  }


  if(c->pos.y > 0 && c->jumping)
    move_force = VectorAdd(move_force, SetVector(0,COW_GRAVITY,0));

  c->force = move_force;

  if(Norm(c->momentum) < 1 && !keyIsDown('p') && !keyIsDown('u')
     && !keyIsDown('e') && !keyIsDown('i') && !c->jumping)
  {
    c->force = SetVector(0,0,0);
    c->momentum = SetVector(0,0,0);
  }

//c->matrix = Mult(T(c->pos.x, c->pos.y, c->pos.z), S(.1,.1,.1));


/*
  if(c->pos.y < 0)
  {
    c->pos.y=0;// = SetVector(c->pos.x, 0, c->pos.z);
    c->jumping = 0;
    c->force.y=0;// = SetVector(c->force.x, 0, c->force.z);
  }
*/
}

void update_floor(floor_s * f, cow_s * c)
{
  //f->matrix = T(-c->pos.x, -c->pos.y, -c->pos.z);
  f->matrix = T(0,0,0);
//                   ArbRotate(SetVector(0,1,0), c->angle));

}

int index_of_lowest_vertex(vec3 * list)
{
  int i, lowest=0, lowestnum=0;
  vec3 tmp = {0,10000,0};
  lowest = tmp.y;
  for(i=0;i<8;i++)
  {
    if(list[i].y < lowest)
    {
      lowest=list[i].y;
      lowestnum = i;
    }
  }

  return lowestnum;
}

void update_wall(wall_s * w, cow_s * c, GLfloat dT)
{
  vec3 dP, dX, dL, dO;
  mat4 Rd;
  vec3 torque_tmp = {0,0,0};
  vec3 force_tmp = {0,0,0};
  vec3 curr_r = SetVector(w->bb.pos.x,0,w->bb.pos.z);
  int no_under = 1;


  force_tmp = ScalarMult(SetVector(-w->momentum.x, 0, -w->momentum.z), .2);
  int i;

  //torque_tmp = SetVector(0,.01,0);
  for(i=0;i<8;i++)
  {

    if(w->bb.vertices[i].y > 1)
    {
      //torque from sliding over the floor
      //torque_tmp = VectorAdd(torque_tmp, ScalarMult(
      //CrossProduct(VectorSub(w->bb.vertices[i], curr_r), SetVector(w->momentum.x, 0, w->momentum.z) ), .0005));
      //torque_tmp = ScalarMult(
      //CrossProduct(VectorSub(w->bb.vertices[i], curr_r), SetVector(w->momentum.x, 0, w->momentum.z) ), .0005);
    }


    if(w->bb.vertices[i].y < 0)
    {
      no_under = 0;
      vec3 rr = VectorSub(w->bb.vertices[i], w->bb.center);
      vec3 FF = SetVector(0,.02,0);

      curr_r = w->bb.vertices[i];



      //w->momentum.y *= -1;
      //w->bb.vertices[i].y += 2;
    }
    if(no_under)
      force_tmp = VectorAdd(force_tmp, SetVector(0,-4,0));

    //add gravity to all vertices
    vec3 r = VectorSub(w->bb.vertices[i], curr_r);
    vec3 F = SetVector(0,-.1,0); //gravity

    //torque_tmp = VectorAdd(torque_tmp, CrossProduct(r, F));
  }

  if(!no_under)
  {
    w->momentum.y *= -.5;

    int ind = index_of_lowest_vertex(w->bb.vertices);
    w->pos.y += -w->bb.vertices[ind].y+.001;
    update_vertices(&w->bb, w->pos, w->R);

    //printf("%i\n", test);
    //printf("%f\n", w->bb.vertices[ind].y);

  }

  w->torque = torque_tmp;
  w->force = force_tmp;

  dP = ScalarMult(w->force, dT);
  dL = ScalarMult(w->torque, dT);
  dX = ScalarMult(w->speed, dT);

  w->momentum = VectorAdd(w->momentum, dP);
  w->speed = ScalarMult(w->momentum, 1.0/w->mass);

  w->angular_momentum = VectorAdd(w->angular_momentum, dL);
  //w->omega = MultVec3(S(1,1,1), w->angular_momentum);

  dO = ScalarMult(w->omega, dT);
  Rd = CrossMatrix(dO);
  Rd = Mult(Rd, w->R);
  w->R = Mult(T(0, 2.5, 0), Mult(MatrixAdd(w->R, Rd), T(0, -2.5, 0)));

  //if(no_under)
  w->pos = VectorAdd(w->pos, dX);
  //w->R = ArbRotate(w->omega, Norm(w->omega));

  w->matrix = Mult(Mult(w->orig_matrix, T(w->pos.x, w->pos.y, (w->pos.z))), w->R);

  update_bb(&w->bb, SetVector(w->pos.x-1, w->pos.y, w->pos.z-1), w->omega);
  update_vertices(&w->bb, w->pos, w->R);
  OrthoNormalizeMatrix(&w->R);

  if(check_collision_2(&w->bb, &c->bb))
  {
    w->momentum = SetVector(c->momentum.x, 0, c->momentum.z);
    c->momentum = SetVector(-c->momentum.x, c->momentum.y, -c->momentum.z);

    vec3 n = SetVector(1,0,0);
    vec3 r = VectorSub(c->bb.center, w->bb.center);
    w->omega = VectorSub(w->omega, MultVec3(IdentityMatrix(), CrossProduct(r,n)));

    while(check_collision_2(&w->bb, &c->bb))
    {
      w->pos = VectorSub(w->pos, ScalarMult(Normalize(r), .05));
      update_bb(&w->bb, SetVector(w->pos.x-1, w->pos.y, w->pos.z-1), w->omega);
      update_vertices(&w->bb, w->pos, w->R);

    }
  }

}

void create_floor(floor_s * f)
{
  f->model = LoadModelPlus("./res/ground.obj");
  LoadTGATextureSimple("./res/grass.tga", &(f->tex));
}

void draw_floor(floor_s * f, GLuint program)
{
  //f->matrix = S(2,2,2);
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, f->matrix.m);
  glBindTexture(GL_TEXTURE_2D, f->tex);
  DrawModel(f->model, program, "inPosition", "inNormal", "inTexCoord");

}

void turn_cow(cow_s * c, float angle)
{
  float diff = c->angle - angle;
  //c->matrix = Mult(S(.1,.1,.1), ArbRotate(SetVector(0,1,0), c->angle));
  c->angle -= diff/10;

}

void create_bone(bone_s * b, vec3 pos)
{
  b->pos = pos;
  b->body = LoadModelPlus("./res/groundsphere.obj");
  //b->body = LoadModelPlus("./res/ko_rough.obj");
  b->matrix = Mult(T(pos.x, pos.y, pos.z), S(.1, .1, .1));
  //b->matrix = S(3, 3, 3);
  //b->matrix = T(0,0,0);

}

void create_joint(joint_s * j, vec3 pos, char * Mvar, char * posvar, char * boneposvar, int id)
{
  j->isnull = 1;
  j->Mvar = Mvar;
  j->posvar = posvar;
  j->boneposvar = boneposvar;

  j->pos = pos;
  j->id = id;
  j->body = LoadModelPlus("./res/groundsphere.obj");
  j->body_matrix = Mult(T(pos.x, pos.y, pos.z), S(.1, .1, .1));
  j->T = T(pos.x, pos.y, pos.z);
  j->M = T(pos.x, pos.y, pos.z);
  j->Mp = T(pos.x, pos.y, pos.z);
  j->Mtot = T(pos.x, pos.y, pos.z);
  //j->Minv = InvertMat4(j->M);
  j->Minv = InvertMat4(T(pos.x, pos.y, pos.z));
  j->R = T(0,0,0);
  //j->body_matrix = S(.1, .1, .1);
  //j->body_matrix = IdentityMatrix();
}

void create_ragdoll_joint(joint_s * j, vec3 pos)
{
  j->angle = 0;
  j->max_angle = 360;
  j->pos = pos;
  j->constraint = 0;
  j->orig_vec = SetVector(0,0,0);
  j->body = LoadModelPlus("./res/groundsphere.obj");
  j->body_matrix = Mult(T(pos.x, pos.y, pos.z), S(.02, .02, .02));
  j->T = T(pos.x, pos.y, pos.z);
}

void create_ragdoll(ragdoll_s * r)
{
  //head
  create_ragdoll_joint(&r->joints[0], SetVector(0,10,0));

  create_ragdoll_joint(&r->joints[1], SetVector(0,8,0));
  create_ragdoll_joint(&r->joints[2], SetVector(0,8,2));
  create_ragdoll_joint(&r->joints[3], SetVector(0,8,-2));
  create_ragdoll_joint(&r->joints[4], SetVector(0,8,4));
  create_ragdoll_joint(&r->joints[5], SetVector(0,8,-4));

  r->joints[4].max_angle = 25;
  r->joints[4].orig_vec = SetVector(0,1,0);
  r->joints[4].constraint = 1;

  r->joints[2].constraint = 1;
  r->joints[2].max_angle = 20;
  r->joints[2].orig_vec = SetVector(0,0,-1);

  r->joints[3].max_angle = 25;
  r->joints[3].orig_vec = SetVector(0,0,1);
  r->joints[3].constraint = 1;




  //neck and shoulders, arms
/*
  create_ragdoll_joint(&r->joints[1], SetVector(0,8,0));
  create_ragdoll_joint(&r->joints[4], SetVector(0,8,2));
  create_ragdoll_joint(&r->joints[5], SetVector(0,8,-2));
  create_ragdoll_joint(&r->joints[6], SetVector(0,8,4));
  create_ragdoll_joint(&r->joints[7], SetVector(0,8,-4));


  //stomach
  create_ragdoll_joint(&r->joints[2], SetVector(0,6,0));
  //groin
  create_ragdoll_joint(&r->joints[3], SetVector(0,4,0));

  //legs
  create_ragdoll_joint(&r->joints[8], SetVector(0,4,7));
  create_ragdoll_joint(&r->joints[9], SetVector(0,4,-2));
  create_ragdoll_joint(&r->joints[10], SetVector(0,2,2));
  create_ragdoll_joint(&r->joints[11], SetVector(0,2,-2));
*/
/*
  r->joints[0].parent = NULL;
  r->joints[1].parent = &r->joints[0];
  r->joints[2].parent = &r->joints[1];
  r->joints[3].parent = &r->joints[2];
  r->joints[4].parent = &r->joints[1];
  r->joints[5].parent = &r->joints[1];
  r->joints[6].parent = &r->joints[4];
  r->joints[7].parent = &r->joints[5];

  r->joints[8].parent = &r->joints[3];
  r->joints[9].parent = &r->joints[3];
  r->joints[10].parent = &r->joints[8];
  r->joints[11].parent = &r->joints[9];
*/

  r->joints[0].child[0] = &r->joints[1];


  r->joints[0].parent = NULL;
  r->joints[1].parent = &r->joints[0];
  r->joints[2].parent = &r->joints[1];
  r->joints[3].parent = &r->joints[1];
  r->joints[4].parent = &r->joints[2];
  r->joints[5].parent = &r->joints[3];


  r->joints[1].dist_to_parent = 2;
  r->joints[2].dist_to_parent = 2;
  r->joints[3].dist_to_parent = 2;
  r->joints[4].dist_to_parent = 2;
  r->joints[5].dist_to_parent = 2;
/*
  r->joints[5].dist_to_parent = 2;
  r->joints[6].dist_to_parent = 2;
  r->joints[7].dist_to_parent = 2;
  r->joints[8].dist_to_parent = 2;
  r->joints[9].dist_to_parent = 2;
  r->joints[10].dist_to_parent = 2;
  r->joints[11].dist_to_parent = 2;
*/
  //set a force on the head
  //r->joints[0].force = SetVector(0,0,.1);
}

void update_ragdoll(ragdoll_s * r, GLfloat dT)
{
  if(keyIsDown('a'))
    r->joints[0].force.z = 20;
  else if(keyIsDown('q'))
    r->joints[0].force.x = 20;
  else if(keyIsDown('o'))
    r->joints[0].force.z = -20;
  else if(keyIsDown('.'))
    r->joints[0].force.y = 80;
  else
  {
    r->joints[0].force.z = 0;
    r->joints[0].force.y = 0;
    r->joints[0].force.x = 0;
  }

  int test=0;
  vec3 dP, dX, calculated_force = {0,0,0};
  vec3 n;
  vec3 real_dist;
  int i;
  for(i=0;i<6;i++)
  {

    joint_s * j = &r->joints[i];
    joint_s * parent = r->joints[i].parent;
    if(parent != NULL && j != NULL)
    {
      //parent = parent->parent;
      real_dist = VectorSub(j->pos, parent->pos);
      n = Normalize(real_dist);
      float dist_diff = (Norm(real_dist)-j->dist_to_parent);
      vec3 speed_diff = VectorSub(j->speed, parent->speed);
      calculated_force = VectorSub(ScalarMult(n, -dist_diff*40),
      ScalarMult(speed_diff, 5));
      j->pos = VectorSub(j->pos, ScalarMult(n, dist_diff));

      float ang = 180*acos(DotProduct(Normalize(real_dist), j->orig_vec))/M_PI;
      while(ang > j->max_angle && j->constraint)
      {
        printf("%f\n", ang);

        real_dist = VectorSub(parent->pos, j->pos);
        dist_diff = (Norm(real_dist)-j->dist_to_parent);

        n = Normalize(real_dist);
        vec3 align_vec = Normalize(VectorSub(real_dist, j->orig_vec));
        if(dist_diff < 0)
          j->pos = VectorAdd(j->pos, ScalarMult(n, -.01));
        else
          j->pos = VectorAdd(j->pos, ScalarMult(n, .01));
        j->pos = VectorAdd(j->pos, ScalarMult(align_vec, .01));
        ang = 180*acos(DotProduct(Normalize(real_dist), j->orig_vec))/M_PI;
      }

      joint_s * jp = j;
      while(jp != NULL)
      {
      while(jp->pos.y <= 0)
      {
        jp->pos.y += .2;
        jp->speed.y *= -.6;
      }
        jp = jp->parent;
      }

    }

      //j = j->parent;
      //parent = parent->parent;
    //}

    dP = ScalarMult(VectorAdd(SetVector(0,-40,0), VectorAdd(j->force, calculated_force)), dT);
    j->speed = VectorAdd(j->speed, dP);
    dX = ScalarMult(j->speed, dT);
    j->pos = VectorAdd(j->pos, dX);
    j->T = Mult(T(j->pos.x, j->pos.y, j->pos.z), S(.5,.5,.5));
    //parent = parent->parent;

  }
}

void draw_cow(cow_s * c, GLuint program)
{
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, c->matrix.m);
  glBindTexture(GL_TEXTURE_2D, c->tex);
  DrawModel(c->main_body, program, "inPosition", "inNormal", "inTexCoord");
}

void draw_bone(bone_s * b, GLuint program)
{
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, b->matrix.m);
  DrawModel(b->body, program, "inPosition", "inNormal", "inTexCoord");

}

void draw_joint(joint_s *j, GLuint program)
{
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, j->body_matrix.m);
  DrawModel(j->body, program, "inPosition", "inNormal", "inTexCoord");
  
}

void create_ball(ball_s * b, vec3 pos)
{
  b->model = LoadModelPlus("./res/groundsphere.obj");
  b->pos = SetVector(pos.x,pos.y,pos.z);
  b->matrix = T(pos.x, pos.y, pos.z);
  b->speed = SetVector(0,0,0);
  b->acc = SetVector(0,0,0);
  b->force = SetVector(0,0,0);
  b->torque = SetVector(0,0,0);
  b->momentum = SetVector(2,0,0);
  b->mass = 2;
}


void update_ball(ball_s * b, cow_s * c, GLfloat dT)
{
  vec3 dP, dX;

  //b->momentum = SetVector(1,0,0);

  float dist = Norm(VectorSub(b->pos,c->pos));

  //if collision between cow and ball
  if(dist < 3)
  {
    vec3 p = {0, 0, 0};
    //point of impact
    p.x = (c->pos.x + b->pos.x)/2;
    p.z = (c->pos.z + b->pos.z)/2;
    vec3 nA = Normalize(VectorSub(c->pos, p));
    float vrel = DotProduct(VectorSub(c->speed, b->speed), nA);
    float eps = .1;
    float jj = vrel*(-(eps+1))/(1/c->mass + 1/b->mass);
    c->momentum = VectorAdd(c->momentum, ScalarMult(nA,jj));
    b->momentum = VectorAdd(b->momentum, ScalarMult(nA,-jj));

    while(dist < 3)
    {
      dist = Norm(VectorSub(b->pos,c->pos));
      c->pos = VectorAdd(c->pos, ScalarMult(Normalize(nA), .01));
      b->pos = VectorSub(b->pos, ScalarMult(Normalize(nA), .01));
    }
  }


    vec3 moment = SetVector(-b->momentum.x, 0, -b->momentum.z);
    if(Norm(moment) != 0)
      b->force = ScalarMult(Normalize(SetVector(-b->momentum.x,0,-b->momentum.z)), FLOOR_FRICTION/10);

  //b->matrix = T(-c->pos.x + b->pos.x, -c->pos.y + b->pos.y, -c->pos.z+b->pos.z);
  b->matrix = T(b->pos.x, b->pos.y, b->pos.z);

  dP = ScalarMult(b->force, dT);

  b->momentum = VectorAdd(b->momentum, dP);
  b->speed = ScalarMult(b->momentum, 1.0/b->mass);
  //printf("%f\n", b->mass);

  dX = ScalarMult(b->speed, dT);

  b->pos = VectorAdd(b->pos, dX);
}

void draw_ball(ball_s * b, GLuint program)
{
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, b->matrix.m);
  DrawModel(b->model, program, "inPosition", "inNormal", "inTexCoord");
}
