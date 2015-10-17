#include "objects.h"
#include <stdlib.h>

void create_cow(cow_s * c)
{
  c->main_body = LoadModelPlus("./res/ko_fine.obj");
  c->matrix = S(.1,.1,.1); //makes it roughly 5 high
  c->pos = SetVector(0,0,0);
  c->speed = SetVector(0,0,0);
  c->acc = SetVector(0,0,0);
  c->force = SetVector(0,0,0);
  c->torque = SetVector(0,0,0);
  c->momentum = SetVector(0,0,0);
  c->mass = 2.0;
  c->angle = 0;
  c->d_angle = 0;
  c->jumping = 0;
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
    c->pos.y=0;// = SetVector(c->pos.x, 0, c->pos.z);
    c->jumping = 0;
    c->momentum.y = 0;
    c->force.y=0;// = SetVector(c->force.x, 0, c->force.z);
  }

  //printf("%f\n", c->pos.y);

  //c->matrix = Mult(S(.1, .1, .1), T(c->pos.x, c->pos.y, c->pos.z));
}

void move_cow(cow_s * c, float angle)
{
  vec3 move_force = SetVector(0,0,0);

  if(keyIsDown(0x20))
  {
    c->momentum.y = 20;
    c->jumping = 1;
  }

  //if(Norm(c->speed) < 15)
  //{

  if(keyIsDown('p'))
    move_force = VectorAdd(move_force, SetVector(-COW_FORCE*cos(angle),0,-COW_FORCE*sin(angle)));
  if(keyIsDown('u'))
    move_force = VectorAdd(move_force, SetVector(COW_FORCE*cos(angle),0,COW_FORCE*sin(angle)));
  if(keyIsDown('e'))
    move_force = VectorAdd(move_force, SetVector(-COW_FORCE*sin(angle),0,COW_FORCE*cos(angle)));
  if(keyIsDown('i'))
    move_force = VectorAdd(move_force, SetVector(COW_FORCE*sin(angle),0,-COW_FORCE*cos(angle)));
  //}


  //if(!keyIsDown('p') && !keyIsDown('u') && !keyIsDown('e') && !keyIsDown('i') && !c->jumping)
  //{
    //c->force = SetVector(0,0,0);
    vec3 moment = SetVector(-c->momentum.x, 0, -c->momentum.z);
    if(Norm(moment) == 0)
      move_force = VectorAdd(move_force, moment);
    else
      move_force = VectorAdd(move_force, ScalarMult(Normalize(SetVector(-c->momentum.x,0,-c->momentum.z)), FLOOR_FRICTION));
  //}



  if(c->pos.y > 0 && c->jumping)
    move_force = VectorAdd(move_force, SetVector(0,-40,0));


  c->force = move_force;

  //if(Norm(SetVector(c->momentum.x, 0, c->momentum.z)) > 30)
  //  c->momentum = ScalarMult(Normalize(c->momentum), 15);
/*
  if(c->pos.y < 0)
  {
    c->pos.y=0;// = SetVector(c->pos.x, 0, c->pos.z);
    c->jumping = 0;
    c->force.y=0;// = SetVector(c->force.x, 0, c->force.z);
  }
*/
  //printf("%f t\n", c->pos.y);

  //printf("%f\n", c->speed.x);
}

void update_floor(floor_s * f, cow_s * c)
{
  f->matrix = T(-c->pos.x, -c->pos.y, -c->pos.z);
//                   ArbRotate(SetVector(0,1,0), c->angle));

}

void create_floor(floor_s * f)
{
  f->model = LoadModelPlus("./res/ground.obj");
  LoadTGATextureSimple("./res/golv.tga", &(f->tex));
}

void draw_floor(floor_s * f, GLuint program)
{
  //f->matrix = S(2,2,2);
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, f->matrix.m);
  //glBindTexture(GL_TEXTURE_2D, f->tex);
  DrawModel(f->model, program, "inPosition", "inNormal", "inTexCoord");

}

void turn_cow(cow_s * c, float angle)
{
  float diff = c->angle - angle;
  c->matrix = Mult(S(.1,.1,.1), ArbRotate(SetVector(0,1,0), c->angle));
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

void draw_cow(cow_s * c, GLuint program)
{
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, c->matrix.m);
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
