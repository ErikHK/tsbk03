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

void draw_debug_sphere(ball_s * b, vec3 pos, GLuint program)
{
  mat4 mat = Mult(T(pos.x, pos.y, pos.z), S(.4,.4,.4));
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, mat.m);
  //glBindTexture(GL_TEXTURE_2D, f->tex);
  DrawModel(b->model, program, "inPosition", "inNormal", "inTexCoord");
}

void update_vertices(bounding_box_s * bb, vec3 pos, vec3 angle)
{
  mat4 R = Mult(T(pos.x, pos.y, pos.z) ,Mult(ArbRotate(angle, Norm(angle)), T(-pos.x, -pos.y, -pos.z)));
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

  update_vertices(bb, pos, SetVector(0,0,0));
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
  update_vertices(&c->bb, c->pos, SetVector(0,c->angle, 0));
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

void update_wall(wall_s * w, cow_s * c, GLfloat dT)
{
  vec3 dP, dX, dL;
  vec3 torque_tmp = {0,0,0};
  vec3 force_tmp = {0,0,0};
  vec3 curr_r = SetVector(w->bb.pos.x+1,0,w->bb.pos.z+1);
  int no_under = 1;


  //set friction force
  //if(Norm(w->momentum) > 1)
  //{
  force_tmp = ScalarMult(SetVector(-w->momentum.x, 0, -w->momentum.z), 2);
  //}
  //else
  //{
    //force_tmp = SetVector(0,0,0);
  //}


  //printf("%f %f %f\n", curr_r.x, curr_r.y, curr_r.z);

  int i;
  for(i=0;i<8;i++)
  {

    //torque from sliding over the floor
    torque_tmp = ScalarMult(
    CrossProduct(VectorSub(w->bb.vertices[i], curr_r), SetVector(w->momentum.x, 0, w->momentum.z) ), .01);

    if(w->bb.vertices[i].y <= 0)
    {
      no_under = 0;
      vec3 rr = VectorSub(w->bb.vertices[i], w->bb.center);
      vec3 FF = SetVector(0,.02,0);
      //torque_tmp = VectorAdd(torque_tmp, CrossProduct(rr,FF));
      //force_tmp = VectorAdd(force_tmp, SetVector(0,-w->bb.vertices[i].y,0));
      //w->momentum.y = 0;
      //w->omega = VectorAdd(w->omega, CrossProduct(rr,FF));
      //w->omega = CrossProduct(rr,FF);

      curr_r = w->bb.vertices[i];

      //force_tmp.y += 4;
      w->pos.y += .03;
    }
    
    if(no_under)
      force_tmp = VectorAdd(force_tmp, SetVector(0,-4,0));
    else
      force_tmp.y=4;// = SetVector(0,4,0);

    //add gravity to all vertices
    //float cdist = Norm(VectorSub(
    //   SetVector(w->bb.vertices[i].x, 0, w->bb.vertices[i].z),
    //   SetVector(w->bb.center.x, 0, w->bb.center.z) ));
    vec3 r = VectorSub(w->bb.vertices[i], curr_r);
    vec3 F = SetVector(0,-.1,0); //gravity
    //if(w->bb.vertices[i].y > w->bb.center.y)
      torque_tmp = VectorAdd(torque_tmp, CrossProduct(r, F));
  }

  dP = ScalarMult(w->force, dT);
  dL = ScalarMult(w->torque, dT);

  w->torque = torque_tmp;
  w->force = force_tmp;

  w->momentum = VectorAdd(w->momentum, dP);
  w->speed = ScalarMult(w->momentum, 1.0/w->mass);

  w->angular_momentum = VectorAdd(w->angular_momentum, dL);

  dX = ScalarMult(w->speed, dT);

  w->pos = VectorAdd(w->pos, dX);
  w->R = ArbRotate(w->omega, Norm(w->omega));

  //w->matrix = Mult(Mult(w->orig_matrix, T(-c->pos.x+w->pos.x, -c->pos.y+w->pos.y, (-c->pos.z+w->pos.z))), w->R);
  w->matrix = Mult(Mult(w->orig_matrix, T(w->pos.x, w->pos.y, (w->pos.z))), w->R);

  update_bb(&w->bb, SetVector(w->pos.x-1, w->pos.y, w->pos.z-1), w->omega);
  update_vertices(&w->bb, w->pos, w->omega);

  w->omega = MultVec3(S(1,1,1), w->angular_momentum);

  if(check_collision_2(&w->bb, &c->bb))
  {
    vec3 r1 = SetVector(0, c->pos.y, 0);
    vec3 r2 = SetVector(0, w->pos.y+.5, 0);
    vec3 r = VectorSub(r1,r2);
    //w->angular_momentum = CrossProduct(r, w->momentum);
    //w->torque = ScalarMult(CrossProduct(r, c->momentum), .05);
    //w->torque = ScalarMult(c->momentum, 1000);
    //w->force = ScalarMult(c->momentum, .5);

    w->momentum = SetVector(c->momentum.x, 0, c->momentum.z);
    c->momentum = SetVector(-c->momentum.x, c->momentum.y, -c->momentum.z);
    //w->torque
  }

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
