#include "objects.h"
#include <stdlib.h>
#include "signal.h"

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

int check_sphere_collision(vec3 pos1, vec3 pos2, float r1, float r2)
{
  vec3 diff = VectorSub(pos1, pos2);
  if(Norm(diff) < r1+r2)
    return 1;
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

void create_farmer(farmer_s * f, vec3 pos)
{
  f->animate = 1;
  f->body = LoadModelPlus("./res/farmergirl1_small.obj");
  f->orig_body = LoadModelPlus("./res/farmergirl1_small.obj");
  LoadTGATextureSimple("./res/farmergirl.tga", &(f->tex));
  f->pos = pos;
  f->body_R = IdentityMatrix();

  f->verts = malloc(sizeof(GLfloat)*f->body->numVertices*3);
  memcpy(f->verts, f->body->vertexArray, sizeof(GLfloat)*f->body->numVertices*3);

  f->skeleton.num_joints = 16;

  //head joint
  create_joint(&f->skeleton.joints[0],
	VectorAdd(f->pos, SetVector(0,6,0)),
	"farmer_head", "farmer_head_pos", "farmer_head_bone_pos", 0);

  //f->skeleton.joints[0].pos = SetVector(-1.4,6,.1);

  //neck joint
  create_joint(&f->skeleton.joints[1],
  	VectorAdd(f->pos, SetVector(-.1, 3.7, 0)),
	"farmer_neck", "farmer_neck_pos", "farmer_neck_bone_pos", 0);

  f->skeleton.joints[1].constraint = 1;
  f->skeleton.joints[1].max_angle = 40;
  f->skeleton.joints[1].orig_vec = SetVector(0,1,0);

  f->skeleton.joints[1].weight = 2;


  //shoulder joints
  create_joint(&f->skeleton.joints[2],
	VectorAdd(f->pos, SetVector(0, 3.7, .6)),
	"farmer_lshoulder", "farmer_lshoulder_pos", "farmer_lshoulder_bone_pos", 0);
  create_joint(&f->skeleton.joints[3],
	VectorAdd(f->pos, SetVector(0, 3.7, -.6)),
	"farmer_rshoulder", "farmer_rshoulder_pos", "farmer_rshoulder_bone_pos", 0);

  //elbow joints
  create_joint(&f->skeleton.joints[4],
	VectorAdd(f->pos, SetVector(0, 3.8, 1.7)),
	NULL, NULL, NULL, 0);
  create_joint(&f->skeleton.joints[5],
	VectorAdd(f->pos, SetVector(0, 3.4, -1.7)),
	NULL, NULL, NULL, 0);

  f->skeleton.joints[4].orig_pos = SetVector(0,3.7,1.7);
  f->skeleton.joints[5].orig_pos = SetVector(0,3.7,-1.7);


  f->skeleton.joints[2].constraint = 1;
  f->skeleton.joints[2].max_angle = 10;
  f->skeleton.joints[2].orig_vec = SetVector(0,0,-1);

  f->skeleton.joints[3].constraint = 1;
  f->skeleton.joints[3].max_angle = 10;
  f->skeleton.joints[3].orig_vec = SetVector(0,0,1);


  f->skeleton.joints[4].constraint = 1;
  f->skeleton.joints[4].max_angle = 80;
  f->skeleton.joints[4].orig_vec = SetVector(0,0,-1);

  f->skeleton.joints[5].constraint = 1;
  f->skeleton.joints[5].max_angle = 80;
  f->skeleton.joints[5].orig_vec = SetVector(0,0,1);

  f->skeleton.joints[4].weight = 1.5;
  f->skeleton.joints[5].weight = 1.5;


  //hand joints
  create_joint(&f->skeleton.joints[6],
	VectorAdd(f->pos, SetVector(0, 3.7, 2.8)),
	NULL, NULL, NULL, 0);
  create_joint(&f->skeleton.joints[7],
	VectorAdd(f->pos, SetVector(0, 3.7, -2.8)),
	NULL, NULL, NULL, 0);

  f->skeleton.joints[6].weight = 0;
  f->skeleton.joints[7].weight = 0;

  f->skeleton.joints[6].orig_pos = SetVector(0,3.7,2.8);
  f->skeleton.joints[7].orig_pos = SetVector(0,3.7,-2.8);


  f->skeleton.joints[6].constraint = 1;
  f->skeleton.joints[6].max_angle = 90;
  f->skeleton.joints[6].orig_vec = SetVector(0,0,-1);

  f->skeleton.joints[7].constraint = 1;
  f->skeleton.joints[7].max_angle = 90;
  f->skeleton.joints[7].orig_vec = SetVector(0,0,1);



  //f->skeleton.joints[7].T = T(0, 3.7, -2.6);
  //f->skeleton.joints[7].pos = SetVector(0, 3.7, -2.6);

  //f->skeleton.joints[5].T = T(0, 3.2, -1.5);
  //f->skeleton.joints[5].pos = SetVector(0, 3.2, -1.5);


  //stomach joint
  create_joint(&f->skeleton.joints[8],
	VectorAdd(f->pos, SetVector(0, 2.7, .05)),
	"farmer_stomach", "farmer_stomach_pos", "farmer_stomach_bone_pos", 0);

  f->skeleton.joints[8].orig_pos = SetVector(0,2.7,0);

  f->skeleton.joints[8].constraint = 1;
  f->skeleton.joints[8].max_angle = 10;
  f->skeleton.joints[8].orig_vec = SetVector(0,1,0);


  //groin joint
  create_joint(&f->skeleton.joints[9],
	VectorAdd(f->pos, SetVector(0, 2, 0)),
	"farmer_groin", "farmer_groin_pos", "farmer_groin_bone_pos", 0);

  f->skeleton.joints[9].constraint = 1;
  f->skeleton.joints[9].max_angle = 10;
  f->skeleton.joints[9].orig_vec = SetVector(0,1,0);

  f->skeleton.joints[9].weight = 0;


  //hip joints
  create_joint(&f->skeleton.joints[10],
	VectorAdd(f->pos, SetVector(0, 2, .4)),
	"farmer_lhip", "farmer_lhip_pos", "farmer_lhip_bone_pos", 0);
  create_joint(&f->skeleton.joints[11],
	VectorAdd(f->pos, SetVector(0, 2, -.4)),
	"farmer_rhip", "farmer_rhip_pos", "farmer_rhip_bone_pos", 0);

  f->skeleton.joints[10].weight = .8;
  f->skeleton.joints[11].weight = .8;


  f->skeleton.joints[10].constraint = 1;
  f->skeleton.joints[10].max_angle = 5;
  f->skeleton.joints[10].orig_vec = SetVector(0,0,-1);


  f->skeleton.joints[11].constraint = 1;
  f->skeleton.joints[11].max_angle = 5;
  f->skeleton.joints[11].orig_vec = SetVector(0,0,1);

  //knee joints
  create_joint(&f->skeleton.joints[12],
	VectorAdd(f->pos, SetVector(0, .8, .35)),
	NULL, NULL, NULL, 0);
//	"farmer_lhip", "farmer_lhip_pos", "farmer_lhip_bone_pos", 0);


  create_joint(&f->skeleton.joints[13],
	VectorAdd(f->pos, SetVector(0, .8, -.35)),
	NULL, NULL, NULL, 0);
//	"farmer_rhip", "farmer_rhip_pos", "farmer_rhip_bone_pos", 0);


  f->skeleton.joints[12].constraint = 1;
  f->skeleton.joints[12].max_angle = 40;
  f->skeleton.joints[12].orig_vec = SetVector(0,1,0);

  f->skeleton.joints[13].constraint = 1;
  f->skeleton.joints[13].max_angle = 40;
  f->skeleton.joints[13].orig_vec = SetVector(0,1,0);

  f->skeleton.joints[12].weight = 1.6;
  f->skeleton.joints[13].weight = 1.6;

  //foot joints
  create_joint(&f->skeleton.joints[14],
	VectorAdd(f->pos, SetVector(0, 0, .3)),
	NULL, NULL, NULL, 0);
  create_joint(&f->skeleton.joints[15],
	VectorAdd(f->pos, SetVector(0, 0, -.3)),
	NULL, NULL, NULL, 0);

  f->skeleton.joints[14].weight = 0;
  f->skeleton.joints[15].weight = 0;


  f->skeleton.joints[14].constraint = 1;
  f->skeleton.joints[14].max_angle = 25;
  f->skeleton.joints[14].orig_vec = SetVector(0,1,0);

  f->skeleton.joints[15].constraint = 1;
  f->skeleton.joints[15].max_angle = 25;
  f->skeleton.joints[15].orig_vec = SetVector(0,1,0);


  mat4 rrr = Mult(T(0,6,0), Mult(Rz(M_PI/2), T(0,-6,0)));
  int i;
  for(i=0;i<16;i++)
  {
    //f->skeleton.joints[i].pos.z += 0.0;
    //f->skeleton.joints[i].orig_pos = MultVec3(rrr, f->skeleton.joints[i].orig_pos);
    //f->skeleton.joints[i].pos = MultVec3(rrr, f->skeleton.joints[i].pos);

  }




  //set children

  //head
  f->skeleton.joints[0].child[0] = &f->skeleton.joints[1];
  //neck
  //f->skeleton.joints[1].child[0] = &f->skeleton.joints[2];
  //f->skeleton.joints[1].child[1] = &f->skeleton.joints[3];
  //f->skeleton.joints[1].child[2] = &f->skeleton.joints[8];
  f->skeleton.joints[1].child[0] = &f->skeleton.joints[8];
  //shoulders
  f->skeleton.joints[2].child[0] = &f->skeleton.joints[4];
  f->skeleton.joints[3].child[0] = &f->skeleton.joints[5];
  //elbows
  f->skeleton.joints[4].child[0] = &f->skeleton.joints[6];
  f->skeleton.joints[5].child[0] = &f->skeleton.joints[7];
  //stomach
  f->skeleton.joints[8].child[0] = &f->skeleton.joints[9];
  //groin
  //f->skeleton.joints[9].child[0] = &f->skeleton.joints[10];
  f->skeleton.joints[10].child[1] = &f->skeleton.joints[11];
  //hips
  f->skeleton.joints[10].child[0] = &f->skeleton.joints[12];
  f->skeleton.joints[11].child[0] = &f->skeleton.joints[13];
  //knees
  f->skeleton.joints[12].child[0] = &f->skeleton.joints[14];
  f->skeleton.joints[13].child[0] = &f->skeleton.joints[15];

  f->skeleton.joints[6].child[0] = NULL;
  f->skeleton.joints[6].child[1] = NULL;


  //set parents
  f->skeleton.joints[0].parent = NULL;
  f->skeleton.joints[1].parent = &f->skeleton.joints[0];
  //shoulders
  f->skeleton.joints[2].parent = &f->skeleton.joints[1];
  f->skeleton.joints[3].parent = &f->skeleton.joints[1];
  //elbows
  f->skeleton.joints[4].parent = &f->skeleton.joints[2];
  f->skeleton.joints[5].parent = &f->skeleton.joints[3];
  //hands
  f->skeleton.joints[6].parent = &f->skeleton.joints[4];
  f->skeleton.joints[7].parent = &f->skeleton.joints[5];
  //stomach
  f->skeleton.joints[8].parent = &f->skeleton.joints[1];
  //groin
  f->skeleton.joints[9].parent = &f->skeleton.joints[8];
  //hips
  f->skeleton.joints[10].parent = &f->skeleton.joints[9];
  f->skeleton.joints[11].parent = &f->skeleton.joints[9];
  //knees
  f->skeleton.joints[12].parent = &f->skeleton.joints[10];
  f->skeleton.joints[13].parent = &f->skeleton.joints[11];
  //feet
  f->skeleton.joints[14].parent = &f->skeleton.joints[12];
  f->skeleton.joints[15].parent = &f->skeleton.joints[13];

  //int kk;
  //for(kk=0;kk < 15;kk++)
  //  f->skeleton.joints[kk].dist_to_parent = 2;
  f->skeleton.joints[0].dist_to_parent = 2.3;
  f->skeleton.joints[1].dist_to_parent = 2.3;
  f->skeleton.joints[2].dist_to_parent = .6;
  f->skeleton.joints[3].dist_to_parent = .6;
  f->skeleton.joints[4].dist_to_parent = 1.1;
  f->skeleton.joints[5].dist_to_parent = 1.1;
  f->skeleton.joints[6].dist_to_parent = 1.1;
  f->skeleton.joints[7].dist_to_parent = 1.1;
  f->skeleton.joints[8].dist_to_parent = 1;
  f->skeleton.joints[9].dist_to_parent = .7;
  f->skeleton.joints[10].dist_to_parent = .4;
  f->skeleton.joints[11].dist_to_parent = .4;
  f->skeleton.joints[12].dist_to_parent = 1.2;
  f->skeleton.joints[13].dist_to_parent = 1.2;
  f->skeleton.joints[14].dist_to_parent = .8;
  f->skeleton.joints[15].dist_to_parent = .8;
}


void calc_bone_transform(farmer_s * f, joint_s * j, int acc, int start_deg)
{
  joint_s * jc;
  joint_s * rootj = j;
  mat4 tmp, tmptrans, invtrans;

  if(acc)
    tmp = j->parent->tmp;
  else
    tmp = IdentityMatrix();
    //tmp = T(farmer.pos.x,
    //farmer.pos.y-6,
    //farmer.pos.z);

  GLfloat Ms[8][16];
  int i=0,ii=0, k=0;
  float currpos[8*3] = {0};
  float bonepos[8*3] = {0};

  while(j->child[0] != NULL)
  {
    vec3 tmp_bonepos;
    vec3 tmp_bonepos_orig;
    jc = j->child[0];

    tmp = Mult(tmp, j->T);

    invtrans = InvertMat4(j->T);

    f->body_R = IdentityMatrix();

    //OrthoNormalizeMatrix(&RR);
    //tmp = Mult(tmp, Mult(RR, invtrans));

    tmp = Mult(tmp, Mult(j->R, invtrans));

    j->tmp = tmp;
    j->isnull = 0;

    //j->R = RR;
    j = j->child[0];
    i++;
  }

}

void calc_bone_transform_ragdoll(farmer_s * f, joint_s * j, int acc, int start_deg)
{
  joint_s * jc;
  joint_s * rootj = j;
  mat4 tmp, tmptrans, invtrans;

  if(acc)
    tmp = j->parent->tmp;
  else
    tmp = IdentityMatrix();
    //tmp = T(farmer.pos.x,
    //farmer.pos.y-6,
    //farmer.pos.z);

  GLfloat Ms[8][16];
  int i=0, ii=0, k=0;
  float currpos[8*3] = {0};
  float bonepos[8*3] = {0};

  while(j->child[0] != NULL)
  {
    jc = j->child[0];
/*
    for(k=1; k < MAX_CHILDREN; k++)
    {
      if(j->child[k] != NULL)
        calc_bone_transform_ragdoll(f, j->child[k], 1, 0);
    }
*/

    invtrans = InvertMat4(j->orig_T);

    //tmp = Mult(invtrans, tmp);
    tmp = Mult(tmp, j->orig_T);

    vec3 a = Normalize(VectorSub(jc->orig_pos, j->orig_pos));
    vec3 b;
    b = Normalize(VectorSub(jc->pos, j->pos));

    vec3 v = Normalize(CrossProduct(a,b));
    float s = Norm(v);
    float c = DotProduct(a, b);

    float deg;

    deg = acos(c);

    mat4 RR;
      RR = ArbRotate(Normalize(v), deg);

    //tmp = Mult(tmp, Mult(RR, invtrans));
    //tmp = Mult(j->orig_T, Mult(RR,tmp));
    //tmp = Mult(tmp, Mult(RR, invtrans));
    //tmp = Mult(j->orig_T, Mult(RR, invtrans));
    tmp = Mult(tmp, Mult(RR, invtrans));

    //if(i!=0)
    //tmp = Mult(j->parent->T, tmp);


    j->tmp = tmp;
    j->isnull = 0;
    j->R = RR;

    //if(i==0)
    //printf("%i %f\n", i, 180*deg/M_PI);

    j = j->child[0];
    i++;
  }

}


void update_farmer(farmer_s * f, GLfloat t)
{

  if(f->animate)
  {
  int i;
  for(i=0;i<16;i++)
  {
    joint_s * j = &f->skeleton.joints[i];
    j->T = T(j->pos.x, j->pos.y, j->pos.z);
    //j->pos.x += .02;
  }

    f->pos.x += .02;

    mat4 test = Mult(T(0,6,0), Mult(f->body_R, T(0,-6,0)));
    f->matrix = Mult(T(f->pos.x, f->pos.y, f->pos.z), test);

    joint_s * j;
    joint_s * jc;
    j = &f->skeleton.joints[2];
    //j->R = Rx(cos(4*t));
    j->R = Mult(Rx(M_PI/2.2 + sin(5*t)/11), Ry(cos(-5*t)/2));

    j = &f->skeleton.joints[3];
    //j->R = Rx(cos(4*t));
    j->R = Mult(Rx(-M_PI/2.2 + sin(5*t)/11), Ry(cos(5*t)/2));

    //left shoulder (her right)
    //jc = &farmer.skeleton.joints[3];
    //jc->R = Rx(-sin(3*t));

    jc = &f->skeleton.joints[4];
    jc->R = Ry(M_PI/3 + cos(5*t)/8);

    jc = &f->skeleton.joints[5];
    jc->R = Ry(-M_PI/3 + cos(5*t)/8);

    jc = &f->skeleton.joints[10];
    //jc->R = Rz(M_PI/3);
    jc->R = Rz(.5-cos(5*t));

    jc = &f->skeleton.joints[12];
    jc->R = Rz((-1-cos(5*t))/2);


    jc = &f->skeleton.joints[11];
    //jc->R = Rz(M_PI/3);
    jc->R = Rz(.5+cos(5*t));

    jc = &f->skeleton.joints[13];
    jc->R = Rz((-1-cos(5*t))/2);

    calc_bone_transform(f, &f->skeleton.joints[0], 0,0);
    calc_bone_transform(f, &f->skeleton.joints[2], 0,0);
    calc_bone_transform(f, &f->skeleton.joints[3], 0,0);
    calc_bone_transform(f, &f->skeleton.joints[1], 0,0);
    calc_bone_transform(f, &f->skeleton.joints[10], 0,90);
    calc_bone_transform(f, &f->skeleton.joints[11], 0,90);
/*
    j = &f->skeleton.joints[0];
    update_skinning(f,j,0);
    j = &f->skeleton.joints[2];
    update_skinning(f,j,0);
    j = &f->skeleton.joints[3];
    update_skinning(f,j,1);
    j = &f->skeleton.joints[1];
    update_skinning(f,j,0);
    j = &f->skeleton.joints[10];
    update_skinning(f,j,0);
    j = &f->skeleton.joints[11];
    update_skinning(f,j,1);

    f->body->vertexArray = f->verts;
    glBindVertexArray(f->body->vao);
    glBindBuffer(GL_ARRAY_BUFFER, f->body->vb);
    glBufferData(GL_ARRAY_BUFFER, f->body->numVertices*3*sizeof(GLfloat), f->body->vertexArray, GL_STATIC_DRAW);
*/
  }

  else if(!f->animate)
  {

  int i;
  for(i=0;i<16;i++)
  {
    joint_s * j = &f->skeleton.joints[i];
    j->T = T(j->pos.x, j->pos.y, j->pos.z);
  }

  //printf("%f \n", f->skeleton.joints[0].force.z);

  f->pos.x = f->skeleton.joints[0].pos.x;
  f->pos.y = f->skeleton.joints[0].pos.y-6;
  f->pos.z = f->skeleton.joints[0].pos.z;
  mat4 test = Mult(T(0,6,0), Mult(f->body_R, T(0,-6,0)));
  f->matrix = Mult(T(f->pos.x, f->pos.y, f->pos.z), test);

  calc_bone_transform_ragdoll(f, &f->skeleton.joints[0], 0,0);
  //calc_bone_transform_ragdoll(f, &f->skeleton.joints[2], 0,0);
//  calc_bone_transform_ragdoll(f, &f->skeleton.joints[3], 0,0);
//  calc_bone_transform_ragdoll(f, &f->skeleton.joints[1], 0,0);
//  calc_bone_transform_ragdoll(f, &f->skeleton.joints[10], 0,90);
//  calc_bone_transform_ragdoll(f, &f->skeleton.joints[11], 0,90);



  //SKINNING HERE
  //mat4 R = Rx(M_PI/2.3);
/*
  joint_s * j = &f->skeleton.joints[2];
  update_skinning(f,j,-1);
  j = &f->skeleton.joints[3];
  update_skinning(f,j,-1);
  j = &f->skeleton.joints[10];
  update_skinning(f,j,-1);
  j = &f->skeleton.joints[11];
  update_skinning(f,j,-1);
  j = &f->skeleton.joints[1];
  update_skinning(f,j,-1);
  j = &f->skeleton.joints[0];
  update_skinning(f,j,-1);
*/

  joint_s * j = &f->skeleton.joints[0];
  update_skinning(f,j,-1);

  f->body->vertexArray = f->verts;
  glBindVertexArray(f->body->vao);
  glBindBuffer(GL_ARRAY_BUFFER, f->body->vb);
  glBufferData(GL_ARRAY_BUFFER, f->body->numVertices*3*sizeof(GLfloat), f->body->vertexArray, GL_STATIC_DRAW);
  }

}



void update_skinning(farmer_s * f, joint_s * j, int side)
{
  int jj=0;
  int i;
  vec3 midd;
  float len;
  while(j->child[0] != NULL)
  {
    //if(j->child[1] != NULL)
    //  update_skinning(f, j->child[1], side);
    //if(j->child[2] != NULL)
    //  update_skinning(f,j->child[2], side);


  for(i=0;i < f->body->numVertices;i++)
  {
    vec3 orig_vert = SetVector(f->orig_body->vertexArray[i*3+0], 
	f->orig_body->vertexArray[i*3+1],
	f->orig_body->vertexArray[i*3+2]);

    vec3 to = j->orig_pos;
    vec3 to2;
    to2 = j->child[0]->orig_pos;

    midd = ScalarMult(VectorAdd(to, to2), .5);

    len = fabs(Norm(VectorSub(to, midd)));

    vec3 dist = VectorSub(orig_vert, midd);

    int boolses;

    if(side == 0)
      boolses = Norm(dist) < len*j->weight && orig_vert.z > 0;
    else if(side == 1)
      boolses = Norm(dist) < len*j->weight && orig_vert.z < 0;
    else
      boolses = Norm(dist) < len*j->weight;

    if(boolses)
    {
      //vec3 res = MultVec3(Mult(T(to.x, to.y, to.z),
	//Mult(j->R,T(-to.x, -to.y, -to.z))), orig_vert);
      vec3 res = MultVec3(j->tmp, orig_vert);

      f->verts[i*3 + 0] = res.x;
      f->verts[i*3 + 1] = res.y;
      f->verts[i*3 + 2] = res.z;
    }

  }

  //printf("hej %f %f %f\n", midd.x, midd.y, midd.z);

  //printf("%i %f\n", jj, len*j->weight);

  j = j->child[0];
  //printf("%i\n", jj);
  jj++;
  }

}

void draw_farmer(farmer_s * f, GLuint program)
{
  glUniform1i(glGetUniformLocation(program, "draw_farmer"), 1);
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, f->matrix.m);
  glBindTexture(GL_TEXTURE_2D, f->tex);
  if(keyIsDown('m'))
  DrawModel(f->body, program, "inPosition", "inNormal", "inTexCoord");
  glUniform1i(glGetUniformLocation(program, "draw_farmer"), 0);


  int i;
  for(i=0;i<16;i++)
  {
    //mat4 tmp = Mult(f->skeleton.joints[i].T, S(.2,.2,.2));
    mat4 tmp = Mult(T(f->skeleton.joints[i].pos.x,
    	f->skeleton.joints[i].pos.y, f->skeleton.joints[i].pos.z), S(.15,.15,.15));
    //mat4 tmp = T(f->skeleton.joints[i].pos.x,
//	f->skeleton.joints[i].pos.y, f->skeleton.joints[i].pos.z);
    glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, tmp.m);
    DrawModel(f->skeleton.joints[i].body, program, "inPosition", "inNormal", "inTexCoord");
  }



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
  for(i=0;i<4;i++)
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
  c->debug_sphere = LoadModelPlus("./res/groundsphere.obj");
  c->matrix = S(.1,.1,.1); //makes it roughly 5 high
  c->pos = SetVector(5,0,5);
  c->speed = SetVector(0,0,0);
  c->acc = SetVector(0,0,0);
  c->force = SetVector(0,0,0);
  c->torque = SetVector(0,0,0);
  c->momentum = SetVector(0,0,0);
  c->mass = 2.0;
  c->angle = 8*M_PI;
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


  //update head pos
  c->head_pos = VectorAdd(c->pos, SetVector(-3,3,0));

  mat4 R = Mult(Mult(T(c->pos.x, c->pos.y, c->pos.z), 
  c->R), T(-c->pos.x, -c->pos.y, -c->pos.z));
  c->head_pos = MultVec3(R, c->head_pos);
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
  j->leader = 0;
  j->isnull = 1;
  j->weight = 1;
  j->Mvar = Mvar;
  j->posvar = posvar;
  j->boneposvar = boneposvar;

  j->orig_pos = pos;
  j->orig_vec = SetVector(0,0,1);
  j->constraint = 0;
  j->max_angle = 360;
  j->pos = pos;
  j->id = id;
  j->body = LoadModelPlus("./res/groundsphere.obj");
  //j->body_matrix = Mult(T(pos.x, pos.y, pos.z), S(.1, .1, .1));
  j->body_matrix = T(pos.x, pos.y, pos.z);
  j->T = T(pos.x, pos.y, pos.z);
  j->orig_T = T(j->orig_pos.x, j->orig_pos.y, j->orig_pos.z);
  j->M = T(pos.x, pos.y, pos.z);
  j->Mp = T(pos.x, pos.y, pos.z);
  j->Mtot = T(pos.x, pos.y, pos.z);
  //j->Minv = InvertMat4(j->M);
  j->Minv = InvertMat4(T(pos.x, pos.y, pos.z));
  j->R = IdentityMatrix();
  //j->body_matrix = S(.1, .1, .1);
  //j->body_matrix = IdentityMatrix();
}

void create_ragdoll_joint(joint_s * j, vec3 pos)
{
  j->pos = pos;
  j->orig_pos = pos;
  j->body = LoadModelPlus("./res/groundsphere.obj");
  j->body_matrix = Mult(T(pos.x, pos.y, pos.z), S(.02, .02, .02));
  j->T = T(pos.x, pos.y, pos.z);
  j->force = SetVector(0,0,0);
  j->calculated_force = SetVector(0,0,0);

}

void create_ragdoll(ragdoll_s * r)
{
  r->num_joints = 4;

  //head
  create_ragdoll_joint(&r->joints[0], SetVector(0,8,0));

  //neck and shoulders, arms
  create_ragdoll_joint(&r->joints[1], SetVector(0,10,0));
/*
  create_ragdoll_joint(&r->joints[4], SetVector(0,8,2));
  create_ragdoll_joint(&r->joints[5], SetVector(0,8,-2));
  create_ragdoll_joint(&r->joints[6], SetVector(0,8,4));
  create_ragdoll_joint(&r->joints[7], SetVector(0,8,-4));

*/
  //stomach
  create_ragdoll_joint(&r->joints[2], SetVector(0,8,-2));
  //groin
  create_ragdoll_joint(&r->joints[3], SetVector(0,10,-2));
/*
  //legs
  create_ragdoll_joint(&r->joints[8], SetVector(0,4, 2));
  create_ragdoll_joint(&r->joints[9], SetVector(0,4,-2));
  create_ragdoll_joint(&r->joints[10], SetVector(0,2,2));
  create_ragdoll_joint(&r->joints[11], SetVector(0,2,-2));
  create_ragdoll_joint(&r->joints[12], SetVector(0,0,2));
  create_ragdoll_joint(&r->joints[13], SetVector(0,0,-2));
*/

  r->joints[0].parent = NULL;
  r->joints[1].parent = &r->joints[0];
  r->joints[2].parent = &r->joints[1];
  r->joints[3].parent = &r->joints[2];
/*
  r->joints[4].parent = &r->joints[1];
  r->joints[5].parent = &r->joints[1];
  r->joints[6].parent = &r->joints[4];
  r->joints[7].parent = &r->joints[5];

  r->joints[8].parent = &r->joints[3];
  r->joints[9].parent = &r->joints[3];
  r->joints[10].parent = &r->joints[8];
  r->joints[11].parent = &r->joints[9];
  r->joints[12].parent = &r->joints[10];
  r->joints[13].parent = &r->joints[11];
*/

  r->joints[1].dist_to_parent = 2;
  r->joints[2].dist_to_parent = 2.5;
  r->joints[3].dist_to_parent = 2;
/*
  r->joints[4].dist_to_parent = 2;
  r->joints[5].dist_to_parent = 2;
  r->joints[6].dist_to_parent = 2;
  r->joints[7].dist_to_parent = 2;
  r->joints[8].dist_to_parent = 2;
  r->joints[9].dist_to_parent = 2;
  r->joints[10].dist_to_parent = 2;
  r->joints[11].dist_to_parent = 2;
  r->joints[12].dist_to_parent = 2;
  r->joints[13].dist_to_parent = 2;
*/
  //set a force on the head
  //r->joints[0].force = SetVector(0,0,.1);
}



void update_ragdoll(ragdoll_s * r, GLfloat dT)
{
  if(keyIsDown('a'))
    r->joints[0].force.z = 20;
  else if(keyIsDown('q'))
    r->joints[0].force.x = -80;
  else if(keyIsDown('o'))
    r->joints[0].force.z = -20;
  else if(keyIsDown('.'))
    r->joints[0].force.y = 80;
  else
  {
    //r->joints[0].force.z = 0;
    //r->joints[0].force.y = 0;
    //r->joints[0].force.x = 0;
  }

  int test=0;
  vec3 dP, dX, calculated_force = {0,-30,0};
  vec3 n;
  vec3 real_dist;
  int i;

  //update orig_vec for 1 and 8
  vec3 orig = Normalize(VectorSub(r->joints[0].pos, r->joints[8].pos));

  r->joints[0].orig_vec = orig;
  r->joints[1].orig_vec = orig;
  r->joints[8].orig_vec = orig;
  r->joints[9].orig_vec = orig;

  r->joints[12].orig_vec = orig;
  r->joints[13].orig_vec = orig;
  r->joints[14].orig_vec = orig;
  r->joints[15].orig_vec = orig;


  //printf("hejs %f %f %f\n", orig.x, orig.y, orig.z);

  for(i=0;i < r->num_joints;i++)
  {
    joint_s * j = &r->joints[i];
    joint_s * parent = r->joints[i].parent;

      if(j->pos.y <= 0)
      {
        //raise(SIGABRT);
        float ang;
        vec3 real_dist;
        j->pos.y += .01;
        joint_s * jp = j->parent;
        while(jp != NULL)
        {
          if(jp->parent != NULL)
          {
          real_dist = VectorSub(jp->parent->pos, jp->pos);
          ang = 180*acos(DotProduct(Normalize(real_dist), jp->orig_vec))/M_PI;
          }
          else ang=0;
          if(ang < jp->max_angle)
          {
          jp->pos.y += .01;
          jp->speed.y *= -.8;
          }
          jp = jp->parent;


        }

        j->speed.y *= -.8;
      }

    //if(i==8)
    //  printf("%f\n", r->joints[8].force.y);

    dP = ScalarMult(VectorAdd(j->force, calculated_force), dT);
    j->speed = VectorAdd(j->speed, dP);
    dX = ScalarMult(j->speed, dT);
    j->pos = VectorAdd(j->pos, dX);
    //parent = parent->parent;

    if(parent != NULL && j != NULL)
    {

      //parent = parent->parent;
      real_dist = VectorSub(j->pos, parent->pos);
      n = Normalize(real_dist);
      float dist_diff = (Norm(real_dist)-j->dist_to_parent);

      //if(i==1)
      //  printf("%f %f %f \n", n.x, n.y, n.z);

      vec3 speed_diff = VectorSub(j->speed, parent->speed);
      //calculated_force = VectorSub(ScalarMult(n, -dist_diff*40),
      //ScalarMult(speed_diff, 5));
      //calculated_force = VectorSub(ScalarMult(n, -dist_diff*10),
      //ScalarMult(speed_diff, 5));
      //calculated_force = SetVector(0,-10,0);
      j->pos = VectorSub(j->pos, ScalarMult(n, dist_diff));

      //j->pos = VectorSub(j->pos, ScalarMult(n, dist_diff));

      float ang = 180*acos(DotProduct(Normalize(real_dist), j->orig_vec))/M_PI;


      while((ang > j->max_angle && j->constraint) && j->pos.y > 0)
      {
        //printf("%f\n", ang);

        real_dist = VectorSub(parent->pos, j->pos);
        dist_diff = (Norm(real_dist) - j->dist_to_parent);

        n = Normalize(real_dist);
        vec3 align_vec = Normalize(VectorSub(real_dist, j->orig_vec));

        if(dist_diff < 0)
          j->pos = VectorAdd(j->pos, ScalarMult(n, -.01));
        else
          j->pos = VectorAdd(j->pos, ScalarMult(n, .01));

        j->pos = VectorAdd(j->pos, ScalarMult(align_vec, .01));
        ang = 180*acos(DotProduct(Normalize(real_dist), j->orig_vec))/M_PI;

	//if(i==4)
	//  printf("hej %f %f %f , %f %f\n", real_dist.x, real_dist.y, real_dist.z, j->dist_to_parent, Norm(real_dist));

      }

    j->T = Mult(T(j->pos.x, j->pos.y, j->pos.z), S(.5,.5,.5));


    }

    else if(j->leader)
    {
      
    }


  if(j->pos.y <= .3)
    j->speed = ScalarMult(j->speed, .95);
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


void create_fence(fence_s * f, int width, vec3 pos)
{

  f->width = width;
  f->pos = pos;

  int i;
  for(i=0;i<width;i++)
  {
    //create upright (type 0)
    create_plank(&f->planks[i*3][0], 
	VectorAdd(pos, SetVector(i*10.8,0,0)), 0);

    create_plank(&f->planks[i*3+1][0],
	 VectorAdd(pos, SetVector(i*10.8,4,0)), 1);

    create_plank(&f->planks[i*3+2][0],
	 VectorAdd(pos, SetVector(i*10.8,2,0)), 1);
  }

  create_plank(&f->planks[(width-1)*3+3][0],
	VectorAdd(pos, SetVector(width*10.8,0,0)), 0);

}

void draw_fence(fence_s * f, GLuint program)
{

  int i;
  for(i=0;i < f->width*3+1;i++)
  {
    if(f->planks[i][0].destroyed_at == -1 || f->planks[i][0].type == 2 || f->planks[i][0].type == 3)
      draw_plank(&f->planks[i][0], program);
    if(f->planks[i][1].destroyed_at == -1 || f->planks[i][1].type == 2 || f->planks[i][1].type == 3)
      draw_plank(&f->planks[i][1], program);
    //glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, f->planks[i].T.m);
    //DrawModel(f->planks[i].body, program, "inPosition", "inNormal", "inTexCoord");
  }


  //draw_plank(&f->planks[0], program);

}


void create_plank(plank_s * p, vec3 pos, int type)
{

  int i,j,x,y,z;
  //int tri_count = (size.x-1)*(size.y-1)*(size.z-1)*3;

  p->pos = pos;
  p->type = type;

  p->mass = 2*(float)rand()/(float)RAND_MAX;
  printf("%f\n", p->mass);

  p->force = SetVector(0,-15,0);
  //not destroyed
  if(type==0 || type==1)
    p->destroyed_at = -1;

  if(type==0)
    //p->T = Mult(Mult(T(pos.x, pos.y, pos.z), Rz(M_PI/2)), T(-pos.x, -pos.y, -pos.z));
    //p->T = Mult( T(pos.x, pos.y, pos.z), Rz(M_PI/2) );
    p->R = Rz(M_PI/2);
    //p->T = Rz(M_PI/2);
    //p->T = T(0,0,0);
  else if(type==2)
    //p->T = Mult(Mult(T(pos.x, pos.y, pos.z), Ry(M_PI/1.1)), T(-p->destroyed_at*5/3.0,0,-.5));
    p->R = Mult(Ry(M_PI/1.1), T(-p->destroyed_at*5/6.0,0,-.5));
  //else
  else
    p->R = IdentityMatrix();

  p->T = T(pos.x, pos.y, pos.z);

  p->M = Mult(p->T, p->R);

  int fine=1;

#define XS 12
#define YS 10
#define ZS 2

GLfloat verts[XS*YS*ZS*3*2];
GLfloat norms[XS*YS*ZS*3*2];
i=0;
  for(x=0;x < XS;x++)
  {
    for(y = 0;y < YS;y++)
    {
      for(z = 0;z < ZS;z++)
      {
        //printf("%i %i %i\n", x,y,z);
        //if upright
        if(type==0)
          verts[i*3 + 0] = x*.5;
        else if(type==1)
          verts[i*3 + 0] = x*1.0;
        else if(type==2)
          verts[i*3 + 0] = (p->destroyed_at)*x/12.0;
        else
        {
          verts[i*3 + 0] = (12-p->destroyed_at)*x/12.0;
          //printf("BAAAJS\n");
        }


        verts[i*3 + 1] = y/10.0;
        verts[i*3 + 2] = z/2.0;

	i++;
      }
    }

  }


if(type==3)
{
verts[0] = -.5;
verts[3] = -.2;

verts[12] = -.9;
verts[15] = -.3;

verts[24] = -1.1;
verts[27] = -.3;

verts[36] = -0.4;
verts[39] = -.9;

verts[48] = -.2;
verts[51] = -.4;
}
else if(type==2)
{
verts[0] = -1;
verts[3] = -.2;

verts[12] = -.9;
verts[15] = -.3;

verts[24] = -.5;
verts[27] = -.3;

verts[36] = -.4;
verts[39] = -.9;

verts[48] = -.2;
verts[51] = -1.2;



}

Point3D tmp_normal;
i=0;
  for(x=0;x < XS;x++)
  {
    for(y = 0;y < YS;y++)
    {
      for(z = 0;z < ZS;z++)
      {

      //calc_normal(verts, x, y, YS, &tmp_normal);
      if(z>0)
        norms[i*3 + 2] = 1;
      else
        norms[i*3 + 2] = -1;
      norms[i*3 + 1] = 0;
      norms[i*3 + 0] = 0;
      i++;
      }
    }
  }


GLuint indic[XS*YS*ZS*3*3] = {0};
//int j;
i=0;
j=0;
//printf("\n\n");
//printf("%i\n\n", y);
for(i=0;i < XS-1;i++)
{
  /*
  indic[i*9 + 0] = i*8;
  indic[i*9 + 1] = i*8+2;
  indic[i*9 + 2] = i*8 + 6+2;
  indic[i*9  + 3] = i*8 + 2;
  indic[i*9  + 4] = i*8+ 2 + 2;
  indic[i*9  + 5] = i*8 + 2 + 6+2;
  indic[i*9  + 6] = i*8 + 2+2;
  indic[i*9  + 7] = i*8+ 2 + 2+2;
  indic[i*9  + 8] = i*8 + 2 + 6+2+2;
  */

  for(j=0;j < YS-1;j++)
  {
    //one side (z=0)
    indic[i*(y-1)*3 + 0+ 3*j] = (y)*2*i + 2*j;
    indic[i*(y-1)*3 + 1+ 3*j] = (y)*2*i + 2*j + 2;
    indic[i*(y-1)*3 + 2+ 3*j] = (y)*2*i + 2*j + (y)*2;

    indic[i*(y-1)*3 + 0+ 3*j + (XS-1)*(YS-1)*3] = (y)*2*i + 2*j + 2;
    indic[i*(y-1)*3 + 1+ 3*j + (XS-1)*(YS-1)*3] = (y)*2*i + 2*j +(y)*2;
    indic[i*(y-1)*3 + 2+ 3*j + (XS-1)*(YS-1)*3] = (y)*2*i + 2*j + (y)*2 +2;

    //other side (z=1)
    indic[i*(y-1)*3 + 0+ 3*j + (XS-1)*(YS-1)*3*2] = (y)*2*i + 2*j + 1;
    indic[i*(y-1)*3 + 1+ 3*j + (XS-1)*(YS-1)*3*2] = (y)*2*i + 2*j + 2 + 1;
    indic[i*(y-1)*3 + 2+ 3*j + (XS-1)*(YS-1)*3*2] = (y)*2*i + 2*j + (y)*2 + 1;

    indic[i*(y-1)*3 + 0+ 3*j + (XS-1)*(YS-1)*3*3] = (y)*2*i + 2*j + 2 + 1;
    indic[i*(y-1)*3 + 1+ 3*j + (XS-1)*(YS-1)*3*3] = (y)*2*i + 2*j +(y)*2 + 1;
    indic[i*(y-1)*3 + 2+ 3*j + (XS-1)*(YS-1)*3*3] = (y)*2*i + 2*j + (y)*2 +2 + 1;

  }


  //bottom
  indic[(XS-1)*(YS-1)*3*4 + i*6 + 0] = y*2*i;
  indic[(XS-1)*(YS-1)*3*4 + i*6 + 1] = y*2*i + 1;
  indic[(XS-1)*(YS-1)*3*4 + i*6 + 2] = y*2*i + 2*(y);

  indic[(XS-1)*(YS-1)*3*4 + i*6 + 3] = y*2*i + 1;
  indic[(XS-1)*(YS-1)*3*4 + i*6 + 4] = y*2*i + y*2;
  indic[(XS-1)*(YS-1)*3*4 + i*6 + 5] = y*2*i + 2*(y) + 1;


  //top
  indic[(XS-1)*(YS-1)*3*5  + i*6 + 0] = y*2*i + (YS-1)*2;
  indic[(XS-1)*(YS-1)*3*5  + i*6 + 1] = y*2*i + 1 + (YS-1)*2;
  indic[(XS-1)*(YS-1)*3*5  + i*6 + 2] = y*2*i + 2*(y) + (YS-1)*2;

  indic[(XS-1)*(YS-1)*3*5  + i*6 + 3] = y*2*i + (YS-1)*2+1;
  indic[(XS-1)*(YS-1)*3*5  + i*6 + 4] = y*2*i + 1 + (YS-1)*4+1;
  indic[(XS-1)*(YS-1)*3*5  + i*6 + 5] = y*2*i + 3 + (YS-1)*4;

}

for(i=0;i<2*YS;i++)
{
  //left side
  indic[(XS-1)*(YS-1)*3*6  + i*6 + 0] = i + 0;
  indic[(XS-1)*(YS-1)*3*6  + i*6 + 1] = i + 1;
  indic[(XS-1)*(YS-1)*3*6  + i*6 + 2] = i + 2;

  //right side
  indic[(XS-1)*(YS-1)*3*6  + i*6 + 3] = i + 0 + 2*(YS)*(XS-1);
  indic[(XS-1)*(YS-1)*3*6  + i*6 + 4] = i + 1 + 2*(YS)*(XS-1);
  indic[(XS-1)*(YS-1)*3*6  + i*6 + 5] = i + 2 + 2*(YS)*(XS-1);


//  indic[72*2+18+18 + i*6 + 0] = i + 0;
//  indic[72*2+18+18 + i*6 + 1] = i + 1;
//  indic[72*2+18+18 + i*6 + 2] = i + 2;

/*
  //right side
  indic[72*2+18+18 + i*6 + 3] = i + 0+30;
  indic[72*2+18+18 + i*6 + 4] = i + 1+30;
  indic[72*2+18+18 + i*6 + 5] = i + 2+30;
*/
}


for(i=0;i < 32*6*4;i++)
{
  //printf("%i ", indic[i]);
  //if((i+1)%3==0)
  //  printf("\n");
}

  p->body = LoadDataToModel(
	verts,
	norms,
	NULL,
	NULL,
	indic,
	(XS)*(YS)*ZS,
	(XS)*(YS)*2*2*2*2-40);

  p->debug_sphere = LoadModelPlus("./res/groundsphere.obj");
}

void draw_plank(plank_s * p, GLuint program)
{
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, p->M.m);
  DrawModel(p->body, program, "inPosition", "inNormal", "inTexCoord");

  int i;
  for(i=0;i < 6;i++)
  {
/*
  if(p->type==1)
  {
  mat4 tmp = T(p->pos.x + 2*i, p->pos.y, p->pos.z);
  glUniformMatrix4fv(glGetUniformLocation(program, "mdl_matrix"), 1, GL_TRUE, tmp.m);
  DrawModel(p->debug_sphere, program, "inPosition", "inNormal", "inTexCoord");
  }
*/
  }
}


void update_fence(fence_s * f, cow_s *c, GLfloat dT)
{
  int i,j;
  for(j=1;j < 40;j+=1)
  {
  float dist = Norm(VectorSub(c->head_pos, f->planks[j][0].pos));
  //first check
  if(dist < 10)
  {
    //printf("COLLISION %f\n", dist);
    //alright, now test for every vertex along x axis
    for(i=0;i < 12;i++)
    {
      dist = Norm(VectorSub(c->head_pos, 
	VectorAdd(f->planks[j][0].pos, SetVector(i*1, 0, 0))));
      if(dist < 1.5)
      {
        //set plank to state "destroyed"
        printf("COLL2, %f %i \n", dist, i);
        if(f->planks[j][0].destroyed_at == -1)
        {
          f->planks[j][0].destroyed_at = i;
          f->planks[j+1][0].destroyed_at = i;
          f->planks[j][1].destroyed_at = i;
          f->planks[j+1][1].destroyed_at = i;

          create_plank(&f->planks[j][0], VectorAdd(f->pos, SetVector((j-1)*3.6, 4, 0)), 2);
          create_plank(&f->planks[j][1], VectorAdd(f->pos, SetVector((j-1)*3.6+i*1.8/2, 4, 0)), 3);
          create_plank(&f->planks[j+1][0], VectorAdd(f->pos, SetVector((j-1)*3.6, 2, 0)), 2);
          create_plank(&f->planks[j+1][1], VectorAdd(f->pos, SetVector((j-1)*3.6+i*1.8/2, 2, 0)), 3);


          f->planks[j][0].momentum = c->momentum;
          f->planks[j][1].momentum = c->momentum;

          f->planks[j+1][0].momentum = c->momentum;
          f->planks[j+1][1].momentum = c->momentum;

          f->planks[j][0].angular_momentum = SetVector(.8,2,.2);
          f->planks[j+1][0].angular_momentum = SetVector(.8,2,.2);
          f->planks[j][1].angular_momentum = SetVector(.8,2,.2);
          f->planks[j+1][1].angular_momentum = SetVector(.8,2,.2);


          c->momentum = ScalarMult(c->momentum, -1);
          c->pos = VectorAdd(c->pos, SetVector(0,0,.1));


          //printf("destroyed at %i\n", i);
        }
        //break;

      }
    }

  }

  }

  for(i=0;i < 20;i++)
  {
    if(f->planks[i][0].type == 2 || f->planks[i][0].type==3)
    {
    update_plank(&f->planks[i][0], c, dT);
    update_plank(&f->planks[i][1], c, dT);
    }
  }

}

void update_plank(plank_s * p, cow_s * c, GLfloat dT)
{
  //p->pos = VectorAdd(p->pos, SetVector(0,0,-dT*20));

  vec3 dP = ScalarMult(p->force, dT);
  vec3 dL = ScalarMult(p->torque, dT);
  p->angular_momentum = VectorAdd(p->angular_momentum, dL);
  p->omega = p->angular_momentum;
  vec3 dO = ScalarMult(p->omega, dT);
  mat4 Rd = CrossMatrix(dO);
  Rd = Mult(Rd, p->R);
  p->R = MatrixAdd(p->R, Rd);

  p->momentum = VectorAdd(p->momentum, dP);
  p->speed = ScalarMult(p->momentum, 1.0/p->mass);
  vec3 dX = ScalarMult(p->speed, dT);

  p->pos = VectorAdd(p->pos, dX);
  p->M = Mult(T(p->pos.x, p->pos.y, p->pos.z), p->R);

  OrthoNormalizeMatrix(&p->R);
}

void reset_farmer_matrices(farmer_s * f)
{
  int i;
  for(i=0;i < f->skeleton.num_joints;i++)
  {
    f->skeleton.joints[i].T = IdentityMatrix();
    f->skeleton.joints[i].R = IdentityMatrix();
    f->skeleton.joints[i].tmp = IdentityMatrix();


  }

}
