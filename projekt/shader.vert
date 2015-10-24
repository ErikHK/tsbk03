#version 150

//in vec3 in_Color;
in vec3 inPosition;
in vec3 inNormal;
in vec2 inTexCoord;

out vec4 colorr;
out vec2 outTexCoord;
out vec4 inPos;
uniform mat4 mdl_matrix;
uniform mat4 cam_matrix;
uniform mat4 proj_matrix;

uniform float time;
uniform int draw_cow;
uniform int draw_floor;
uniform int draw_ball;
uniform int collision;

uniform mat4 foot_joint[4];
uniform mat4 testjoint[8];

uniform mat4 legjoint0[8];
uniform mat4 legjoint1[8];
uniform mat4 legjoint2[8];
uniform mat4 legjoint3[8];

uniform mat4 testjoint2;

uniform vec3 currpos[8];
uniform vec3 bonepos[8];

uniform vec3 legcurrpos0[8];
uniform vec3 legcurrpos1[8];
uniform vec3 legcurrpos2[8];
uniform vec3 legcurrpos3[8];

uniform vec3 legbonepos0[8];
uniform vec3 legbonepos1[8];
uniform vec3 legbonepos2[8];
uniform vec3 legbonepos3[8];

uniform vec3 tailcurrpos[8];
uniform vec3 tailbonepos[8];
uniform mat4 tailjoint[8];

uniform vec3 headcurrpos[8];
uniform vec3 headbonepos[8];
uniform mat4 headjoint[8];

uniform vec3 bodycurrpos[8];
uniform vec3 bodybonepos[8];
uniform mat4 bodyjoint[8];

uniform vec3 rightearcurrpos[8];
uniform vec3 rightearbonepos[8];
uniform mat4 rightearjoint[8];


mat4 rotationMatrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}

mat4 translationMatrix(vec3 pos)
{
return mat4(1,0,0,pos.x,
	    0,1,0,pos.y,
	    0,0,1,pos.z,
	    0,0,0,1);
}


//const vec3 lightDir = normalize(vec3(0.3, 0.5, 1.0));

void main(void)
{
  vec3 lightv = vec3(.6, .6, .6);

  float shade = dot(lightv, inNormal);

  vec3 lightv2 = vec3(-.6, -.6, -.6);

  float shade2 = dot(lightv2, inNormal);

  if(draw_cow==1)
  {
  //calc distance to the bones
  float dist, leng;

  //head
  for(int i=0;i<4;i++)
  {
    //length of bone (joint to middle of bone times 2)
    leng = 2*distance(headcurrpos[i], headbonepos[i]);

    //distance to middle of bone
    dist = distance(vec3(inPosition), headbonepos[i]);

    if(dist < leng*1.1)
      gl_Position += proj_matrix*cam_matrix*headjoint[i]*mdl_matrix*vec4(inPosition, 1);
  }

  //body
  for(int i=0;i<7;i++)
  {
    //length of bone (joint to middle of bone times 2)
    leng = 2*distance(bodycurrpos[i], bodybonepos[i]);

    //distance to middle of bone
    dist = distance(vec3(inPosition), bodybonepos[i]);

    if(dist < leng*2.2 && i<=1)
      gl_Position += proj_matrix*cam_matrix*bodyjoint[i]*mdl_matrix*vec4(inPosition, 1);

    else if(dist < leng/2 && i<=2)
      gl_Position += proj_matrix*cam_matrix*bodyjoint[i]*mdl_matrix*vec4(inPosition, 1);

    else if(dist < leng*1.3 && i>3)
      gl_Position += proj_matrix*cam_matrix*bodyjoint[i]*mdl_matrix*vec4(inPosition, 1);
  }

  

  //leg 0
  for(int i=0;i<4;i++)
  {
    //length of bone (joint to joint)
    leng = 2*distance(legcurrpos0[i], legbonepos0[i]);

    //distance to middle of bone
    dist = distance(vec3(inPosition), legbonepos0[i]);

    if(abs(dist) < leng/1.5 && i<2)
      gl_Position += proj_matrix*cam_matrix*legjoint0[i]*mdl_matrix*vec4(inPosition, 1);
    else if(dist < leng/1.5)
      gl_Position += proj_matrix*cam_matrix*legjoint0[i]*mdl_matrix*vec4(inPosition, 1);

  }

  //leg 1
  for(int i=0;i<4;i++)
  {
    //length of bone (joint to joint)
    leng = 2*distance(legcurrpos1[i], legbonepos1[i]);

    //distance to middle of bone
    dist = distance(vec3(inPosition), legbonepos1[i]);

    if(abs(dist) < leng/2+1.8)
      gl_Position += proj_matrix*cam_matrix*legjoint1[i]*mdl_matrix*vec4(inPosition, 1);
  }


  //leg 2 (hind)
  for(int i=0;i<4;i++)
  {
    //length of bone (joint to joint)
    leng = 2*distance(legcurrpos2[i], legbonepos2[i]);

    //distance to middle of bone
    dist = distance(vec3(inPosition), legbonepos2[i]);

    if(abs(dist) < leng/2+1.7)
      gl_Position += proj_matrix*cam_matrix*legjoint2[i]*mdl_matrix*vec4(inPosition, 1);
  }

  //leg 3 (hind)
  for(int i=0;i<4;i++)
  {
    //length of bone (joint to joint)
    leng = 2*distance(legcurrpos3[i], legbonepos3[i]);

    //distance to middle of bone
    dist = distance(vec3(inPosition), legbonepos3[i]);

    if(abs(dist) < leng/2+1.7)
      gl_Position += proj_matrix*cam_matrix*legjoint3[i]*mdl_matrix*vec4(inPosition, 1);
  }

  //tail
  for(int i=0;i<5;i++)
  {
    //length of bone (joint to joint)
    leng = 2*distance(tailcurrpos[i], tailbonepos[i]);

    //distance to middle of bone
    dist = distance(vec3(inPosition), tailbonepos[i]);

    if(dist < leng*2)
      gl_Position += proj_matrix*cam_matrix*tailjoint[i]*mdl_matrix*vec4(inPosition, 1);
  }



  //right ear
  for(int i=0;i<2;i++)
  {
    //length of bone (joint to middle of bone times 2)
    leng = 2*distance(rightearcurrpos[i], rightearbonepos[i]);

    //distance to middle of bone
    dist = distance(vec3(inPosition), rightearbonepos[i]);

    if(dist < leng*3)
      gl_Position += proj_matrix*cam_matrix*rightearjoint[i]*mdl_matrix*vec4(inPosition, 1);
  }



  }

  else
    gl_Position = proj_matrix*cam_matrix*mdl_matrix*vec4(inPosition,1);


  if(draw_cow == 1)
    //colorr = vec4(shade, shade, shade, .6);
    colorr = vec4(.7,.7,.7,1) + (vec4(shade, shade, shade, 1))*.3;// + vec4(shade2,shade2,shade2,1)*1);
  else if(draw_ball == 1)
    //colorr = vec4(shade,shade,shade,1)*.6+vec4(shade2,shade2,shade2,1)*.4;
    colorr = vec4(shade,shade,shade,1)*.6;
  else if(draw_floor == 1)
    colorr = vec4(50*sin(inPosition.x*200), 40*cos(inPosition.y*400), 20*sin(inPosition.z*250), 1);
  else if(collision == 1)
    colorr = vec4(shade, shade, shade, 1) + vec4(.6,.3,.3,0);
  else
    colorr = vec4(shade, shade, shade, 1);
    
  outTexCoord = inTexCoord;
  inPos = vec4(inPosition,1);
}

