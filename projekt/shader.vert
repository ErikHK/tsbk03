#version 150

//in vec3 in_Color;
in vec3 inPosition;
in vec3 inNormal;
in vec2 inTexCoord;

out vec4 colorr;
out vec2 outTexCoord;
uniform mat4 mdl_matrix;
uniform mat4 cam_matrix;
uniform mat4 proj_matrix;

uniform float time;
uniform int draw_cow;

uniform mat4 foot_joint[4];
uniform mat4 testjoint[8];
uniform mat4 legtestjoint[8];
uniform mat4 testjoint2;

uniform vec3 currpos[8];
uniform vec3 bonepos[8];

uniform vec3 legcurrpos[8];
uniform vec3 legbonepos[8];

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

  if(draw_cow==1)
  {
  //calc distance to the bones
  float dist, dist2, len1, len2;

  for(int i=0;i<4;i++)
  {
    //length of bone (joint to joint)
    float leng = 2*distance(currpos[i], bonepos[i]);

    //distance to middle of bone
    dist = distance(vec3(inPosition), bonepos[i]);

    if(abs(dist) < leng/2+9)
      gl_Position += proj_matrix*cam_matrix*testjoint[i]*mdl_matrix*vec4(inPosition, 1);
  }


  for(int i=0;i<4;i++)
  {
    //length of bone (joint to joint)
    float leng = 2*distance(legcurrpos[i], legbonepos[i]);

    //distance to middle of bone
    dist = distance(vec3(inPosition), legbonepos[i]);

    if(abs(dist) < leng/2+2)
      gl_Position += proj_matrix*cam_matrix*legtestjoint[i]*mdl_matrix*vec4(inPosition, 1);
  }

  }

  else
    gl_Position = proj_matrix*cam_matrix*mdl_matrix*vec4(inPosition,1);


  if(draw_cow == 1)
    //colorr = vec4(shade, shade, shade, .6);
    colorr = vec4(shade, shade, shade, 1) + vec4(.4,.2,.4,0);
  else
    colorr = vec4(1);
    
  outTexCoord = inTexCoord;
}

