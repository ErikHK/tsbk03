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
  //float len = length(inPosition - vec3(0,1,0));
/*
  float len = distance(inPosition, vec3(0,8,0));
  if(len < 11)
    gl_Position = proj_matrix * cam_matrix * mdl_matrix * transpose(translationMatrix(vec3(0,20,0))) * rotationMatrix(vec3(0,0,5), sin(time*5)/5) * vec4(inPosition + vec3(0,-20,0), 1.0);
  else if(distance(inPosition, vec3(-30,8,0)) < 11)
    gl_Position = proj_matrix * cam_matrix * mdl_matrix * transpose(translationMatrix(vec3(-30,20,0))) * rotationMatrix(vec3(0,0,5), cos(time*5)/5) * vec4(inPosition + vec3(30,-20,0), 1.0);
  else
*/
    gl_Position = proj_matrix * cam_matrix * mdl_matrix * vec4(inPosition, 1.0);



  //gl_Position = vec4(in_Position, 1.0);
  if(draw_cow == 1)
    colorr = vec4(shade, shade, shade, .6) + vec4(.2,0,.2,0);
  else
    colorr = vec4(1);
    
  outTexCoord = inTexCoord;
}

