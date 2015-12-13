#version 150

in vec4 colorr;
in vec2 outTexCoord;
in vec4 inPos;
in vec3 ex_normal;
uniform sampler2D tex;
out vec4 fragColor;

uniform int draw_cow;
uniform int draw_floor;
uniform int draw_plank;
uniform int draw_ball;
uniform int collision;

uniform vec3 cow_pos;
uniform mat4 cam_matrix;

vec3 s;
vec3 n;

mat3 light_cam_matrix = mat3(cam_matrix);

vec3 colors;
vec3 reflection;

void main(void)
{
  colors = vec3(0,0,0);
  s = normalize(light_cam_matrix * vec3(inPos));
  n = normalize(ex_normal);
  float lambert = dot(n,s)-.001;  

  if(lambert > 0)
  {
    colors += vec3(1,1,1)*lambert;
    reflection = reflect(s,n);
    
  }

  if(draw_cow==1)
    fragColor = (vec4(.2,.2,.2,0) +  vec4(colors,1))*texture(tex, outTexCoord*1);
  else if(draw_floor==1)
    fragColor = texture(tex, outTexCoord*1);
  else if(draw_plank==1)
    fragColor = vec4(colors/6,0) + vec4(.25,.1,0, 1);
  else
    fragColor = texture(tex, outTexCoord)*vec4(colors+.4, 1);

  if(length(vec3(inPos)-cow_pos) < 2.5)
    fragColor += vec4(-.2, -.2, -.2, 0);

  //fragColor = vec4(colors,1);

  //if(collision==1 && draw_cow==1)
  //if(collision==1)
  //  fragColor = vec4(.5,0,.5,1);

}
