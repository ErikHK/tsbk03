#version 150

in vec4 colorr;
in vec2 outTexCoord;
in vec4 inPos;
uniform sampler2D tex;
out vec4 fragColor;

uniform int draw_cow;
uniform int draw_floor;
uniform int draw_ball;

uniform vec3 cow_pos;

void main(void)
{
	if(draw_cow)
	  fragColor = colorr*texture(tex, outTexCoord*1);
	else if(draw_floor)
	  //fragColor = colorr*texture(tex, outTexCoord*5)*.1 + vec4(.4,.6,.4,1);
	  fragColor = vec4(sin(inPos.x), sin(inPos.y), sin(inPos.z), 1);
	else
	  fragColor = colorr;

	if(length(vec3(inPos)-cow_pos) < 2.5)
	  fragColor += vec4(-.2, -.2, -.2, 0);
}
