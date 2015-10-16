#version 150

in vec4 colorr;
in vec2 outTexCoord;
uniform sampler2D tex;
out vec4 fragColor;

uniform int draw_cow;

void main(void)
{
	if(draw_cow)
	  fragColor = colorr;
	else
	  fragColor = colorr*texture(tex, outTexCoord*10);
}
