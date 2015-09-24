#version 150

//in vec3 in_Color;
in vec3 in_Position;
in vec3 in_Normal;
in vec2 in_TexCoord;
uniform mat4 matrix;

out vec4 g_color;
const vec3 lightDir = normalize(vec3(0.3, 0.5, 1.0));

// Uppgift 3: Soft-skinning på GPU
//
// Flytta över din implementation av soft skinning från CPU-sidan
// till vertexshadern. Mer info finns på hemsidan.

uniform mat4 g_bones_0_rot;
uniform mat4 g_bones_1_rot;
uniform mat4 g_bones_0_pos;
uniform mat4 g_bones_1_pos;



void main(void)
{
vec4 M0_pos0 = in_TexCoord.x*g_bones_0_pos*g_bones_0_rot*inverse(g_bones_0_pos)*vec4(in_Position, 1.0);
vec4 M1_pos1 = in_TexCoord.y*g_bones_1_pos*g_bones_1_rot*inverse(g_bones_1_pos)*vec4(in_Position, 1.0);

vec4 res = M0_pos0 + M1_pos1;

	// transformera resultatet med ModelView- och Projection-matriserna
	gl_Position = matrix*res;

	// sätt röd+grön färgkanal till vertex Weights
	vec4 color = vec4(in_TexCoord.x, in_TexCoord.y, 0.0, 1.0);

	// Lägg på en enkel ljussättning på vertexarna 	
	float intensity = dot(in_Normal, lightDir);
	color.xyz *= intensity;

	g_color = color;
}

