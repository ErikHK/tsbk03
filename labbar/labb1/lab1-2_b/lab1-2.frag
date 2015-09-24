#version 150
// bump mapping should be calculated
// 1) in view coordinates
// 2) in texture coordinates

in vec2 outTexCoord;
in vec3 out_Normal;
in vec3 Ps;
in vec3 Pt;
in vec3 pixPos;  // Needed for specular reflections
uniform sampler2D texUnit;
out vec4 out_Color;

void main(void)
{
	vec3 light = vec3(0.0, 0.7, 0.7); // Light source in view coordinates
	vec3 normal = normalize(out_Normal);
	mat3 Mvt = transpose(mat3(Ps, Pt, normal));

	
	// Calculate gradients here
	float offset = 1.0 / 256.0; // texture size, same in both directions
	vec3 dx = vec3(texture(texUnit, vec2(outTexCoord.x + offset, outTexCoord.y)) - texture(texUnit, outTexCoord));
	vec3 dy = vec3(texture(texUnit, vec2(outTexCoord.x, outTexCoord.y+offset)) - texture(texUnit, outTexCoord));

	normal.x += dx.x;
	normal.y += dy.y;
	normal = normalize(normal);
	
	// Simplified lighting calculation.
	// A full solution would include material, ambient, specular, light sources, multiply by texture.
	out_Color = vec4(dot(normal, Mvt*light));
}
