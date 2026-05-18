#version 450 core
// SimAll Beta — emits barycentric coordinates so a fragment shader
// can draw exact-thickness wireframe over filled triangles.

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

in VS_OUT {
    vec3 posWorld;
    vec3 normalWorld;
    vec2 uv;
    float scalar;
} gs_in[];

out GS_OUT {
    vec3 posWorld;
    vec3 normalWorld;
    vec2 uv;
    float scalar;
    vec3 bary;       // barycentric weights of this vertex
} gs_out;

const vec3 kBary[3] = vec3[3](
    vec3(1.0, 0.0, 0.0),
    vec3(0.0, 1.0, 0.0),
    vec3(0.0, 0.0, 1.0)
);

void main() {
    for (int i = 0; i < 3; ++i) {
        gl_Position           = gl_in[i].gl_Position;
        gs_out.posWorld       = gs_in[i].posWorld;
        gs_out.normalWorld    = gs_in[i].normalWorld;
        gs_out.uv             = gs_in[i].uv;
        gs_out.scalar         = gs_in[i].scalar;
        gs_out.bary           = kBary[i];
        EmitVertex();
    }
    EndPrimitive();
}
