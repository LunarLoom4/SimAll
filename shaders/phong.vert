#version 450 core
// SimAll Beta — Blinn-Phong vertex stage.
layout(location = 0) in vec3 a_position;
layout(location = 1) in vec3 a_normal;
layout(location = 2) in vec2 a_uv;
layout(location = 3) in float a_scalar;

layout(std140, binding = 0) uniform Camera {
    mat4 uView;
    mat4 uProj;
    mat4 uViewProj;
    vec3 uEyeWorld;
    float _pad0;
};

layout(std140, binding = 1) uniform Model {
    mat4 uModel;
    mat4 uNormalMat;     // transpose(inverse(model))
};

out VS_OUT {
    vec3 posWorld;
    vec3 normalWorld;
    vec2 uv;
    float scalar;
} vs_out;

void main() {
    vec4 wp = uModel * vec4(a_position, 1.0);
    vs_out.posWorld    = wp.xyz;
    vs_out.normalWorld = normalize(mat3(uNormalMat) * a_normal);
    vs_out.uv          = a_uv;
    vs_out.scalar      = a_scalar;
    gl_Position        = uViewProj * wp;
}
