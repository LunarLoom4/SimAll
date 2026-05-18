#version 450 core
// SimAll Beta — volume ray-cast proxy-cube vertex stage.
// The proxy is a unit cube spanning uBoxMin..uBoxMax in world space.
// We emit per-vertex the world-space ray origin (eye) and direction
// (vertex - eye); the fragment shader then performs analytic
// ray/box intersection in world coordinates.

layout(location = 0) in vec3 a_position;

layout(std140, binding = 0) uniform Camera {
    mat4 uView;
    mat4 uProj;
    mat4 uViewProj;
    vec3 uEyeWorld;
    float _pad0;
};

layout(std140, binding = 7) uniform Volume {
    mat4  uInvModelViewProj;
    vec3  uEyeWorld_v;     float uStepSize;
    vec3  uBoxMin;         float uMinScalar;
    vec3  uBoxMax;         float uMaxScalar;
    float uMaxSteps;       float uJitter;
    float uDensity;        float _pad;
};

out vec3 vRayOrigin;
out vec3 vRayDir;

void main() {
    // Map a_position in [0,1]^3 onto the world-space box.
    vec3 wp = mix(uBoxMin, uBoxMax, a_position);
    vRayOrigin = uEyeWorld;
    vRayDir    = wp - uEyeWorld;
    gl_Position = uViewProj * vec4(wp, 1.0);
}
