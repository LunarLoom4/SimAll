#version 450 core
// SimAll Beta — Blinn-Phong fragment stage with 3 directional + 1 headlight.

in VS_OUT {
    vec3 posWorld;
    vec3 normalWorld;
    vec2 uv;
    float scalar;
} fs_in;

layout(std140, binding = 0) uniform Camera {
    mat4 uView;
    mat4 uProj;
    mat4 uViewProj;
    vec3 uEyeWorld;
    float _pad0;
};

layout(std140, binding = 2) uniform Material {
    vec3 uKd;    float _p0;
    vec3 uKs;    float uShininess;
    vec3 uKa;    float uOpacity;
};

#define NUM_DIR_LIGHTS 3
layout(std140, binding = 3) uniform Lights {
    vec4 uLightDir [NUM_DIR_LIGHTS]; // xyz dir, w intensity
    vec4 uLightCol [NUM_DIR_LIGHTS]; // rgb
    vec3 uAmbient; float _pad1;
};

layout(location = 0) out vec4 outColor;

vec3 shade(vec3 N, vec3 V, vec3 L, vec3 lc, float lI, vec3 baseKd) {
    float ndotl = max(dot(N, L), 0.0);
    vec3 H      = normalize(V + L);
    float ndoth = max(dot(N, H), 0.0);
    vec3 diff   = baseKd * ndotl;
    vec3 spec   = uKs * pow(ndoth, max(uShininess, 1.0));
    return lI * lc * (diff + spec);
}

void main() {
    vec3 N = normalize(fs_in.normalWorld);
    vec3 V = normalize(uEyeWorld - fs_in.posWorld);

    vec3 color = uAmbient * uKa;
    for (int i = 0; i < NUM_DIR_LIGHTS; ++i) {
        vec3 L = normalize(-uLightDir[i].xyz);
        color += shade(N, V, L, uLightCol[i].rgb, uLightDir[i].w, uKd);
    }
    // headlight
    color += shade(N, V, V, vec3(1.0), 0.35, uKd);

    outColor = vec4(color, uOpacity);
}
