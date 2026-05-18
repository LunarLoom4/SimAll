#version 450 core
// SimAll Beta — front-to-back single-pass volume ray-cast.
// Inputs:
//   uVolume      : 3D scalar field
//   uTransfer    : 1D RGBA transfer function
//   uExitDepth   : depth of back-faces of the proxy cube (in eye-space)

layout(binding = 0) uniform sampler3D uVolume;
layout(binding = 1) uniform sampler1D uTransfer;

layout(std140, binding = 7) uniform Volume {
    mat4  uInvModelViewProj;
    vec3  uEyeWorld;       float uStepSize;
    vec3  uBoxMin;         float uMinScalar;
    vec3  uBoxMax;         float uMaxScalar;
    float uMaxSteps;       float uJitter;
    float uDensity;        float _pad;
};

in vec3 vRayDir;    // world-space ray direction from a fullscreen pass
in vec3 vRayOrigin;

layout(location = 0) out vec4 outColor;

bool intersectBox(vec3 ro, vec3 rd, out float t0, out float t1) {
    vec3 inv = 1.0 / rd;
    vec3 tA  = (uBoxMin - ro) * inv;
    vec3 tB  = (uBoxMax - ro) * inv;
    vec3 tmin = min(tA, tB);
    vec3 tmax = max(tA, tB);
    t0 = max(max(tmin.x, tmin.y), tmin.z);
    t1 = min(min(tmax.x, tmax.y), tmax.z);
    return t1 > max(t0, 0.0);
}

float sampleScalar(vec3 p) {
    vec3 uvw = (p - uBoxMin) / (uBoxMax - uBoxMin);
    return texture(uVolume, uvw).r;
}

float hash(vec2 p) {
    return fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453);
}

void main() {
    vec3 rd = normalize(vRayDir);
    vec3 ro = vRayOrigin;

    float t0, t1;
    if (!intersectBox(ro, rd, t0, t1)) discard;
    t0 = max(t0, 0.0);

    float jitter = uJitter * hash(gl_FragCoord.xy);
    float t  = t0 + jitter * uStepSize;
    vec4 dst = vec4(0.0);

    for (int i = 0; i < int(uMaxSteps); ++i) {
        if (t > t1 || dst.a > 0.99) break;
        vec3  p = ro + rd * t;
        float s = sampleScalar(p);
        float n = clamp((s - uMinScalar) / max(uMaxScalar - uMinScalar, 1e-9), 0.0, 1.0);
        vec4  src = texture(uTransfer, n);
        src.a *= uDensity * uStepSize;
        dst.rgb += (1.0 - dst.a) * src.a * src.rgb;
        dst.a   += (1.0 - dst.a) * src.a;
        t += uStepSize;
    }
    outColor = dst;
}
