#version 450 core
// SimAll Beta — Line-Integral Convolution for 2D vector fields
// projected into screen space (used for slice / streamline previews).

layout(binding = 0) uniform sampler2D uNoise;       // white-noise field
layout(binding = 1) uniform sampler2D uVectorField; // .rg = vx, vy normalised
in vec2 vUV;

layout(std140, binding = 8) uniform Lic {
    float uStep;
    int   uSteps;
    float uContrast;
    float _pad;
};

layout(location = 0) out vec4 outColor;

vec2 sampleVec(vec2 p) {
    vec2 v = texture(uVectorField, p).rg * 2.0 - 1.0;
    float l = length(v);
    return l > 1e-6 ? v / l : vec2(0.0);
}

float advect(vec2 origin, vec2 dir, int steps, float step) {
    float acc = 0.0;
    float w   = 0.0;
    vec2 p = origin;
    for (int i = 0; i < steps; ++i) {
        vec2 v = sampleVec(p);
        if (dot(v, v) < 1e-12) break;
        v *= dir.x;                 // forward/backward sign
        p += v * step;
        float wi = 1.0 - float(i) / float(steps);
        acc += wi * texture(uNoise, p).r;
        w   += wi;
    }
    return w > 0.0 ? acc / w : 0.0;
}

void main() {
    float f = advect(vUV, vec2( 1.0), uSteps, uStep);
    float b = advect(vUV, vec2(-1.0), uSteps, uStep);
    float lic = 0.5 * (f + b);
    lic = pow(clamp(lic, 0.0, 1.0), uContrast);
    outColor = vec4(vec3(lic), 1.0);
}
