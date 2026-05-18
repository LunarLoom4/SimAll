#version 450 core
// SimAll Beta — scalar-to-colour LUT fragment shader.
// Supported colour maps: 0 = jet, 1 = viridis, 2 = plasma, 3 = grayscale.

in VS_OUT {
    vec3 posWorld;
    vec3 normalWorld;
    vec2 uv;
    float scalar;
} fs_in;

layout(std140, binding = 4) uniform Lut {
    float uMin;
    float uMax;
    int   uMap;     // 0=jet, 1=viridis, 2=plasma, 3=gray
    int   uLogScale;
    float uOpacity;
    float _p0; float _p1; float _p2;
};

layout(location = 0) out vec4 outColor;

vec3 jet(float t) {
    float r = clamp(1.5 - abs(4.0 * t - 3.0), 0.0, 1.0);
    float g = clamp(1.5 - abs(4.0 * t - 2.0), 0.0, 1.0);
    float b = clamp(1.5 - abs(4.0 * t - 1.0), 0.0, 1.0);
    return vec3(r, g, b);
}

// Compact polynomial approximations of MPL viridis / plasma.
vec3 viridis(float t) {
    t = clamp(t, 0.0, 1.0);
    return vec3(
        0.267 + t*(0.105 + t*(-0.331 + t*(2.800 + t*(-3.236 + t* 1.385)))),
        0.005 + t*(1.405 + t*(-0.298 + t*(-0.530 + t*( 0.610 + t*(-0.180))))),
        0.330 + t*(0.962 + t*(-0.683 + t*(-1.605 + t*( 3.034 + t*(-1.040)))))
    );
}

vec3 plasma(float t) {
    t = clamp(t, 0.0, 1.0);
    return vec3(
        0.050 + t*(2.158 + t*(-0.830 + t*(-1.230 + t*( 1.480 + t*(-0.625))))),
        0.030 + t*(-0.221 + t*( 0.730 + t*( 0.040 + t*(-0.470 + t*( 0.150))))),
        0.530 + t*(0.640 + t*(-1.265 + t*( 0.250 + t*( 0.490 + t*(-0.265))))));
}

void main() {
    float v = fs_in.scalar;
    if (uLogScale != 0) v = log(max(v, 1e-12));
    float t = clamp((v - uMin) / max(uMax - uMin, 1e-12), 0.0, 1.0);

    vec3 col;
    if      (uMap == 0) col = jet(t);
    else if (uMap == 1) col = viridis(t);
    else if (uMap == 2) col = plasma(t);
    else                col = vec3(t);

    outColor = vec4(col, uOpacity);
}
