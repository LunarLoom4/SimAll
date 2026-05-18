#version 450 core
// SimAll Beta — single layer of dual depth-peeling for order-
// independent transparency. Inputs are bound depth textures of the
// previous front and back layers.

in VS_OUT {
    vec3 posWorld;
    vec3 normalWorld;
    vec2 uv;
    float scalar;
} fs_in;

layout(binding = 0) uniform sampler2D uPrevFrontDepth;
layout(binding = 1) uniform sampler2D uPrevBackDepth;
layout(binding = 2) uniform sampler2D uColorAccum;

layout(std140, binding = 6) uniform Peel {
    vec2 uViewport;
    float uAlpha;
    int   uIsFirstPeel;
};

layout(location = 0) out vec4 outColor;
layout(location = 1) out vec2 outDepth;

void main() {
    vec2 uv = gl_FragCoord.xy / uViewport;
    float z = gl_FragCoord.z;

    if (uIsFirstPeel == 0) {
        float front = texture(uPrevFrontDepth, uv).r;
        float back  = texture(uPrevBackDepth,  uv).r;
        // accept only fragments strictly between the previous front and back
        if (z <= front || z >= back) discard;
    }

    // Material sampling — keep this simple; replace with Phong shading
    // in a real bind by reading the material UBO.
    vec3 N = normalize(fs_in.normalWorld);
    vec3 col = vec3(0.7) * (0.4 + 0.6 * max(N.z, 0.0));
    outColor = vec4(col * uAlpha, uAlpha);
    outDepth = vec2(z, z);
}
