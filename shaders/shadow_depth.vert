#version 450 core
// SimAll Beta — shadow map depth pass (vertex).
layout(location = 0) in vec3 a_position;

layout(std140, binding = 9) uniform Shadow {
    mat4 uLightViewProj;
    mat4 uModel;
};

void main() {
    gl_Position = uLightViewProj * uModel * vec4(a_position, 1.0);
}
