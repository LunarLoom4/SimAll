#version 450 core
// SimAll Beta — feature-edge vertex stage. Renders thin lines along
// feature edges produced by meshing::FeatureEdgeExtractor.

layout(location = 0) in vec3 a_position;

layout(std140, binding = 0) uniform Camera {
    mat4 uView;
    mat4 uProj;
    mat4 uViewProj;
    vec3 uEyeWorld;
    float _pad0;
};

layout(std140, binding = 1) uniform Model {
    mat4 uModel;
    mat4 uNormalMat;
};

void main() {
    gl_Position = uViewProj * uModel * vec4(a_position, 1.0);
}
