#version 450 core
layout(location = 0) out vec4 outColor;

layout(std140, binding = 5) uniform Edge {
    vec4 uEdgeColor;       // rgba
};

void main() {
    outColor = uEdgeColor;
}
