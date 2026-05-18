# GLSL Shader Set

All shaders target **GLSL 4.5 core** (matches VTK 9.2+ on the OpenGL
backend) and are loaded by `visualization::ShaderManager`.

## Pipeline catalogue

| Shader | Stage | Purpose |
|---|---|---|
| `phong.vert` / `phong.frag`   | vert + frag | default Blinn-Phong surface |
| `wireframe.geom`              | geom        | barycentric pass-through |
| `edge.vert` / `edge.frag`     | vert + frag | feature-edge lines |
| `scalar_lut.frag`             | frag        | scalar → colour with LUT |
| `transparency.frag`           | frag        | dual depth-peel layer |
| `volume_raycast.vert/frag`    | vert + frag | front-to-back volume ray-cast |
| `lic.frag`                    | frag        | 2-D Line-Integral Convolution |
| `shadow_depth.vert/frag`      | vert + frag | shadow map pre-pass |

## Uniform block bindings

| Binding | Block       | Owner |
|---:|---|---|
| 0 | Camera   | `vis::Viewport` |
| 1 | Model    | per-actor |
| 2 | Material | `vis::MaterialUniform` |
| 3 | Lights   | `vis::LightingUniform` |
| 4 | Lut      | `vis::LutUniform` |
| 5 | Edge     | `vis::EdgeUniform` |
| 6 | Peel     | `vis::DepthPeelPass` |
| 7 | Volume   | `vis::VolumeRenderer` |
| 8 | Lic      | `vis::LicFilter` |
| 9 | Shadow   | `vis::ShadowPass` |
