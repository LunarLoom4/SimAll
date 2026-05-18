// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/VolumeRaycast.hpp
// Phase  : 3 / Week 14 (filters: CPU reference volume renderer)
//
// Render a regular 3D scalar grid into a 2D image with front-to-back
// compositing through a user-supplied transfer function.  This is a CPU
// reference path used by the GUI for thumbnails, the tests, and as a
// correctness oracle for the future GPU raycaster (W18.5).
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

#include <vector>

namespace simall::visualization {

struct RegularGrid3 {
    std::uint32_t      nx = 0, ny = 0, nz = 0;
    util::Vec3d        origin {0, 0, 0};
    util::Vec3d        spacing{1, 1, 1};
    std::vector<float> values;   // length nx*ny*nz, k-major (z slowest)

    bool valid() const noexcept {
        return nx > 0 && ny > 0 && nz > 0 &&
               static_cast<std::size_t>(nx) * ny * nz == values.size();
    }
    util::BoundingBox bbox() const noexcept {
        return {origin,
                {origin.x + spacing.x * (nx - 1),
                 origin.y + spacing.y * (ny - 1),
                 origin.z + spacing.z * (nz - 1)}};
    }
};

struct RaycastConfig {
    double sampleStep = 0.0;     // <=0 → automatic (0.5 * min spacing)
    double opacityScale = 1.0;   // multiplies TF alpha (1 = stops are absolute)
    Color4 background = {0.0f, 0.0f, 0.0f, 0.0f};
    bool   useEarlyTermination = true;
    double earlyOpacity = 0.98;
};

class VolumeRaycast {
public:
    /// Render `volume` from `camera` using `tf`.  Resolution comes from
    /// camera.viewportWidth/height.  Returns an RGBA8 image suitable for
    /// blitting straight onto the back buffer.
    static Image render(const RegularGrid3& volume,
                        const Camera& camera,
                        const TransferFunction& tf,
                        const RaycastConfig& cfg = {});
};

}  // namespace simall::visualization
