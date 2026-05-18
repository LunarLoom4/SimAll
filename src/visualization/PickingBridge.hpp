// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/PickingBridge.hpp
// Phase  : 3 / Week 14 (interaction: ray casting + picking)
//
// Convert screen-space picks into actor / triangle hits.  Iterates the
// registry's SurfaceMesh actors and runs the Möller–Trumbore ray–triangle
// intersection.  The bridge exists so the GUI can fire selection signals
// through core::events::SelectionChanged without taking a dependency on
// VTK's vtkCellPicker.
// =============================================================================
#pragma once

#include "visualization/ActorRegistry.hpp"
#include "visualization/VisualizationTypes.hpp"

#include <optional>

namespace simall::visualization {

struct PickResult {
    ActorId             actor     = kInvalidActorId;
    std::uint32_t       triangle  = 0;
    util::Vec3d         worldHit;
    double              distance  = 0.0;
};

struct Ray {
    util::Vec3d origin;
    util::Vec3d direction;
};

class PickingBridge {
public:
    /// Convert (pxX, pxY) from the camera viewport into a world-space ray.
    static Ray ray_from_pixel(const Camera& cam, double pxX, double pxY);

    /// Intersect the ray with every pickable Surface actor in the registry.
    /// Returns the closest hit (smallest positive t), or nullopt on miss.
    static std::optional<PickResult> pick(const ActorRegistry& reg,
                                           const Ray& ray);

    /// Convenience: combine ray_from_pixel + pick.
    static std::optional<PickResult> pick_pixel(const ActorRegistry& reg,
                                                 const Camera& cam,
                                                 double pxX, double pxY);

    /// Möller–Trumbore ray–triangle.  Returns the t parameter on hit.
    static std::optional<double> intersect_triangle(const Ray& ray,
                                                     const util::Vec3d& a,
                                                     const util::Vec3d& b,
                                                     const util::Vec3d& c);
};

}  // namespace simall::visualization
