// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/OrientationGizmo.hpp
// Phase  : 3 / Week 14 (overlays: XYZ axis triad)
//
// A small data + helper that produces:
//   * the six canonical axis directions (+X, -X, +Y, -Y, +Z, -Z)
//   * a hit-test against an axis cluster in screen space (returns which
//     axis was clicked so the GUI can snap the camera).
// =============================================================================
#pragma once

#include "visualization/CameraController.hpp"

#include <optional>
#include <utility>

namespace simall::visualization
{

enum class Axis : std::uint8_t
{
    PosX,
    NegX,
    PosY,
    NegY,
    PosZ,
    NegZ
};

struct GizmoLayout
{
    double cornerX = 0.85; // [0..1] normalized viewport
    double cornerY = 0.85;
    double sizePx = 80.0;
};

class OrientationGizmo
{
public:
    /// Convert an Axis into its world-space unit vector.
    static util::Vec3d direction(Axis a) noexcept;

    /// Reasonable up-vector to pair with each axis when snapping.
    static util::Vec3d up_for(Axis a) noexcept;

    /// Hit-test screen pixel against an icon cluster located in the
    /// configured corner of the viewport.  Returns nullopt when the click
    /// missed.  The hit zone is a circle of `sizePx` diameter.
    static std::optional<Axis> hit_test(const Camera& cam,
                                        const GizmoLayout& layout,
                                        double pxX,
                                        double pxY);

    /// Snap the controller to the chosen axis using up_for(axis).
    static void snap(CameraController& controller, Axis axis);
};

} // namespace simall::visualization
