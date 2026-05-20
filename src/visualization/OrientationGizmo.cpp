// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/OrientationGizmo.cpp
// =============================================================================
#include "visualization/OrientationGizmo.hpp"

#include <array>
#include <cmath>

namespace simall::visualization
{

util::Vec3d OrientationGizmo::direction(Axis a) noexcept
{
    switch (a) {
    case Axis::PosX:
        return {1, 0, 0};
    case Axis::NegX:
        return {-1, 0, 0};
    case Axis::PosY:
        return {0, 1, 0};
    case Axis::NegY:
        return {0, -1, 0};
    case Axis::PosZ:
        return {0, 0, 1};
    case Axis::NegZ:
        return {0, 0, -1};
    }
    return {0, 0, 1};
}

util::Vec3d OrientationGizmo::up_for(Axis a) noexcept
{
    switch (a) {
    case Axis::PosY:
    case Axis::NegY:
        return {0, 0, 1};
    default:
        return {0, 1, 0};
    }
}

std::optional<Axis> OrientationGizmo::hit_test(const Camera& cam,
                                               const GizmoLayout& layout,
                                               double pxX,
                                               double pxY)
{
    const double cx = layout.cornerX * static_cast<double>(cam.viewportWidth);
    const double cy = layout.cornerY * static_cast<double>(cam.viewportHeight);
    const double r = layout.sizePx * 0.5;
    const double dx = pxX - cx;
    const double dy = pxY - cy;
    if (dx * dx + dy * dy > r * r)
        return std::nullopt;

    // Pick whichever of the six canonical directions projects closest to the
    // click point under the current camera orientation.
    constexpr std::array<Axis, 6> axes{
        Axis::PosX, Axis::NegX, Axis::PosY, Axis::NegY, Axis::PosZ, Axis::NegZ};
    auto safe_norm = [](util::Vec3d v, util::Vec3d fb) {
        const double n = v.norm();
        return (n > 1e-30) ? util::Vec3d{v.x / n, v.y / n, v.z / n} : fb;
    };
    util::Vec3d fwd = safe_norm(cam.focalPoint - cam.position, {0, 0, -1});
    util::Vec3d right = safe_norm(fwd.cross(cam.up), {1, 0, 0});
    util::Vec3d up = safe_norm(right.cross(fwd), {0, 1, 0});

    Axis best = Axis::PosZ;
    double bestDistSq = 1e300;
    for (Axis a : axes) {
        util::Vec3d d = direction(a);
        const double xs = right.dot(d);
        const double ys = up.dot(d);
        const double sx = cx + xs * r;
        const double sy = cy - ys * r;
        const double sdx = pxX - sx;
        const double sdy = pxY - sy;
        const double dsq = sdx * sdx + sdy * sdy;
        if (dsq < bestDistSq) {
            bestDistSq = dsq;
            best = a;
        }
    }
    return best;
}

void OrientationGizmo::snap(CameraController& controller, Axis axis)
{
    controller.snap_axis(direction(axis), up_for(axis));
}

} // namespace simall::visualization
