// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/PickingBridge.cpp
// =============================================================================
#include "visualization/PickingBridge.hpp"

#include <cmath>

namespace simall::visualization
{

namespace
{
inline util::Vec3d safe_normalize(const util::Vec3d& v, const util::Vec3d& fb)
{
    const double n = v.norm();
    return (n > 1e-30) ? util::Vec3d{v.x / n, v.y / n, v.z / n} : fb;
}
} // namespace

Ray PickingBridge::ray_from_pixel(const Camera& cam, double pxX, double pxY)
{
    util::Vec3d fwd = safe_normalize(cam.focalPoint - cam.position, {0, 0, -1});
    util::Vec3d right = safe_normalize(fwd.cross(cam.up), {1, 0, 0});
    util::Vec3d up = safe_normalize(right.cross(fwd), {0, 1, 0});

    const double sx = (2.0 * (pxX + 0.5) / cam.viewportWidth - 1.0);
    const double sy = (1.0 - 2.0 * (pxY + 0.5) / cam.viewportHeight);
    const double aspect = cam.aspect();
    const double halfH = cam.orthographic
                             ? cam.orthoHeight * 0.5
                             : std::tan(cam.fovYDegrees * 0.5 * 3.14159265358979323846 / 180.0);
    const double halfW = halfH * aspect;

    Ray ray;
    if (cam.orthographic) {
        ray.origin = cam.position + right * (sx * halfW) + up * (sy * halfH);
        ray.direction = fwd;
    } else {
        ray.origin = cam.position;
        ray.direction = safe_normalize(fwd + right * (sx * halfW) + up * (sy * halfH), {0, 0, -1});
    }
    return ray;
}

std::optional<double> PickingBridge::intersect_triangle(const Ray& r,
                                                        const util::Vec3d& a,
                                                        const util::Vec3d& b,
                                                        const util::Vec3d& c)
{
    constexpr double kEps = 1e-12;
    const util::Vec3d e1 = b - a;
    const util::Vec3d e2 = c - a;
    const util::Vec3d p = r.direction.cross(e2);
    const double det = e1.dot(p);
    if (std::abs(det) < kEps)
        return std::nullopt;
    const double invDet = 1.0 / det;
    const util::Vec3d t = r.origin - a;
    const double u = t.dot(p) * invDet;
    if (u < 0.0 || u > 1.0)
        return std::nullopt;
    const util::Vec3d q = t.cross(e1);
    const double v = r.direction.dot(q) * invDet;
    if (v < 0.0 || (u + v) > 1.0)
        return std::nullopt;
    const double tt = e2.dot(q) * invDet;
    if (tt <= kEps)
        return std::nullopt;
    return tt;
}

std::optional<PickResult> PickingBridge::pick(const ActorRegistry& reg, const Ray& ray)
{
    PickResult best;
    double bestT = 1e300;
    bool found = false;

    const auto ids = reg.list(ActorKind::Surface);
    for (ActorId id : ids) {
        const auto recOpt = reg.get(id);
        if (!recOpt)
            continue;
        const ActorRecord& rec = *recOpt;
        if (!rec.display.visible || !rec.display.picking)
            continue;
        const SurfaceMesh* surf = std::get_if<SurfaceMesh>(&rec.payload);
        if (!surf)
            continue;
        const std::size_t nTri = surf->triangle_count();
        for (std::size_t t = 0; t < nTri; ++t) {
            const auto& A = surf->points[surf->triIndex[3 * t + 0]];
            const auto& B = surf->points[surf->triIndex[3 * t + 1]];
            const auto& C = surf->points[surf->triIndex[3 * t + 2]];
            auto hit = intersect_triangle(ray, A, B, C);
            if (hit && *hit < bestT) {
                bestT = *hit;
                best.actor = id;
                best.triangle = static_cast<std::uint32_t>(t);
                best.distance = *hit;
                best.worldHit = ray.origin + ray.direction * (*hit);
                found = true;
            }
        }
    }
    if (!found)
        return std::nullopt;
    return best;
}

std::optional<PickResult> PickingBridge::pick_pixel(const ActorRegistry& reg,
                                                    const Camera& cam,
                                                    double pxX,
                                                    double pxY)
{
    return pick(reg, ray_from_pixel(cam, pxX, pxY));
}

} // namespace simall::visualization
