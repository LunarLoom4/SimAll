// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/ClippingPlane.cpp
// =============================================================================
#include "visualization/ClippingPlane.hpp"

#include <array>
#include <cmath>

namespace simall::visualization
{

namespace
{

struct Vertex
{
    util::Vec3d p;
    double s; // scalar (may be 0 if surface had no scalars)
    double d; // signed distance to plane
};

inline Vertex lerp_vertex(const Vertex& a, const Vertex& b)
{
    const double t = a.d / (a.d - b.d); // |a.d - b.d| > 0 guaranteed at use
    const double tc = std::clamp(t, 0.0, 1.0);
    return Vertex{util::Vec3d{a.p.x + tc * (b.p.x - a.p.x),
                              a.p.y + tc * (b.p.y - a.p.y),
                              a.p.z + tc * (b.p.z - a.p.z)},
                  a.s + tc * (b.s - a.s),
                  0.0};
}

} // namespace

SurfaceMesh ClippingPlane::clip(const SurfaceMesh& surf, const Plane& plane, bool keepBelow)
{
    SurfaceMesh out;
    if (surf.empty())
        return out;

    const util::Vec3d n = plane.normal;
    const double D = -n.dot(plane.point);
    const double sign = keepBelow ? -1.0 : 1.0;
    const bool hasScalars = surf.pointScalars.size() == surf.points.size();

    auto emit_tri = [&](const Vertex& a, const Vertex& b, const Vertex& c) {
        const std::int32_t base = static_cast<std::int32_t>(out.points.size());
        out.points.push_back(a.p);
        out.points.push_back(b.p);
        out.points.push_back(c.p);
        if (hasScalars) {
            out.pointScalars.push_back(a.s);
            out.pointScalars.push_back(b.s);
            out.pointScalars.push_back(c.s);
        }
        out.triIndex.push_back(base);
        out.triIndex.push_back(base + 1);
        out.triIndex.push_back(base + 2);
    };

    const std::size_t nTri = surf.triangle_count();
    for (std::size_t t = 0; t < nTri; ++t) {
        std::array<Vertex, 3> v;
        for (int k = 0; k < 3; ++k) {
            const std::int32_t i = surf.triIndex[3 * t + k];
            v[k].p = surf.points[i];
            v[k].s = hasScalars ? surf.pointScalars[i] : 0.0;
            v[k].d = sign * (n.dot(v[k].p) + D);
        }
        const int code = (v[0].d >= 0 ? 1 : 0) | (v[1].d >= 0 ? 2 : 0) | (v[2].d >= 0 ? 4 : 0);
        switch (code) {
        case 0:
            continue; // all below
        case 7:
            emit_tri(v[0], v[1], v[2]);
            break; // all above
        case 1: {  // only v0 above
            Vertex a = lerp_vertex(v[0], v[1]);
            Vertex b = lerp_vertex(v[0], v[2]);
            emit_tri(v[0], a, b);
            break;
        }
        case 2: { // only v1 above
            Vertex a = lerp_vertex(v[1], v[0]);
            Vertex b = lerp_vertex(v[1], v[2]);
            emit_tri(v[1], b, a);
            break;
        }
        case 4: { // only v2 above
            Vertex a = lerp_vertex(v[2], v[0]);
            Vertex b = lerp_vertex(v[2], v[1]);
            emit_tri(v[2], a, b);
            break;
        }
        case 3: { // v0,v1 above (v2 below)
            Vertex a = lerp_vertex(v[2], v[0]);
            Vertex b = lerp_vertex(v[2], v[1]);
            emit_tri(v[0], v[1], b);
            emit_tri(v[0], b, a);
            break;
        }
        case 5: { // v0,v2 above (v1 below)
            Vertex a = lerp_vertex(v[1], v[0]);
            Vertex b = lerp_vertex(v[1], v[2]);
            emit_tri(v[0], a, b);
            emit_tri(v[0], b, v[2]);
            break;
        }
        case 6: { // v1,v2 above (v0 below)
            Vertex a = lerp_vertex(v[0], v[1]);
            Vertex b = lerp_vertex(v[0], v[2]);
            emit_tri(a, v[1], v[2]);
            emit_tri(a, v[2], b);
            break;
        }
        default:
            continue;
        }
    }
    return out;
}

} // namespace simall::visualization
