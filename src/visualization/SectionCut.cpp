// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/SectionCut.cpp
// =============================================================================
#include "visualization/SectionCut.hpp"

#include <array>
#include <cmath>

namespace simall::visualization {

namespace {

struct Vertex {
    util::Vec3d p;
    double      s;
    double      d;
};

inline Vertex lerp_vertex(const Vertex& a, const Vertex& b) {
    const double t = a.d / (a.d - b.d);
    const double tc = std::clamp(t, 0.0, 1.0);
    return Vertex{
        util::Vec3d{a.p.x + tc*(b.p.x-a.p.x),
                    a.p.y + tc*(b.p.y-a.p.y),
                    a.p.z + tc*(b.p.z-a.p.z)},
        a.s + tc * (b.s - a.s),
        0.0};
}

}  // namespace

SurfaceMesh SectionCut::slice(const VolumeMesh& vol, const Plane& plane)
{
    SurfaceMesh out;
    if (vol.empty()) return out;

    const util::Vec3d n = plane.normal;
    const double      D = -n.dot(plane.point);
    const bool hasScalars = vol.pointScalars.size() == vol.points.size();

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

    const std::size_t nTet = vol.tet_count();
    for (std::size_t t = 0; t < nTet; ++t) {
        std::array<Vertex, 4> v;
        for (int k = 0; k < 4; ++k) {
            const std::int32_t i = vol.tetIndex[4*t + k];
            v[k].p = vol.points[i];
            v[k].s = hasScalars ? vol.pointScalars[i] : 0.0;
            v[k].d = n.dot(v[k].p) + D;
        }
        // Partition into "above" and "below" lists preserving order.
        std::array<int, 4> above{}, below{};
        int nA = 0, nB = 0;
        for (int k = 0; k < 4; ++k) {
            if (v[k].d >= 0.0) above[nA++] = k;
            else                below[nB++] = k;
        }
        if (nA == 0 || nB == 0) continue;

        // 1 above + 3 below  → triangle on 3 edges from the single above vertex.
        if (nA == 1) {
            const Vertex& A = v[above[0]];
            Vertex i0 = lerp_vertex(A, v[below[0]]);
            Vertex i1 = lerp_vertex(A, v[below[1]]);
            Vertex i2 = lerp_vertex(A, v[below[2]]);
            emit_tri(i0, i1, i2);
        } else if (nA == 3) {
            const Vertex& B = v[below[0]];
            Vertex i0 = lerp_vertex(B, v[above[0]]);
            Vertex i1 = lerp_vertex(B, v[above[1]]);
            Vertex i2 = lerp_vertex(B, v[above[2]]);
            emit_tri(i0, i1, i2);
        } else {  // nA == 2  → quad
            const Vertex& A0 = v[above[0]];
            const Vertex& A1 = v[above[1]];
            const Vertex& B0 = v[below[0]];
            const Vertex& B1 = v[below[1]];
            Vertex e00 = lerp_vertex(A0, B0);
            Vertex e01 = lerp_vertex(A0, B1);
            Vertex e10 = lerp_vertex(A1, B0);
            Vertex e11 = lerp_vertex(A1, B1);
            // Quad order around the slice: e00 → e10 → e11 → e01.
            emit_tri(e00, e10, e11);
            emit_tri(e00, e11, e01);
        }
    }
    return out;
}

}  // namespace simall::visualization
