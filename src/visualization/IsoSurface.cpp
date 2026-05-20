// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/IsoSurface.cpp
// =============================================================================
#include "visualization/IsoSurface.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace simall::visualization
{

namespace
{

inline util::Vec3d lerp(
    const util::Vec3d& a, const util::Vec3d& b, double sa, double sb, double iso)
{
    const double d = sb - sa;
    if (std::abs(d) < 1e-30)
        return a;
    const double t = std::clamp((iso - sa) / d, 0.0, 1.0);
    return {a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), a.z + t * (b.z - a.z)};
}

} // namespace

SurfaceMesh IsoSurface::extract(const VolumeMesh& vol, double iso)
{
    SurfaceMesh out;
    if (vol.pointScalars.size() != vol.points.size())
        return out;

    const std::size_t nTet = vol.tet_count();
    for (std::size_t t = 0; t < nTet; ++t) {
        std::array<std::int32_t, 4> idx{vol.tetIndex[4 * t + 0],
                                        vol.tetIndex[4 * t + 1],
                                        vol.tetIndex[4 * t + 2],
                                        vol.tetIndex[4 * t + 3]};
        std::array<double, 4> s{vol.pointScalars[idx[0]],
                                vol.pointScalars[idx[1]],
                                vol.pointScalars[idx[2]],
                                vol.pointScalars[idx[3]]};
        std::array<util::Vec3d, 4> p{
            vol.points[idx[0]], vol.points[idx[1]], vol.points[idx[2]], vol.points[idx[3]]};

        // Sort vertices so that scalars are ascending.  This collapses the
        // 16 sign-pattern cases to four with a uniform output.
        std::array<int, 4> order{0, 1, 2, 3};
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                if (s[order[j]] < s[order[i]])
                    std::swap(order[i], order[j]);
            }
        }
        const double s0 = s[order[0]], s1 = s[order[1]], s2 = s[order[2]], s3 = s[order[3]];
        const util::Vec3d& v0 = p[order[0]];
        const util::Vec3d& v1 = p[order[1]];
        const util::Vec3d& v2 = p[order[2]];
        const util::Vec3d& v3 = p[order[3]];

        const auto emit_tri =
            [&](const util::Vec3d& a, const util::Vec3d& b, const util::Vec3d& c) {
                const std::int32_t base = static_cast<std::int32_t>(out.points.size());
                out.points.push_back(a);
                out.points.push_back(b);
                out.points.push_back(c);
                out.pointScalars.push_back(iso);
                out.pointScalars.push_back(iso);
                out.pointScalars.push_back(iso);
                out.triIndex.push_back(base);
                out.triIndex.push_back(base + 1);
                out.triIndex.push_back(base + 2);
            };

        if (iso < s0 || iso > s3)
            continue; // no crossing
        if (iso == s0 && iso == s1 && iso == s2 && iso == s3)
            continue;

        if (iso <= s1) {
            // One vertex (v0) below or equal, three above → 1 triangle on
            // edges (v0,v1), (v0,v2), (v0,v3).
            if (iso <= s0)
                continue;
            emit_tri(
                lerp(v0, v1, s0, s1, iso), lerp(v0, v2, s0, s2, iso), lerp(v0, v3, s0, s3, iso));
        } else if (iso <= s2) {
            // Two vertices below (v0,v1), two above (v2,v3) → quad on edges
            // (v0,v2),(v1,v2),(v1,v3),(v0,v3) → emit as two triangles.
            const util::Vec3d e02 = lerp(v0, v2, s0, s2, iso);
            const util::Vec3d e12 = lerp(v1, v2, s1, s2, iso);
            const util::Vec3d e13 = lerp(v1, v3, s1, s3, iso);
            const util::Vec3d e03 = lerp(v0, v3, s0, s3, iso);
            emit_tri(e02, e12, e13);
            emit_tri(e02, e13, e03);
        } else {
            // Three vertices below (v0..v2), one above (v3) → 1 triangle on
            // edges (v0,v3),(v1,v3),(v2,v3).
            emit_tri(
                lerp(v0, v3, s0, s3, iso), lerp(v1, v3, s1, s3, iso), lerp(v2, v3, s2, s3, iso));
        }
    }
    return out;
}

} // namespace simall::visualization
