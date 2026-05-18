// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/ContourFilter.cpp
// =============================================================================
#include "visualization/ContourFilter.hpp"

#include <cmath>

namespace simall::visualization {

namespace {

inline util::Vec3d edge_lerp(const util::Vec3d& a, const util::Vec3d& b,
                              double sa, double sb, double iso) {
    const double denom = sb - sa;
    if (std::abs(denom) < 1e-30) return a;
    const double t = (iso - sa) / denom;
    const double tc = std::clamp(t, 0.0, 1.0);
    return util::Vec3d{
        a.x + tc * (b.x - a.x),
        a.y + tc * (b.y - a.y),
        a.z + tc * (b.z - a.z)};
}

}  // namespace

LineSet ContourFilter::extract(const SurfaceMesh& surface,
                               const std::vector<double>& isovalues)
{
    LineSet out;
    if (surface.pointScalars.size() != surface.points.size()) return out;
    out.lineOffsets.push_back(0);

    const std::size_t nTri = surface.triangle_count();

    for (double iso : isovalues) {
        for (std::size_t t = 0; t < nTri; ++t) {
            const std::int32_t a = surface.triIndex[3*t + 0];
            const std::int32_t b = surface.triIndex[3*t + 1];
            const std::int32_t c = surface.triIndex[3*t + 2];
            const double sa = surface.pointScalars[a];
            const double sb = surface.pointScalars[b];
            const double sc = surface.pointScalars[c];

            // 3-bit case code: bit i set when scalar > iso.
            const int code = (sa > iso ? 1 : 0)
                           | (sb > iso ? 2 : 0)
                           | (sc > iso ? 4 : 0);
            if (code == 0 || code == 7) continue;        // no crossing

            const util::Vec3d& A = surface.points[a];
            const util::Vec3d& B = surface.points[b];
            const util::Vec3d& C = surface.points[c];

            // Each crossing yields exactly one segment AB-AC, AB-BC, or BC-AC.
            util::Vec3d p0{}, p1{};
            switch (code) {
                case 1: case 6:   // vertex A on one side
                    p0 = edge_lerp(A, B, sa, sb, iso);
                    p1 = edge_lerp(A, C, sa, sc, iso);
                    break;
                case 2: case 5:   // vertex B on one side
                    p0 = edge_lerp(A, B, sa, sb, iso);
                    p1 = edge_lerp(B, C, sb, sc, iso);
                    break;
                case 3: case 4:   // vertex C on one side
                    p0 = edge_lerp(A, C, sa, sc, iso);
                    p1 = edge_lerp(B, C, sb, sc, iso);
                    break;
                default: continue;
            }

            const std::int32_t base =
                static_cast<std::int32_t>(out.points.size());
            out.points.push_back(p0);
            out.points.push_back(p1);
            out.scalars.push_back(iso);
            out.scalars.push_back(iso);
            out.lineOffsets.push_back(base + 2);
        }
    }
    return out;
}

LineSet ContourFilter::extract_uniform(const SurfaceMesh& surface,
                                       double sMin, double sMax,
                                       std::size_t count)
{
    std::vector<double> levels;
    if (count == 0 || sMax <= sMin) return ContourFilter::extract(surface, levels);
    levels.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const double t = (count == 1) ? 0.5
            : static_cast<double>(i) / static_cast<double>(count - 1);
        levels.push_back(sMin + t * (sMax - sMin));
    }
    return extract(surface, levels);
}

}  // namespace simall::visualization
