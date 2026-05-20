// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/SlidingInterface.cpp
//
// Implementation: Sutherland-Hodgman polygon clipping in the common plane
// for each (target, source) pair after rigid transformation. The clip is
// done in 2-D coordinates obtained by projecting onto the local face frame
// (e1 = first edge, e2 = n × e1, n = unit area normal).
//
// Complexity is O(N_T × N_S) in the worst case; for production meshes the
// solver wraps this with an axis-aligned bounding-box prefilter on the
// target side (built once per pairing).
// =============================================================================
#include "solver/SlidingInterface.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace simall::solver
{

namespace
{

struct V3
{
    double x, y, z;
};
inline V3 sub(V3 a, V3 b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline V3 add(V3 a, V3 b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline V3 mul(V3 a, double s)
{
    return {a.x * s, a.y * s, a.z * s};
}
inline double dot(V3 a, V3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline V3 cross(V3 a, V3 b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(V3 a)
{
    return std::sqrt(dot(a, a));
}
inline V3 normalize(V3 a)
{
    double n = norm(a);
    return n > 1e-30 ? mul(a, 1.0 / n) : V3{0, 0, 0};
}

/// Rodrigues rotation of v around unit axis k by angle θ.
V3 rodrigues(V3 v, V3 k, double theta)
{
    const double c = std::cos(theta), s = std::sin(theta);
    const V3 kxv = cross(k, v);
    const double kdv = dot(k, v);
    return {v.x * c + kxv.x * s + k.x * kdv * (1 - c),
            v.y * c + kxv.y * s + k.y * kdv * (1 - c),
            v.z * c + kxv.z * s + k.z * kdv * (1 - c)};
}

/// Apply rigid transform of patch to point p.
V3 apply_xform(const SlidingPatch& P, V3 p)
{
    V3 q = sub(p, {P.rotOrigin.x, P.rotOrigin.y, P.rotOrigin.z});
    V3 axis = normalize({P.rotAxis.x, P.rotAxis.y, P.rotAxis.z});
    V3 r = rodrigues(q, axis, P.rotAngle);
    r = add(r, {P.rotOrigin.x, P.rotOrigin.y, P.rotOrigin.z});
    return add(r, {P.translation.x, P.translation.y, P.translation.z});
}

/// Collect face vertices (transformed) into a vector.
void collect_face(const meshing::Mesh& m,
                  meshing::FaceId f,
                  const SlidingPatch& P,
                  std::vector<V3>& out)
{
    const auto& F = m.faces();
    const auto& N = m.nodes();
    const int s = F.nodeOffsets[f], e = F.nodeOffsets[f + 1];
    out.clear();
    out.reserve(e - s);
    for (int k = s; k < e; ++k) {
        const auto v = F.nodeIndices[k];
        out.push_back(apply_xform(P, {N.x[v], N.y[v], N.z[v]}));
    }
}

/// 2-D Sutherland-Hodgman polygon clipping.
/// subject and clip are CCW polygons in the (u, v) plane. Returns intersection.
std::vector<std::pair<double, double>> clip_sh(
    const std::vector<std::pair<double, double>>& subject,
    const std::vector<std::pair<double, double>>& clip)
{
    if (subject.empty() || clip.empty())
        return {};
    auto inside =
        [](std::pair<double, double> p, std::pair<double, double> a, std::pair<double, double> b) {
            // Left of edge a→b (CCW) → inside.
            return (b.first - a.first) * (p.second - a.second)
                       - (b.second - a.second) * (p.first - a.first)
                   >= -1e-15;
        };
    auto intersect = [](std::pair<double, double> p1,
                        std::pair<double, double> p2,
                        std::pair<double, double> a,
                        std::pair<double, double> b) {
        const double x1 = p1.first, y1 = p1.second;
        const double x2 = p2.first, y2 = p2.second;
        const double x3 = a.first, y3 = a.second;
        const double x4 = b.first, y4 = b.second;
        const double den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        if (std::abs(den) < 1e-30)
            return p2;
        const double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
        return std::make_pair(x1 + t * (x2 - x1), y1 + t * (y2 - y1));
    };

    std::vector<std::pair<double, double>> out = subject;
    for (std::size_t i = 0; i < clip.size(); ++i) {
        if (out.empty())
            break;
        auto a = clip[i];
        auto b = clip[(i + 1) % clip.size()];
        std::vector<std::pair<double, double>> next;
        next.reserve(out.size() * 2);
        for (std::size_t j = 0; j < out.size(); ++j) {
            auto p1 = out[(j + out.size() - 1) % out.size()];
            auto p2 = out[j];
            const bool in1 = inside(p1, a, b);
            const bool in2 = inside(p2, a, b);
            if (in2) {
                if (!in1)
                    next.push_back(intersect(p1, p2, a, b));
                next.push_back(p2);
            } else if (in1) {
                next.push_back(intersect(p1, p2, a, b));
            }
        }
        out.swap(next);
    }
    return out;
}

double poly_area_2d(const std::vector<std::pair<double, double>>& p)
{
    if (p.size() < 3)
        return 0.0;
    double a = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const auto& p1 = p[i];
        const auto& p2 = p[(i + 1) % p.size()];
        a += p1.first * p2.second - p2.first * p1.second;
    }
    return 0.5 * std::abs(a);
}

} // namespace

SlidingPairTable SlidingInterface::build_pairs(const meshing::Mesh& m,
                                               const SlidingPatch& tgt,
                                               const SlidingPatch& src,
                                               double tol)
{
    const auto& F = m.faces();
    SlidingPairTable T;
    // Collect candidate faces per side.
    std::vector<meshing::FaceId> tFaces, sFaces;
    for (std::size_t f = 0; f < F.size(); ++f) {
        if (F.boundaryZone[f] == tgt.boundaryZone)
            tFaces.push_back(f);
        if (F.boundaryZone[f] == src.boundaryZone)
            sFaces.push_back(f);
    }
    T.targetFaces = tFaces;
    T.offsets.assign(tFaces.size() + 1, 0);

    std::vector<V3> tVerts, sVerts;
    std::vector<std::pair<double, double>> t2d, s2d;
    for (std::size_t ti = 0; ti < tFaces.size(); ++ti) {
        T.offsets[ti] = static_cast<int>(T.entries.size());
        const meshing::FaceId fT = tFaces[ti];
        collect_face(m, fT, tgt, tVerts);
        if (tVerts.size() < 3)
            continue;
        // Build local frame on target face.
        const V3 e0 = normalize(sub(tVerts[1], tVerts[0]));
        const V3 n = normalize(cross(sub(tVerts[1], tVerts[0]), sub(tVerts[2], tVerts[0])));
        const V3 e1 = normalize(cross(n, e0));
        const V3 org = tVerts[0];
        auto to2d = [&](V3 p) {
            const V3 d = sub(p, org);
            return std::make_pair(dot(d, e0), dot(d, e1));
        };
        t2d.clear();
        for (auto& v : tVerts)
            t2d.push_back(to2d(v));
        const double At = poly_area_2d(t2d);
        if (At < tol)
            continue;

        for (meshing::FaceId fS : sFaces) {
            collect_face(m, fS, src, sVerts);
            if (sVerts.size() < 3)
                continue;
            // Project source vertices onto target plane (drop component along n).
            s2d.clear();
            bool anyClose = false;
            for (auto& v : sVerts) {
                const V3 d = sub(v, org);
                const double zoff = dot(d, n);
                if (std::abs(zoff) < 0.5)
                    anyClose = true; // crude band filter
                s2d.push_back({dot(d, e0), dot(d, e1)});
            }
            (void) anyClose; // band filter is advisory; clipping is exact.
            const auto inter = clip_sh(t2d, s2d);
            const double Aint = poly_area_2d(inter);
            if (Aint > tol) {
                T.entries.push_back({fS, Aint / At});
            }
        }
    }
    T.offsets.back() = static_cast<int>(T.entries.size());
    // Normalize weights per target face (clipping may slightly exceed unity
    // due to non-planar source faces).
    for (std::size_t ti = 0; ti < tFaces.size(); ++ti) {
        const int a = T.offsets[ti], b = T.offsets[ti + 1];
        double sum = 0.0;
        for (int k = a; k < b; ++k)
            sum += T.entries[k].w;
        if (sum > 1e-30)
            for (int k = a; k < b; ++k)
                T.entries[k].w /= sum;
    }
    SIMALL_LOG_INFO("Sliding",
                    "built table: ",
                    T.targetFaces.size(),
                    " target faces, ",
                    T.entries.size(),
                    " (src, w) entries");
    return T;
}

double SlidingInterface::interpolate_scalar(const SlidingPairTable& T,
                                            std::size_t ti,
                                            const meshing::Mesh& m,
                                            const std::vector<double>& phi)
{
    if (ti + 1 >= T.offsets.size())
        return 0.0;
    const auto& F = m.faces();
    const int a = T.offsets[ti], b = T.offsets[ti + 1];
    double v = 0.0;
    for (int k = a; k < b; ++k) {
        const meshing::FaceId fS = T.entries[k].sourceFace;
        const meshing::CellId c = F.owner[fS];
        if (c < phi.size())
            v += T.entries[k].w * phi[c];
    }
    return v;
}

} // namespace simall::solver
