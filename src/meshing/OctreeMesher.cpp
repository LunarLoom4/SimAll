// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/OctreeMesher.cpp
// =============================================================================
#include "meshing/OctreeMesher.hpp"
#include "meshing/Connectivity.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

namespace simall::meshing {

namespace {

struct AABB {
    double xmin, ymin, zmin, xmax, ymax, zmax;
    double cx() const { return 0.5 * (xmin + xmax); }
    double cy() const { return 0.5 * (ymin + ymax); }
    double cz() const { return 0.5 * (zmin + zmax); }
    double dx() const { return xmax - xmin; }
    double dy() const { return ymax - ymin; }
    double dz() const { return zmax - zmin; }
};

// Triangle vs AABB SAT test (Akenine-Möller 2001).
bool tri_aabb_intersect(const util::Vec3d& v0, const util::Vec3d& v1,
                        const util::Vec3d& v2, const AABB& b) {
    const double cx = b.cx(), cy = b.cy(), cz = b.cz();
    const double hx = 0.5 * b.dx(), hy = 0.5 * b.dy(), hz = 0.5 * b.dz();
    const double V0[3] = { v0.x - cx, v0.y - cy, v0.z - cz };
    const double V1[3] = { v1.x - cx, v1.y - cy, v1.z - cz };
    const double V2[3] = { v2.x - cx, v2.y - cy, v2.z - cz };
    // AABB axis tests
    for (int a = 0; a < 3; ++a) {
        const double mn = std::min({V0[a], V1[a], V2[a]});
        const double mx = std::max({V0[a], V1[a], V2[a]});
        const double h  = (a==0)?hx:(a==1)?hy:hz;
        if (mn > h || mx < -h) return false;
    }
    // Triangle normal plane test
    const double e0[3] = { V1[0]-V0[0], V1[1]-V0[1], V1[2]-V0[2] };
    const double e1[3] = { V2[0]-V1[0], V2[1]-V1[1], V2[2]-V1[2] };
    const double n[3]  = { e0[1]*e1[2]-e0[2]*e1[1],
                           e0[2]*e1[0]-e0[0]*e1[2],
                           e0[0]*e1[1]-e0[1]*e1[0] };
    const double r = hx*std::abs(n[0]) + hy*std::abs(n[1]) + hz*std::abs(n[2]);
    const double s = n[0]*V0[0] + n[1]*V0[1] + n[2]*V0[2];
    if (std::abs(s) > r) return false;
    return true;
}

// Möller-Trumbore intersection: ray from p along +x with the triangle (v0,v1,v2).
// Returns whether t > 0 and the hit lies on the triangle.
bool ray_x_tri(const util::Vec3d& p,
               const util::Vec3d& v0, const util::Vec3d& v1, const util::Vec3d& v2,
               double& t) {
    const double d[3] = { 1.0, 0.0, 0.0 };
    const double e1[3] = { v1.x-v0.x, v1.y-v0.y, v1.z-v0.z };
    const double e2[3] = { v2.x-v0.x, v2.y-v0.y, v2.z-v0.z };
    const double pvec[3] = { d[1]*e2[2]-d[2]*e2[1],
                             d[2]*e2[0]-d[0]*e2[2],
                             d[0]*e2[1]-d[1]*e2[0] };
    const double det = e1[0]*pvec[0] + e1[1]*pvec[1] + e1[2]*pvec[2];
    if (std::abs(det) < 1e-20) return false;
    const double inv = 1.0 / det;
    const double tvec[3] = { p.x-v0.x, p.y-v0.y, p.z-v0.z };
    const double u = (tvec[0]*pvec[0]+tvec[1]*pvec[1]+tvec[2]*pvec[2]) * inv;
    if (u < 0 || u > 1) return false;
    const double qvec[3] = { tvec[1]*e1[2]-tvec[2]*e1[1],
                             tvec[2]*e1[0]-tvec[0]*e1[2],
                             tvec[0]*e1[1]-tvec[1]*e1[0] };
    const double v = (d[0]*qvec[0]+d[1]*qvec[1]+d[2]*qvec[2]) * inv;
    if (v < 0 || u + v > 1) return false;
    t = (e2[0]*qvec[0]+e2[1]*qvec[1]+e2[2]*qvec[2]) * inv;
    return t > 0.0;
}

bool point_inside(const util::Vec3d& p, const StlSurface& s) {
    int crossings = 0;
    for (const auto& tri : s.triangles) {
        double t;
        if (ray_x_tri(p, s.vertices[tri[0]], s.vertices[tri[1]], s.vertices[tri[2]], t))
            ++crossings;
    }
    return (crossings & 1) == 1;
}

struct Leaf {
    AABB box;
    bool boundary = false;
};

void subdivide(const AABB& b, std::array<AABB, 8>& out) {
    const double mx = b.cx(), my = b.cy(), mz = b.cz();
    int k = 0;
    for (int iz = 0; iz < 2; ++iz)
    for (int iy = 0; iy < 2; ++iy)
    for (int ix = 0; ix < 2; ++ix) {
        out[k++] = {
            ix ? mx : b.xmin, iy ? my : b.ymin, iz ? mz : b.zmin,
            ix ? b.xmax : mx, iy ? b.ymax : my, iz ? b.zmax : mz
        };
    }
}

void refine(const AABB& box, int depth, int maxDepth, int minDepth,
            const StlSurface& s, std::vector<Leaf>& leaves) {
    bool hitsSurface = false;
    for (const auto& t : s.triangles) {
        if (tri_aabb_intersect(s.vertices[t[0]], s.vertices[t[1]], s.vertices[t[2]], box)) {
            hitsSurface = true; break;
        }
    }
    if (depth >= maxDepth || (!hitsSurface && depth >= minDepth)) {
        leaves.push_back({box, hitsSurface});
        return;
    }
    std::array<AABB, 8> kids;
    subdivide(box, kids);
    for (auto& k : kids) refine(k, depth + 1, maxDepth, minDepth, s, leaves);
}

}  // namespace

void OctreeMesher::mesh(const StlSurface& surface,
                        OctreeMeshOptions opt, Mesh& out) {
    // Surface bbox
    AABB bb{ std::numeric_limits<double>::max(),
             std::numeric_limits<double>::max(),
             std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max() };
    for (const auto& v : surface.vertices) {
        bb.xmin = std::min(bb.xmin, v.x); bb.xmax = std::max(bb.xmax, v.x);
        bb.ymin = std::min(bb.ymin, v.y); bb.ymax = std::max(bb.ymax, v.y);
        bb.zmin = std::min(bb.zmin, v.z); bb.zmax = std::max(bb.zmax, v.z);
    }
    // Pad bbox slightly
    const double pad = 0.01 * std::max({bb.dx(), bb.dy(), bb.dz()});
    bb.xmin-=pad; bb.ymin-=pad; bb.zmin-=pad;
    bb.xmax+=pad; bb.ymax+=pad; bb.zmax+=pad;

    std::vector<Leaf> leaves;
    refine(bb, 0, opt.maxDepth, opt.minDepthGlobal, surface, leaves);

    // Classify non-boundary leaves and emit hex cells for inside + boundary.
    NodeStorage ns;
    std::vector<CellDescriptor> cells;
    cells.reserve(leaves.size());

    for (const Leaf& L : leaves) {
        if (!L.boundary) {
            util::Vec3d c{ L.box.cx(), L.box.cy(), L.box.cz() };
            if (!point_inside(c, surface)) continue;
        }
        const std::size_t n0 = ns.size();
        // 8 corner nodes (i,j,k) — order matches the hex face convention used
        // by ConnectivityBuilder elsewhere in the code base.
        const double xs[2] = { L.box.xmin, L.box.xmax };
        const double ys[2] = { L.box.ymin, L.box.ymax };
        const double zs[2] = { L.box.zmin, L.box.zmax };
        for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i) {
            ns.x.push_back(xs[i]);
            ns.y.push_back(ys[j]);
            ns.z.push_back(zs[k]);
        }
        auto P = [&](int i, int j, int k) -> NodeId {
            return NodeId(n0 + k*4 + j*2 + i);
        };
        CellDescriptor cd;
        cd.faces = {
            {P(0,0,0), P(0,1,0), P(0,1,1), P(0,0,1)},  // x-min
            {P(1,0,0), P(1,0,1), P(1,1,1), P(1,1,0)},  // x-max
            {P(0,0,0), P(0,0,1), P(1,0,1), P(1,0,0)},  // y-min
            {P(0,1,0), P(1,1,0), P(1,1,1), P(0,1,1)},  // y-max
            {P(0,0,0), P(1,0,0), P(1,1,0), P(0,1,0)},  // z-min
            {P(0,0,1), P(0,1,1), P(1,1,1), P(1,0,1)}   // z-max
        };
        cells.push_back(std::move(cd));
    }
    ConnectivityBuilder::build(out, ns, cells);
    SIMALL_LOG_INFO("Octree", leaves.size(), " leaves → ", cells.size(),
        " hex cells (depth ", opt.maxDepth, ")");
}

}  // namespace simall::meshing
