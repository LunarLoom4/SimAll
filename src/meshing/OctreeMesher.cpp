// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/OctreeMesher.cpp
// =============================================================================
#include "meshing/OctreeMesher.hpp"

#include "core/Logger.hpp"
#include "meshing/Connectivity.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

namespace simall::meshing
{

namespace
{

struct AABB
{
    double xmin, ymin, zmin, xmax, ymax, zmax;
    double cx() const { return 0.5 * (xmin + xmax); }
    double cy() const { return 0.5 * (ymin + ymax); }
    double cz() const { return 0.5 * (zmin + zmax); }
    double dx() const { return xmax - xmin; }
    double dy() const { return ymax - ymin; }
    double dz() const { return zmax - zmin; }
};

// Triangle vs AABB SAT test (Akenine-Möller 2001).
bool tri_aabb_intersect(const util::Vec3d& v0,
                        const util::Vec3d& v1,
                        const util::Vec3d& v2,
                        const AABB& b)
{
    const double cx = b.cx(), cy = b.cy(), cz = b.cz();
    const double hx = 0.5 * b.dx(), hy = 0.5 * b.dy(), hz = 0.5 * b.dz();
    const double V0[3] = {v0.x - cx, v0.y - cy, v0.z - cz};
    const double V1[3] = {v1.x - cx, v1.y - cy, v1.z - cz};
    const double V2[3] = {v2.x - cx, v2.y - cy, v2.z - cz};
    // AABB axis tests
    for (int a = 0; a < 3; ++a) {
        const double mn = std::min({V0[a], V1[a], V2[a]});
        const double mx = std::max({V0[a], V1[a], V2[a]});
        const double h = (a == 0) ? hx : (a == 1) ? hy : hz;
        if (mn > h || mx < -h)
            return false;
    }
    // Triangle normal plane test
    const double e0[3] = {V1[0] - V0[0], V1[1] - V0[1], V1[2] - V0[2]};
    const double e1[3] = {V2[0] - V1[0], V2[1] - V1[1], V2[2] - V1[2]};
    const double n[3] = {e0[1] * e1[2] - e0[2] * e1[1],
                         e0[2] * e1[0] - e0[0] * e1[2],
                         e0[0] * e1[1] - e0[1] * e1[0]};
    const double r = hx * std::abs(n[0]) + hy * std::abs(n[1]) + hz * std::abs(n[2]);
    const double s = n[0] * V0[0] + n[1] * V0[1] + n[2] * V0[2];
    if (std::abs(s) > r)
        return false;
    return true;
}

// Möller-Trumbore intersection: ray from p along +x with the triangle (v0,v1,v2).
// Returns whether t > 0 and the hit lies on the triangle.
bool ray_x_tri(const util::Vec3d& p,
               const util::Vec3d& v0,
               const util::Vec3d& v1,
               const util::Vec3d& v2,
               double& t)
{
    const double d[3] = {1.0, 0.0, 0.0};
    const double e1[3] = {v1.x - v0.x, v1.y - v0.y, v1.z - v0.z};
    const double e2[3] = {v2.x - v0.x, v2.y - v0.y, v2.z - v0.z};
    const double pvec[3] = {
        d[1] * e2[2] - d[2] * e2[1], d[2] * e2[0] - d[0] * e2[2], d[0] * e2[1] - d[1] * e2[0]};
    const double det = e1[0] * pvec[0] + e1[1] * pvec[1] + e1[2] * pvec[2];
    if (std::abs(det) < 1e-20)
        return false;
    const double inv = 1.0 / det;
    const double tvec[3] = {p.x - v0.x, p.y - v0.y, p.z - v0.z};
    const double u = (tvec[0] * pvec[0] + tvec[1] * pvec[1] + tvec[2] * pvec[2]) * inv;
    if (u < 0 || u > 1)
        return false;
    const double qvec[3] = {tvec[1] * e1[2] - tvec[2] * e1[1],
                            tvec[2] * e1[0] - tvec[0] * e1[2],
                            tvec[0] * e1[1] - tvec[1] * e1[0]};
    const double v = (d[0] * qvec[0] + d[1] * qvec[1] + d[2] * qvec[2]) * inv;
    if (v < 0 || u + v > 1)
        return false;
    t = (e2[0] * qvec[0] + e2[1] * qvec[1] + e2[2] * qvec[2]) * inv;
    return t > 0.0;
}

bool point_inside(const util::Vec3d& p, const StlSurface& s)
{
    int crossings = 0;
    for (const auto& tri : s.triangles) {
        double t;
        if (ray_x_tri(p, s.vertices[tri[0]], s.vertices[tri[1]], s.vertices[tri[2]], t))
            ++crossings;
    }
    return (crossings & 1) == 1;
}

struct Leaf
{
    AABB box;
    bool boundary = false;
};

void subdivide(const AABB& b, std::array<AABB, 8>& out)
{
    const double mx = b.cx(), my = b.cy(), mz = b.cz();
    int k = 0;
    for (int iz = 0; iz < 2; ++iz)
        for (int iy = 0; iy < 2; ++iy)
            for (int ix = 0; ix < 2; ++ix) {
                out[k++] = {ix ? mx : b.xmin,
                            iy ? my : b.ymin,
                            iz ? mz : b.zmin,
                            ix ? b.xmax : mx,
                            iy ? b.ymax : my,
                            iz ? b.zmax : mz};
            }
}

void refine(const AABB& box,
            int depth,
            int maxDepth,
            int minDepth,
            const StlSurface& s,
            std::vector<Leaf>& leaves)
{
    bool hitsSurface = false;
    for (const auto& t : s.triangles) {
        if (tri_aabb_intersect(s.vertices[t[0]], s.vertices[t[1]], s.vertices[t[2]], box)) {
            hitsSurface = true;
            break;
        }
    }
    if (depth >= maxDepth || (!hitsSurface && depth >= minDepth)) {
        leaves.push_back({box, hitsSurface});
        return;
    }
    std::array<AABB, 8> kids;
    subdivide(box, kids);
    for (auto& k : kids)
        refine(k, depth + 1, maxDepth, minDepth, s, leaves);
}

// =============================================================================
// Pass 13 - true cut-cell emission for BOUNDARY leaves.
// =============================================================================

// Hex corner layout (matches the VTK convention used elsewhere):
//   0=(0,0,0) 1=(1,0,0) 2=(1,1,0) 3=(0,1,0)
//   4=(0,0,1) 5=(1,0,1) 6=(1,1,1) 7=(0,1,1)
constexpr int kEdgeCorners[12][2] = {
    {0, 1},
    {2, 3},
    {4, 5},
    {6, 7}, // x-edges
    {0, 3},
    {1, 2},
    {4, 7},
    {5, 6}, // y-edges
    {0, 4},
    {1, 5},
    {2, 6},
    {3, 7} // z-edges
};

// 6 outward CCW face cycles (corner indices).
constexpr int kFaceCycles[6][4] = {
    {0, 4, 7, 3}, // x-min
    {1, 2, 6, 5}, // x-max
    {0, 1, 5, 4}, // y-min
    {3, 7, 6, 2}, // y-max
    {0, 3, 2, 1}, // z-min
    {4, 5, 6, 7}  // z-max
};

// Find the edge index (0..11) for an unordered corner pair.
inline int find_edge_index(int a, int b) noexcept
{
    if (a > b)
        std::swap(a, b);
    for (int e = 0; e < 12; ++e)
        if (kEdgeCorners[e][0] == a && kEdgeCorners[e][1] == b)
            return e;
    return -1;
}

inline util::Vec3d corner_of(const AABB& b, int c) noexcept
{
    const double x = (c & 1) ? b.xmax : b.xmin;
    const double y = (c & 2) ? b.ymax : b.ymin;
    const double z = (c & 4) ? b.zmax : b.zmin;
    return {x, y, z};
}

// Binary-search the STL crossing along the segment whose first endpoint is
// guaranteed inside the surface and second endpoint guaranteed outside.
util::Vec3d bisect_edge(util::Vec3d inP, util::Vec3d outP, const StlSurface& s, int iters)
{
    for (int k = 0; k < iters; ++k) {
        util::Vec3d m{0.5 * (inP.x + outP.x), 0.5 * (inP.y + outP.y), 0.5 * (inP.z + outP.z)};
        if (point_inside(m, s))
            inP = m;
        else
            outP = m;
    }
    return {0.5 * (inP.x + outP.x), 0.5 * (inP.y + outP.y), 0.5 * (inP.z + outP.z)};
}

// Newell-formula polygon normal (un-normalized, encodes area * direction).
util::Vec3d polygon_newell_normal(const std::vector<NodeId>& poly, const NodeStorage& ns)
{
    util::Vec3d n{0, 0, 0};
    const std::size_t k = poly.size();
    for (std::size_t i = 0; i < k; ++i) {
        const NodeId a = poly[i];
        const NodeId b = poly[(i + 1) % k];
        n.x += (ns.y[a] - ns.y[b]) * (ns.z[a] + ns.z[b]);
        n.y += (ns.z[a] - ns.z[b]) * (ns.x[a] + ns.x[b]);
        n.z += (ns.x[a] - ns.x[b]) * (ns.y[a] + ns.y[b]);
    }
    return {0.5 * n.x, 0.5 * n.y, 0.5 * n.z};
}

// Emit a fully-inside hex (helper reused below and by the legacy path).
void emit_full_hex(const AABB& box, NodeStorage& ns, std::vector<CellDescriptor>& cells)
{
    const std::size_t n0 = ns.size();
    const double xs[2] = {box.xmin, box.xmax};
    const double ys[2] = {box.ymin, box.ymax};
    const double zs[2] = {box.zmin, box.zmax};
    for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 2; ++i) {
                ns.x.push_back(xs[i]);
                ns.y.push_back(ys[j]);
                ns.z.push_back(zs[k]);
            }
    auto P = [&](int i, int j, int k) -> NodeId { return NodeId(n0 + k * 4 + j * 2 + i); };
    CellDescriptor cd;
    cd.faces = {
        {P(0, 0, 0), P(0, 1, 0), P(0, 1, 1), P(0, 0, 1)}, // x-min
        {P(1, 0, 0), P(1, 0, 1), P(1, 1, 1), P(1, 1, 0)}, // x-max
        {P(0, 0, 0), P(0, 0, 1), P(1, 0, 1), P(1, 0, 0)}, // y-min
        {P(0, 1, 0), P(1, 1, 0), P(1, 1, 1), P(0, 1, 1)}, // y-max
        {P(0, 0, 0), P(1, 0, 0), P(1, 1, 0), P(0, 1, 0)}, // z-min
        {P(0, 0, 1), P(0, 1, 1), P(1, 1, 1), P(1, 0, 1)}  // z-max
    };
    cells.push_back(std::move(cd));
}

// Process one BOUNDARY leaf, possibly emitting one cut-cell polyhedron.
void emit_cut_cell(const AABB& box,
                   const StlSurface& s,
                   int bisectIters,
                   NodeStorage& ns,
                   std::vector<CellDescriptor>& cells)
{
    // 1. Inside-classify the 8 corners.
    bool inside[8] = {false};
    util::Vec3d cpos[8];
    int nInside = 0;
    for (int c = 0; c < 8; ++c) {
        cpos[c] = corner_of(box, c);
        inside[c] = point_inside(cpos[c], s);
        if (inside[c])
            ++nInside;
    }
    if (nInside == 0)
        return;
    if (nInside == 8) {
        emit_full_hex(box, ns, cells);
        return;
    }

    // 2. Pre-register inside corners as NodeStorage entries.
    std::array<NodeId, 8> cornerId{};
    for (int c = 0; c < 8; ++c) {
        if (!inside[c])
            continue;
        cornerId[c] = static_cast<NodeId>(ns.size());
        ns.x.push_back(cpos[c].x);
        ns.y.push_back(cpos[c].y);
        ns.z.push_back(cpos[c].z);
    }

    // 3. Edge intersections via bisection.
    std::array<NodeId, 12> edgeNode;
    edgeNode.fill(std::numeric_limits<NodeId>::max());
    for (int e = 0; e < 12; ++e) {
        const int a = kEdgeCorners[e][0];
        const int b = kEdgeCorners[e][1];
        if (inside[a] == inside[b])
            continue;
        const util::Vec3d inP = inside[a] ? cpos[a] : cpos[b];
        const util::Vec3d outP = inside[a] ? cpos[b] : cpos[a];
        const util::Vec3d mid = bisect_edge(inP, outP, s, bisectIters);
        edgeNode[e] = static_cast<NodeId>(ns.size());
        ns.x.push_back(mid.x);
        ns.y.push_back(mid.y);
        ns.z.push_back(mid.z);
    }

    // 4. Per hex face: build the clipped polygon (Sutherland-Hodgman over
    //    the 4-cycle).  Record each cut-face's two intersection ids so we
    //    can assemble the cap polygon.
    std::vector<std::vector<NodeId>> facePolys;
    facePolys.reserve(6);
    struct SegRec
    {
        int face;
        NodeId a, b;
    }; // segment per cut face
    std::vector<SegRec> segments;
    segments.reserve(6);

    for (int f = 0; f < 6; ++f) {
        std::vector<NodeId> poly;
        poly.reserve(6);
        std::vector<NodeId> segEnds;
        segEnds.reserve(2);
        for (int i = 0; i < 4; ++i) {
            const int cur = kFaceCycles[f][i];
            const int nxt = kFaceCycles[f][(i + 1) & 3];
            if (inside[cur])
                poly.push_back(cornerId[cur]);
            if (inside[cur] != inside[nxt]) {
                const int e = find_edge_index(cur, nxt);
                const NodeId nid = edgeNode[e];
                poly.push_back(nid);
                segEnds.push_back(nid);
            }
        }
        if (poly.size() >= 3)
            facePolys.push_back(std::move(poly));
        if (segEnds.size() == 2) {
            segments.push_back({f, segEnds[0], segEnds[1]});
        }
    }

    // 5. Cap polygon - walk the segment graph.  Each edgeNode is shared by
    //    exactly 2 cut faces (the 2 hex faces meeting at that edge), so the
    //    graph is a single cycle.
    if (!segments.empty()) {
        std::unordered_map<NodeId, std::array<NodeId, 2>> adj;
        std::unordered_map<NodeId, int> deg;
        auto addLink = [&](NodeId u, NodeId v) {
            int d = deg[u];
            if (d < 2)
                adj[u][d] = v;
            deg[u] = d + 1;
        };
        for (const auto& s_ : segments) {
            addLink(s_.a, s_.b);
            addLink(s_.b, s_.a);
        }

        // Walk.
        std::vector<NodeId> cap;
        cap.reserve(segments.size());
        const NodeId start = segments.front().a;
        cap.push_back(start);
        NodeId prev = std::numeric_limits<NodeId>::max();
        NodeId cur = start;
        while (true) {
            const auto it = adj.find(cur);
            if (it == adj.end())
                break;
            const NodeId nxt = (it->second[0] == prev) ? it->second[1] : it->second[0];
            if (nxt == start)
                break;
            cap.push_back(nxt);
            prev = cur;
            cur = nxt;
            if (cap.size() > segments.size() + 1)
                break; // safety
        }

        // 6. Orient cap so its outward normal points AWAY from the inside region.
        if (cap.size() >= 3) {
            // Centroid of all inside-corner+edge-intersection points (good
            // inside-region anchor).
            double xC = 0, yC = 0, zC = 0;
            std::size_t nC = 0;
            for (int c = 0; c < 8; ++c)
                if (inside[c]) {
                    xC += cpos[c].x;
                    yC += cpos[c].y;
                    zC += cpos[c].z;
                    ++nC;
                }
            xC /= nC;
            yC /= nC;
            zC /= nC;

            double xCap = 0, yCap = 0, zCap = 0;
            for (NodeId v : cap) {
                xCap += ns.x[v];
                yCap += ns.y[v];
                zCap += ns.z[v];
            }
            xCap /= cap.size();
            yCap /= cap.size();
            zCap /= cap.size();

            const util::Vec3d nrm = polygon_newell_normal(cap, ns);
            const double dotOut = nrm.x * (xCap - xC) + nrm.y * (yCap - yC) + nrm.z * (zCap - zC);
            if (dotOut < 0.0)
                std::reverse(cap.begin(), cap.end());
            facePolys.push_back(std::move(cap));
        }
    }

    if (facePolys.size() < 4)
        return; // degenerate, skip
    CellDescriptor cd;
    cd.faces = std::move(facePolys);
    cells.push_back(std::move(cd));
}

} // namespace

void OctreeMesher::mesh(const StlSurface& surface, OctreeMeshOptions opt, Mesh& out)
{
    // Surface bbox
    AABB bb{std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max()};
    for (const auto& v : surface.vertices) {
        bb.xmin = std::min(bb.xmin, v.x);
        bb.xmax = std::max(bb.xmax, v.x);
        bb.ymin = std::min(bb.ymin, v.y);
        bb.ymax = std::max(bb.ymax, v.y);
        bb.zmin = std::min(bb.zmin, v.z);
        bb.zmax = std::max(bb.zmax, v.z);
    }
    // Pad bbox slightly
    const double pad = 0.01 * std::max({bb.dx(), bb.dy(), bb.dz()});
    bb.xmin -= pad;
    bb.ymin -= pad;
    bb.zmin -= pad;
    bb.xmax += pad;
    bb.ymax += pad;
    bb.zmax += pad;

    std::vector<Leaf> leaves;
    refine(bb, 0, opt.maxDepth, opt.minDepthGlobal, surface, leaves);

    // Classify non-boundary leaves and emit hex cells for inside + boundary.
    NodeStorage ns;
    std::vector<CellDescriptor> cells;
    cells.reserve(leaves.size());

    for (const Leaf& L : leaves) {
        if (!L.boundary) {
            util::Vec3d c{L.box.cx(), L.box.cy(), L.box.cz()};
            if (!point_inside(c, surface))
                continue;
            emit_full_hex(L.box, ns, cells);
            continue;
        }
        // Boundary leaf.  In legacy (stepped) mode we keep the entire hex
        // unconditionally; in cut-cell mode we clip it against the surface.
        if (opt.enableCutCells) {
            emit_cut_cell(L.box, surface, std::max(4, opt.edgeBisectIters), ns, cells);
        } else {
            emit_full_hex(L.box, ns, cells);
        }
    }
    ConnectivityBuilder::build(out, ns, cells);
    SIMALL_LOG_INFO("Octree",
                    leaves.size(),
                    " leaves → ",
                    cells.size(),
                    " cells (depth ",
                    opt.maxDepth,
                    opt.enableCutCells ? ", cut-cell" : ", stepped",
                    ")");
}

} // namespace simall::meshing
