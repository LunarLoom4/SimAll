// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/SnappyHexMesher.cpp
// Phase  : 23 Pass 10
// =============================================================================
#include "meshing/SnappyHexMesher.hpp"
#include "meshing/OctreeMesher.hpp"
#include "meshing/PrismLayerExtruder.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace simall::meshing {

namespace {

using util::Vec3d;

// -----------------------------------------------------------------------------
// Closest point on triangle (Ericson, Real-Time Collision Detection §5.1.5).
// Returns the point on triangle (a,b,c) nearest to p in Euclidean distance.
// -----------------------------------------------------------------------------
Vec3d closest_point_on_triangle(const Vec3d& p,
                                const Vec3d& a,
                                const Vec3d& b,
                                const Vec3d& c) noexcept {
    const Vec3d ab = b - a;
    const Vec3d ac = c - a;
    const Vec3d ap = p - a;

    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0) return a;                       // vertex region A

    const Vec3d bp = p - b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3) return b;                        // vertex region B

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {                  // edge region AB
        const double v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const Vec3d cp = p - c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6) return c;                        // vertex region C

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {                  // edge region AC
        const double w = d2 / (d2 - d6);
        return a + ac * w;
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {    // edge region BC
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }

    // interior of the face
    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    return a + ab * v + ac * w;
}

// -----------------------------------------------------------------------------
// Brute-force nearest-point-on-STL.  Adequate for boundary-node snapping
// in Pass 10 (snap count is small compared to volume cell count); a BVH
// acceleration structure is queued for a later perf-tuning pass.
// -----------------------------------------------------------------------------
Vec3d nearest_on_surface(const Vec3d& p,
                         const StlSurface& s,
                         double& dist2Out) noexcept {
    Vec3d  best{};
    double bestD2 = std::numeric_limits<double>::max();
    for (const auto& tri : s.triangles) {
        const Vec3d q = closest_point_on_triangle(p,
                                                   s.vertices[tri[0]],
                                                   s.vertices[tri[1]],
                                                   s.vertices[tri[2]]);
        const Vec3d d = q - p;
        const double d2 = d.dot(d);
        if (d2 < bestD2) { bestD2 = d2; best = q; }
    }
    dist2Out = bestD2;
    return best;
}

// -----------------------------------------------------------------------------
// Collect the set of node ids that are touched by at least one boundary face,
// along with a list of boundary face indices.
// -----------------------------------------------------------------------------
void collect_boundary(const Mesh& m,
                      std::unordered_set<NodeId>& boundaryNodes,
                      std::vector<std::size_t>&   boundaryFaceIdx) {
    const auto& f = m.faces();
    const std::size_t nF = f.size();
    for (std::size_t fi = 0; fi < nF; ++fi) {
        if (f.neighbor[fi] != kBoundaryCell) continue;
        boundaryFaceIdx.push_back(fi);
        const std::int32_t b = f.nodeOffsets[fi];
        const std::int32_t e = f.nodeOffsets[fi + 1];
        for (std::int32_t k = b; k < e; ++k) {
            boundaryNodes.insert(f.nodeIndices[k]);
        }
    }
}

// -----------------------------------------------------------------------------
// Estimate the finest characteristic cell edge length so the snap-distance
// budget scales with the input geometry.
// -----------------------------------------------------------------------------
double estimate_finest_edge(const StlSurface& surface, int maxDepth) noexcept {
    util::BoundingBox bb;
    for (const auto& v : surface.vertices) bb.expand(v);
    const Vec3d ext = bb.extent();
    const double diag = std::max({ext.x, ext.y, ext.z});
    return diag / static_cast<double>(1 << std::max(1, maxDepth));
}

}  // anonymous namespace

// =============================================================================
SnappyHexStats SnappyHexMesher::mesh_background_only(
        const StlSurface&  surface,
        SnappyHexOptions   opt,
        Mesh&              backgroundOut) {
    SnappyHexStats stats{};

    OctreeMeshOptions oo{};
    oo.maxDepth       = opt.maxDepth;
    oo.minDepthGlobal = opt.minDepth;
    OctreeMesher{}.mesh(surface, oo, backgroundOut);

    stats.backgroundCells = backgroundOut.cells().size();
    stats.backgroundNodes = backgroundOut.nodes().size();
    return stats;
}

// =============================================================================
SnappyHexStats SnappyHexMesher::mesh(const StlSurface&  surface,
                                     SnappyHexOptions   opt,
                                     Mesh&              backgroundOut,
                                     Mesh&              prismLayerOut) {
    // ---- 1. Castellation --------------------------------------------------
    SnappyHexStats stats = mesh_background_only(surface, opt, backgroundOut);

    if (stats.backgroundCells == 0) {
        SIMALL_LOG_INFO("SnappyHex",
                        "castellation produced 0 cells; skipping snap/layer");
        return stats;
    }

    // ---- 2. Snapping ------------------------------------------------------
    std::unordered_set<NodeId> bNodeSet;
    std::vector<std::size_t>   bFaceIdx;
    collect_boundary(backgroundOut, bNodeSet, bFaceIdx);
    stats.boundaryNodes = bNodeSet.size();

    if (opt.enableSnapping && opt.nSnapIters > 0 && !bNodeSet.empty()) {
        const double edge       = estimate_finest_edge(surface, opt.maxDepth);
        const double budget     = edge * std::max(0.0, opt.snapMaxDistFrac);
        const double budget2    = budget * budget;

        std::vector<NodeId> bNodes(bNodeSet.begin(), bNodeSet.end());
        std::sort(bNodes.begin(), bNodes.end());  // deterministic order

        std::unordered_set<NodeId> movedAtLeastOnce;
        auto& nx = backgroundOut.nodes().x;
        auto& ny = backgroundOut.nodes().y;
        auto& nz = backgroundOut.nodes().z;

        for (int iter = 0; iter < opt.nSnapIters; ++iter) {
            const double relax = 1.0 /
                static_cast<double>(opt.nSnapIters - iter);
            for (NodeId nid : bNodes) {
                const Vec3d p{nx[nid], ny[nid], nz[nid]};
                double d2 = 0.0;
                const Vec3d q = nearest_on_surface(p, surface, d2);
                if (d2 > budget2) continue;            // beyond snap budget
                if (d2 < 1e-30)   continue;            // already on surface
                const Vec3d step = (q - p) * relax;
                nx[nid] += step.x;
                ny[nid] += step.y;
                nz[nid] += step.z;
                movedAtLeastOnce.insert(nid);
            }
        }
        stats.snapItersRun = opt.nSnapIters;
        stats.snappedNodes = movedAtLeastOnce.size();

        // Refresh geometry (areas, centroids, volumes) after node motion.
        backgroundOut.compute_geometry();
    }

    // ---- 3. Layer addition (optional) ------------------------------------
    if (opt.nLayers > 0 && !bFaceIdx.empty()) {
        // Build a compact wall-node table and a triangle list (quads split
        // into two triangles).  All other face shapes are skipped with a
        // warning - the octree background produces only quads, but a future
        // poly background could feed mixed n-gons here.
        std::unordered_map<NodeId, std::uint32_t> g2l;
        std::vector<Vec3d>                          wallNodes;
        std::vector<std::array<std::uint32_t, 3>>   wallTris;
        std::vector<Vec3d>                          normalsAccum;
        std::vector<double>                         areaAccum;

        auto local_idx = [&](NodeId nid) -> std::uint32_t {
            auto it = g2l.find(nid);
            if (it != g2l.end()) return it->second;
            const auto idx = static_cast<std::uint32_t>(wallNodes.size());
            g2l.emplace(nid, idx);
            wallNodes.push_back(Vec3d{backgroundOut.nodes().x[nid],
                                       backgroundOut.nodes().y[nid],
                                       backgroundOut.nodes().z[nid]});
            normalsAccum.emplace_back();   // zero-init
            areaAccum.push_back(0.0);
            return idx;
        };

        auto add_tri = [&](NodeId a, NodeId b, NodeId c) {
            const std::uint32_t ia = local_idx(a);
            const std::uint32_t ib = local_idx(b);
            const std::uint32_t ic = local_idx(c);
            const Vec3d va = wallNodes[ia];
            const Vec3d vb = wallNodes[ib];
            const Vec3d vc = wallNodes[ic];
            const Vec3d n  = (vb - va).cross(vc - va);   // unnormalized
            const double area = n.norm() * 0.5;
            if (area <= 0.0) return;
            wallTris.push_back({ia, ib, ic});
            normalsAccum[ia] = normalsAccum[ia] + n;
            normalsAccum[ib] = normalsAccum[ib] + n;
            normalsAccum[ic] = normalsAccum[ic] + n;
            areaAccum[ia]   += area;
            areaAccum[ib]   += area;
            areaAccum[ic]   += area;
        };

        const auto& f = backgroundOut.faces();
        for (std::size_t fi : bFaceIdx) {
            const std::int32_t b = f.nodeOffsets[fi];
            const std::int32_t e = f.nodeOffsets[fi + 1];
            const std::int32_t n = e - b;
            if (n == 3) {
                add_tri(f.nodeIndices[b],
                        f.nodeIndices[b + 1],
                        f.nodeIndices[b + 2]);
            } else if (n == 4) {
                add_tri(f.nodeIndices[b],
                        f.nodeIndices[b + 1],
                        f.nodeIndices[b + 2]);
                add_tri(f.nodeIndices[b],
                        f.nodeIndices[b + 2],
                        f.nodeIndices[b + 3]);
            } else {
                SIMALL_LOG_INFO("SnappyHex",
                                "skipping boundary face with ", n,
                                " nodes (only tri/quad supported in P10)");
            }
        }

        // Normalise the accumulated normals.
        std::vector<Vec3d> wallNormals(wallNodes.size());
        for (std::size_t k = 0; k < wallNodes.size(); ++k) {
            wallNormals[k] = normalsAccum[k].normalized();
        }

        if (!wallTris.empty()) {
            PrismLayerOptions po{};
            po.nLayers          = opt.nLayers;
            po.firstLayerHeight = opt.firstLayerHeight;
            po.growthRatio      = opt.layerGrowthRatio;
            PrismLayerExtruder{}.extrude(wallNodes, wallTris,
                                         wallNormals, po, prismLayerOut);

            stats.prismLayerCells = prismLayerOut.cells().size();
            stats.prismLayerNodes = prismLayerOut.nodes().size();
        }
    }

    SIMALL_LOG_INFO("SnappyHex",
                    "castellated=", stats.backgroundCells,
                    " bndNodes=",   stats.boundaryNodes,
                    " snapped=",    stats.snappedNodes,
                    " prismCells=", stats.prismLayerCells);
    return stats;
}

}  // namespace simall::meshing
