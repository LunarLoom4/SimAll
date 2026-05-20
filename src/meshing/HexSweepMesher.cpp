// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/HexSweepMesher.cpp
// =============================================================================
#include "meshing/HexSweepMesher.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::meshing
{

namespace
{
inline util::Vec3d sub(const util::Vec3d& a, const util::Vec3d& b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline util::Vec3d add(const util::Vec3d& a, const util::Vec3d& b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline util::Vec3d scale(const util::Vec3d& a, double s)
{
    return {a.x * s, a.y * s, a.z * s};
}
inline double dot(const util::Vec3d& a, const util::Vec3d& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline util::Vec3d cross(const util::Vec3d& a, const util::Vec3d& b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double mag(const util::Vec3d& a)
{
    return std::sqrt(dot(a, a));
}
inline util::Vec3d norm(util::Vec3d a)
{
    const double m = mag(a);
    return (m > 1e-30) ? scale(a, 1.0 / m) : util::Vec3d{1, 0, 0};
}

/// Linear-interpolated arc-length sample of `pts` at parameter s ∈ [0,1].
util::Vec3d sample_path(const std::vector<util::Vec3d>& pts, double s)
{
    if (pts.empty())
        return {0, 0, 0};
    if (pts.size() == 1)
        return pts[0];
    s = std::clamp(s, 0.0, 1.0);
    const double seg = s * (pts.size() - 1);
    const std::size_t i = std::min<std::size_t>(static_cast<std::size_t>(seg), pts.size() - 2);
    const double f = seg - i;
    return add(scale(pts[i], 1 - f), scale(pts[i + 1], f));
}
} // namespace

void HexSweepMesher::initialize(SweepProps props)
{
    p_ = props;
    SIMALL_LOG_INFO("Meshing",
                    "HexSweep init: nLayers=",
                    p_.nLayers,
                    " zones src=",
                    p_.sourceZone,
                    " tgt=",
                    p_.targetZone,
                    " side=",
                    p_.sideZone);
}

std::size_t HexSweepMesher::sweep(const SourceQuadMesh& source,
                                  const std::vector<util::Vec3d>& pathPts,
                                  Mesh& outMesh)
{
    if (source.nodes.empty() || source.quads.empty() || pathPts.size() < 2 || p_.nLayers < 2)
        return 0;

    auto& N = outMesh.nodes();
    auto& F = outMesh.faces();
    if (F.nodeOffsets.empty())
        F.nodeOffsets.push_back(0);
    auto& C = outMesh.cells();
    const NodeId baseId = static_cast<NodeId>(N.size());

    // Source-plane normal from first quad.
    const auto& q0 = source.quads[0];
    const util::Vec3d srcN = norm(cross(sub(source.nodes[q0[1]], source.nodes[q0[0]]),
                                        sub(source.nodes[q0[2]], source.nodes[q0[0]])));
    // Source-plane centroid (for translation reference).
    util::Vec3d srcC{0, 0, 0};
    for (const auto& v : source.nodes)
        srcC = add(srcC, v);
    srcC = scale(srcC, 1.0 / source.nodes.size());

    // Parallel-transport frame along path.
    util::Vec3d t_prev = norm(sub(pathPts[1], pathPts[0]));
    // Pick a stable up-vector not parallel to t_prev.
    util::Vec3d up = (std::abs(t_prev.z) < 0.9) ? util::Vec3d{0, 0, 1} : util::Vec3d{1, 0, 0};
    util::Vec3d u0 = norm(cross(t_prev, up));
    util::Vec3d v0 = cross(t_prev, u0);

    // Generate node layers.
    const std::size_t nPerLayer = source.nodes.size();
    for (std::size_t k = 0; k < p_.nLayers; ++k) {
        const double s = static_cast<double>(k) / (p_.nLayers - 1);
        const util::Vec3d centre = sample_path(pathPts, s);
        // Tangent (forward difference).
        const double ds = 1e-4;
        const util::Vec3d tang = norm(sub(sample_path(pathPts, std::min(s + ds, 1.0)),
                                          sample_path(pathPts, std::max(s - ds, 0.0))));
        // Parallel-transport: rotate (u0, v0) into plane orthogonal to tang.
        // Project them and re-orthonormalise.
        util::Vec3d u = u0;
        u = sub(u, scale(tang, dot(u, tang)));
        u = norm(u);
        util::Vec3d v = norm(cross(tang, u));
        u0 = u;
        v0 = v;
        t_prev = tang;

        const double sc = p_.scaleAlongPath ? (p_.startScale * (1 - s) + p_.endScale * s) : 1.0;
        for (const auto& q : source.nodes) {
            // Source-local (u,v) coords with respect to srcC and srcN basis.
            const util::Vec3d r = sub(q, srcC);
            // Build srcU,srcV orthogonal to srcN.
            util::Vec3d srcU = (std::abs(srcN.z) < 0.9) ? norm(cross(srcN, util::Vec3d{0, 0, 1}))
                                                        : norm(cross(srcN, util::Vec3d{1, 0, 0}));
            const util::Vec3d srcV = cross(srcN, srcU);
            const double uc = dot(r, srcU) * sc;
            const double vc = dot(r, srcV) * sc;
            const util::Vec3d pos{centre.x + uc * u.x + vc * v.x,
                                  centre.y + uc * u.y + vc * v.y,
                                  centre.z + uc * u.z + vc * v.z};
            N.x.push_back(pos.x);
            N.y.push_back(pos.y);
            N.z.push_back(pos.z);
        }
    }

    // Build hex cells & faces.
    auto cellOff_start =
        C.faceOffsets.empty() ? 0 : static_cast<std::int32_t>(C.faceOffsets.back());
    if (C.faceOffsets.empty())
        C.faceOffsets.push_back(0);
    std::size_t hexCount = 0;

    auto addFace =
        [&](NodeId a, NodeId b, NodeId c, NodeId d, CellId owner, CellId neigh, ZoneId zone) {
            F.owner.push_back(owner);
            F.neighbor.push_back(neigh);
            F.areaX.push_back(0);
            F.areaY.push_back(0);
            F.areaZ.push_back(0);
            F.centroidX.push_back(0);
            F.centroidY.push_back(0);
            F.centroidZ.push_back(0);
            F.boundaryZone.push_back(zone);
            F.nodeIndices.push_back(baseId + a);
            F.nodeIndices.push_back(baseId + b);
            F.nodeIndices.push_back(baseId + c);
            F.nodeIndices.push_back(baseId + d);
            F.nodeOffsets.push_back(static_cast<std::int32_t>(F.nodeIndices.size()));
        };

    for (std::size_t k = 0; k + 1 < p_.nLayers; ++k) {
        for (const auto& qq : source.quads) {
            const CellId cellId = static_cast<CellId>(C.size());
            // 8 node indices: 4 bottom + 4 top.
            const NodeId b0 = static_cast<NodeId>(k * nPerLayer + qq[0]);
            const NodeId b1 = static_cast<NodeId>(k * nPerLayer + qq[1]);
            const NodeId b2 = static_cast<NodeId>(k * nPerLayer + qq[2]);
            const NodeId b3 = static_cast<NodeId>(k * nPerLayer + qq[3]);
            const NodeId t0 = static_cast<NodeId>((k + 1) * nPerLayer + qq[0]);
            const NodeId t1 = static_cast<NodeId>((k + 1) * nPerLayer + qq[1]);
            const NodeId t2 = static_cast<NodeId>((k + 1) * nPerLayer + qq[2]);
            const NodeId t3 = static_cast<NodeId>((k + 1) * nPerLayer + qq[3]);
            // 6 hex faces (CCW outward).
            const ZoneId zSrc = (k == 0) ? p_.sourceZone : 0;
            const ZoneId zTgt = (k + 1 == p_.nLayers - 1) ? p_.targetZone : 0;
            const ZoneId zSide = p_.sideZone;
            addFace(b0, b3, b2, b1, cellId, kBoundaryCell, zSrc);  // bottom
            addFace(t0, t1, t2, t3, cellId, kBoundaryCell, zTgt);  // top
            addFace(b0, b1, t1, t0, cellId, kBoundaryCell, zSide); // side
            addFace(b1, b2, t2, t1, cellId, kBoundaryCell, zSide);
            addFace(b2, b3, t3, t2, cellId, kBoundaryCell, zSide);
            addFace(b3, b0, t0, t3, cellId, kBoundaryCell, zSide);
            // Cell volumes/centroids left for compute_geometry().
            C.volume.push_back(0.0);
            C.centroidX.push_back(0.0);
            C.centroidY.push_back(0.0);
            C.centroidZ.push_back(0.0);
            cellOff_start += 6;
            C.faceOffsets.push_back(cellOff_start);
            for (int ff = 0; ff < 6; ++ff)
                C.faceIndices.push_back(static_cast<FaceId>(F.size() - 6 + ff));
            ++hexCount;
        }
    }
    SIMALL_LOG_INFO(
        "Meshing", "HexSweep generated ", hexCount, " hexes across ", p_.nLayers, " layers");
    return hexCount;
}

} // namespace simall::meshing
