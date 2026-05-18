// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MeshSmoother.cpp
// =============================================================================
#include "meshing/MeshSmoother.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <unordered_map>

namespace simall::meshing {

void MeshSmoother::build_node_neighbours(const Mesh& m,
                                         std::vector<std::int32_t>& offsets,
                                         std::vector<NodeId>&       indices) {
    const auto& F = m.faces();
    const std::size_t nN = m.nodes().size();
    std::vector<std::unordered_set<NodeId>> adj(nN);
    for (std::size_t f = 0; f < F.size(); ++f) {
        const int s = F.nodeOffsets[f], e = F.nodeOffsets[f+1];
        for (int k = s; k < e; ++k) {
            const NodeId v1 = F.nodeIndices[k];
            const NodeId v2 = F.nodeIndices[s + ((k - s + 1) % (e - s))];
            adj[v1].insert(v2);
            adj[v2].insert(v1);
        }
    }
    offsets.assign(nN + 1, 0);
    for (std::size_t i = 0; i < nN; ++i)
        offsets[i+1] = offsets[i] + static_cast<int>(adj[i].size());
    indices.assign(offsets.back(), 0);
    for (std::size_t i = 0; i < nN; ++i) {
        int idx = offsets[i];
        for (NodeId j : adj[i]) indices[idx++] = j;
    }
}

bool MeshSmoother::smooth(Mesh& m, const SmoothParams& p) {
    auto& N = m.nodes();
    const std::size_t nN = N.size();
    if (nN == 0) return true;

    // Collect pinned (boundary) nodes: a node is pinned if any of its
    // incident faces sits on a boundary zone in pinnedBoundaryZones, OR if
    // pinnedBoundaryZones is empty (default: all boundary nodes pinned).
    std::vector<std::uint8_t> pinned(nN, 0);
    const auto& F = m.faces();
    const bool pinAllBoundary = p.pinnedBoundaryZones.empty();
    std::unordered_set<ZoneId> pinSet(p.pinnedBoundaryZones.begin(),
                                       p.pinnedBoundaryZones.end());
    for (std::size_t f = 0; f < F.size(); ++f) {
        const ZoneId z = F.boundaryZone[f];
        if (z == 0) continue;
        const bool pin = pinAllBoundary || pinSet.count(z);
        if (!pin) continue;
        const int s = F.nodeOffsets[f], e = F.nodeOffsets[f+1];
        for (int k = s; k < e; ++k) pinned[F.nodeIndices[k]] = 1;
    }

    std::vector<std::int32_t> offs; std::vector<NodeId> idx;
    build_node_neighbours(m, offs, idx);

    util::aligned_vector<double> xn(nN), yn(nN), zn(nN);

    for (int iter = 0; iter < p.iterations; ++iter) {
        for (std::size_t i = 0; i < nN; ++i) {
            if (pinned[i]) { xn[i] = N.x[i]; yn[i] = N.y[i]; zn[i] = N.z[i]; continue; }
            double sx = 0, sy = 0, sz = 0, sw = 0;
            const int a = offs[i], b = offs[i+1];
            for (int k = a; k < b; ++k) {
                const NodeId j = idx[k];
                double w = 1.0;
                if (p.mode == SmoothMode::LengthWeighted
                 || p.mode == SmoothMode::Optimization) {
                    const double dx = N.x[j] - N.x[i];
                    const double dy = N.y[j] - N.y[i];
                    const double dz = N.z[j] - N.z[i];
                    w = 1.0 / (std::sqrt(dx*dx + dy*dy + dz*dz) + 1e-30);
                }
                sx += w * N.x[j]; sy += w * N.y[j]; sz += w * N.z[j];
                sw += w;
            }
            if (sw < 1e-30) { xn[i] = N.x[i]; yn[i] = N.y[i]; zn[i] = N.z[i]; continue; }
            const double xT = sx / sw, yT = sy / sw, zT = sz / sw;
            const double r = p.relax;
            xn[i] = (1 - r) * N.x[i] + r * xT;
            yn[i] = (1 - r) * N.y[i] + r * yT;
            zn[i] = (1 - r) * N.z[i] + r * zT;
        }
        N.x = xn; N.y = yn; N.z = zn;
        if (p.mode == SmoothMode::Optimization) {
            // Refresh geometry mid-pass so the next iteration sees the new
            // face areas (used implicitly through neighbour distances).
            m.compute_geometry();
        }
    }
    m.compute_geometry();

    // Validity check: no zero/negative volumes left behind.
    const auto& C = m.cells();
    for (std::size_t c = 0; c < C.size(); ++c) {
        if (C.volume[c] <= 0.0) {
            SIMALL_LOG_WARN("Smoother", "cell ", c, " volume=", C.volume[c]);
            return false;
        }
    }
    return true;
}

}  // namespace simall::meshing
