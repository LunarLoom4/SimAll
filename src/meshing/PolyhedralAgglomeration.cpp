// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/PolyhedralAgglomeration.cpp
// =============================================================================
#include "meshing/PolyhedralAgglomeration.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <queue>
#include <unordered_map>
#include <vector>

namespace simall::meshing {

void PolyhedralAgglomeration::initialize(AgglomerationProps props) {
    p_ = props;
    SIMALL_LOG_INFO("Meshing",
        "Agglomeration init: maxCluster=", p_.maxClusterSize,
        " quality=", p_.qualityThreshold, " seedByNode=", p_.seedByNode);
}

std::size_t PolyhedralAgglomeration::agglomerate(
        const Mesh& mesh,
        std::vector<std::int32_t>& clusterId) {
    const auto& C  = mesh.cells();
    const auto& Ff = mesh.faces();
    const std::size_t nC = C.size();
    clusterId.assign(nC, -1);

    // Build cell-to-cell neighbour list from face adjacency.
    std::vector<std::vector<CellId>> nbr(nC);
    for (std::size_t f = 0; f < Ff.size(); ++f) {
        const auto o = Ff.owner[f];
        const auto n = Ff.neighbor[f];
        if (n == kBoundaryCell) continue;
        nbr[o].push_back(n);
        nbr[n].push_back(o);
    }

    // Greedy seeded growth: pick an un-clustered cell as a seed, BFS-expand
    // into neighbours until the cluster reaches maxClusterSize.
    std::int32_t nextId = 0;
    for (std::size_t seed = 0; seed < nC; ++seed) {
        if (clusterId[seed] >= 0) continue;
        std::queue<CellId> q;
        q.push(static_cast<CellId>(seed));
        clusterId[seed] = nextId;
        std::size_t grown = 1;
        while (!q.empty() && grown < p_.maxClusterSize) {
            const auto c = q.front(); q.pop();
            for (auto nbCell : nbr[c]) {
                if (clusterId[nbCell] >= 0) continue;
                clusterId[nbCell] = nextId;
                q.push(nbCell);
                if (++grown >= p_.maxClusterSize) break;
            }
        }
        ++nextId;
    }
    SIMALL_LOG_INFO("Meshing",
        "Agglomerated ", nC, " cells into ", nextId, " polyhedral clusters");
    return static_cast<std::size_t>(nextId);
}

std::size_t PolyhedralAgglomeration::emit_coarse_mesh(
        const Mesh& mesh,
        const std::vector<std::int32_t>& clusterId,
        Mesh& outMesh) {
    if (clusterId.size() != mesh.cells().size()) return 0;
    const auto& srcN = mesh.nodes();
    const auto& srcF = mesh.faces();
    const auto& srcC = mesh.cells();
    auto& dN = outMesh.nodes();
    auto& dF = outMesh.faces();
    auto& dC = outMesh.cells();
    if (dF.nodeOffsets.empty()) dF.nodeOffsets.push_back(0);
    if (dC.faceOffsets.empty()) dC.faceOffsets.push_back(0);

    // Copy nodes unchanged (1:1 mapping).
    const NodeId baseId = static_cast<NodeId>(dN.size());
    for (std::size_t i = 0; i < srcN.size(); ++i) {
        dN.x.push_back(srcN.x[i]); dN.y.push_back(srcN.y[i]); dN.z.push_back(srcN.z[i]);
    }
    // Determine number of clusters & seed coarse cells.
    std::int32_t nClusters = 0;
    for (auto id : clusterId) nClusters = std::max(nClusters, id + 1);
    for (std::int32_t k = 0; k < nClusters; ++k) {
        dC.volume.push_back(0);
        dC.centroidX.push_back(0); dC.centroidY.push_back(0); dC.centroidZ.push_back(0);
    }
    std::vector<std::int32_t> coarseFaceOff(nClusters + 1, 0);

    std::int32_t cellOff = dC.faceOffsets.empty() ? 0 : dC.faceOffsets.back();

    // Walk source faces; classify based on owner/neighbour cluster.
    for (std::size_t f = 0; f < srcF.size(); ++f) {
        const auto o = srcF.owner[f];
        const auto n = srcF.neighbor[f];
        const auto co = clusterId[o];
        if (n == kBoundaryCell) {
            // External boundary face → carry into coarse mesh.
            dF.owner.push_back(static_cast<CellId>(co));
            dF.neighbor.push_back(kBoundaryCell);
            dF.areaX.push_back(srcF.areaX[f]); dF.areaY.push_back(srcF.areaY[f]);
            dF.areaZ.push_back(srcF.areaZ[f]);
            dF.centroidX.push_back(srcF.centroidX[f]);
            dF.centroidY.push_back(srcF.centroidY[f]);
            dF.centroidZ.push_back(srcF.centroidZ[f]);
            dF.boundaryZone.push_back(srcF.boundaryZone[f]);
            const auto a = srcF.nodeOffsets[f];
            const auto b = srcF.nodeOffsets[f+1];
            for (std::int32_t k = a; k < b; ++k)
                dF.nodeIndices.push_back(baseId + srcF.nodeIndices[k]);
            dF.nodeOffsets.push_back(static_cast<std::int32_t>(dF.nodeIndices.size()));
            coarseFaceOff[co] += 1;
        } else {
            const auto cn = clusterId[n];
            if (co == cn) continue;       // internal to a cluster → dropped
            dF.owner.push_back(static_cast<CellId>(co));
            dF.neighbor.push_back(static_cast<CellId>(cn));
            dF.areaX.push_back(srcF.areaX[f]); dF.areaY.push_back(srcF.areaY[f]);
            dF.areaZ.push_back(srcF.areaZ[f]);
            dF.centroidX.push_back(srcF.centroidX[f]);
            dF.centroidY.push_back(srcF.centroidY[f]);
            dF.centroidZ.push_back(srcF.centroidZ[f]);
            dF.boundaryZone.push_back(0);
            const auto a = srcF.nodeOffsets[f];
            const auto b = srcF.nodeOffsets[f+1];
            for (std::int32_t k = a; k < b; ++k)
                dF.nodeIndices.push_back(baseId + srcF.nodeIndices[k]);
            dF.nodeOffsets.push_back(static_cast<std::int32_t>(dF.nodeIndices.size()));
            coarseFaceOff[co] += 1;
            coarseFaceOff[cn] += 1;
        }
    }
    (void)cellOff;

    // Per-cell CSR for coarse mesh: scan faces and append.
    std::vector<std::int32_t> off(nClusters + 1, 0);
    for (std::int32_t k = 0; k < nClusters; ++k) off[k+1] = off[k] + coarseFaceOff[k];
    std::vector<FaceId> idxBuf(off.back());
    std::vector<std::int32_t> cursor = off;
    for (std::size_t f = 0; f < dF.size(); ++f) {
        const auto o = dF.owner[f];
        idxBuf[cursor[o]++] = static_cast<FaceId>(f);
        if (dF.neighbor[f] != kBoundaryCell) {
            const auto nb = dF.neighbor[f];
            idxBuf[cursor[nb]++] = static_cast<FaceId>(f);
        }
    }
    dC.faceOffsets.assign(off.begin(), off.end());
    dC.faceIndices.assign(idxBuf.begin(), idxBuf.end());

    SIMALL_LOG_INFO("Meshing",
        "Coarse mesh emitted: cells=", nClusters, " faces=", dF.size(),
        " nodes=", dN.size());
    // Approximate cluster volume = Σ src cell volumes.
    for (std::size_t c = 0; c < srcC.size(); ++c)
        dC.volume[clusterId[c]] += srcC.volume[c];
    return static_cast<std::size_t>(nClusters);
}

}  // namespace simall::meshing
