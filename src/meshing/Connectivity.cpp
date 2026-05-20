// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/Connectivity.cpp
// =============================================================================
#include "meshing/Connectivity.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <unordered_map>

namespace simall::meshing
{

namespace
{

/// Order-independent hash of a face's vertex list.
struct FaceKey
{
    std::vector<NodeId> sorted;
    bool operator==(const FaceKey& o) const noexcept { return sorted == o.sorted; }
};
struct FaceKeyHash
{
    std::size_t operator()(const FaceKey& k) const noexcept
    {
        std::size_t h = 1469598103934665603ull;
        for (NodeId v : k.sorted) {
            h ^= v + 0x9e3779b97f4a7c15ull + (h << 12) + (h >> 4);
        }
        return h;
    }
};

FaceKey make_key(const std::vector<NodeId>& nodes)
{
    FaceKey k{nodes};
    std::sort(k.sorted.begin(), k.sorted.end());
    return k;
}

} // namespace

void ConnectivityBuilder::build(Mesh& mesh,
                                const NodeStorage& nodes,
                                const std::vector<CellDescriptor>& cells)
{
    auto& N = mesh.nodes();
    auto& F = mesh.faces();
    auto& C = mesh.cells();

    // ---- nodes (move-copy) --------------------------------------------------
    N.x = nodes.x;
    N.y = nodes.y;
    N.z = nodes.z;

    // ---- temporary CSR builders --------------------------------------------
    const std::size_t nCells = cells.size();
    std::vector<std::vector<FaceId>> cellFaces(nCells);

    std::unordered_map<FaceKey, FaceId, FaceKeyHash> faceMap;
    faceMap.reserve(nCells * 6);

    auto pushFace = [&](FaceId f, const std::vector<NodeId>& v) {
        F.nodeOffsets.push_back(static_cast<std::int32_t>(F.nodeIndices.size()));
        for (NodeId n : v)
            F.nodeIndices.push_back(n);
        F.owner.push_back(static_cast<CellId>(f)); // overwritten below
        F.neighbor.push_back(kBoundaryCell);
        F.boundaryZone.push_back(0);
    };

    for (std::size_t c = 0; c < nCells; ++c) {
        for (const auto& v : cells[c].faces) {
            FaceKey k = make_key(v);
            auto it = faceMap.find(k);
            if (it == faceMap.end()) {
                FaceId fid = F.owner.size();
                pushFace(fid, v);
                F.owner.back() = c;
                faceMap.emplace(std::move(k), fid);
                cellFaces[c].push_back(fid);
            } else {
                FaceId fid = it->second;
                F.neighbor[fid] = c;
                cellFaces[c].push_back(fid);
            }
        }
    }
    F.nodeOffsets.push_back(static_cast<std::int32_t>(F.nodeIndices.size()));

    // ---- per-cell face CSR --------------------------------------------------
    C.faceOffsets.clear();
    C.faceIndices.clear();
    C.faceOffsets.push_back(0);
    for (std::size_t c = 0; c < nCells; ++c) {
        for (FaceId f : cellFaces[c])
            C.faceIndices.push_back(f);
        C.faceOffsets.push_back(static_cast<std::int32_t>(C.faceIndices.size()));
    }
    C.volume.assign(nCells, 0);
    C.centroidX.assign(nCells, 0);
    C.centroidY.assign(nCells, 0);
    C.centroidZ.assign(nCells, 0);

    // ---- mark boundary zones (interior=0, exterior=1 by default) -----------
    for (std::size_t f = 0; f < F.size(); ++f)
        if (F.neighbor[f] == kBoundaryCell)
            F.boundaryZone[f] = 1;

    mesh.compute_geometry();

    SIMALL_LOG_INFO("Mesh", "Connectivity built: ", nCells, " cells, ", F.size(), " faces");
}

} // namespace simall::meshing
