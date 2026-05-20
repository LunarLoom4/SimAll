// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/BoundaryRecovery.cpp
// =============================================================================
#include "meshing/BoundaryRecovery.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <array>
#include <unordered_set>

namespace simall::meshing
{

namespace
{
struct PairHash
{
    size_t operator()(const std::pair<NodeId, NodeId>& p) const noexcept
    {
        return std::hash<NodeId>{}(p.first) ^ (std::hash<NodeId>{}(p.second) << 1);
    }
};
struct TripleHash
{
    size_t operator()(const std::array<NodeId, 3>& t) const noexcept
    {
        return std::hash<NodeId>{}(t[0]) ^ (std::hash<NodeId>{}(t[1]) << 1)
               ^ (std::hash<NodeId>{}(t[2]) << 2);
    }
};
} // namespace

void BoundaryRecovery::initialize(BoundaryRecoveryProps props)
{
    p_ = props;
    SIMALL_LOG_INFO("Meshing", "BoundaryRecovery init: maxSteiner=", p_.maxSteinerInsertions);
}

bool BoundaryRecovery::edge_present(const Mesh& m, const RequiredEdge& e) const
{
    // Check face nodeIndices for both endpoints (sufficient for tetra mesh:
    // every edge lies on at least one face).
    const auto& F = m.faces();
    for (std::size_t f = 0; f < F.size(); ++f) {
        const auto a = F.nodeOffsets[f];
        const auto b = F.nodeOffsets[f + 1];
        bool hasA = false, hasB = false;
        for (std::int32_t k = a; k < b; ++k) {
            if (F.nodeIndices[k] == e.a)
                hasA = true;
            if (F.nodeIndices[k] == e.b)
                hasB = true;
        }
        if (hasA && hasB)
            return true;
    }
    return false;
}

bool BoundaryRecovery::face_present(const Mesh& m, const RequiredFace& f) const
{
    const auto& F = m.faces();
    std::array<NodeId, 3> target{f.a, f.b, f.c};
    std::sort(target.begin(), target.end());
    for (std::size_t i = 0; i < F.size(); ++i) {
        const auto a = F.nodeOffsets[i];
        const auto b = F.nodeOffsets[i + 1];
        if (b - a < 3)
            continue;
        std::vector<NodeId> nodes;
        for (std::int32_t k = a; k < b; ++k)
            nodes.push_back(F.nodeIndices[k]);
        std::sort(nodes.begin(), nodes.end());
        if (std::includes(nodes.begin(), nodes.end(), target.begin(), target.end()))
            return true;
    }
    return false;
}

NodeId BoundaryRecovery::insert_steiner(Mesh& mesh, const RequiredEdge& e)
{
    auto& N = mesh.nodes();
    const double mx = 0.5 * (N.x[e.a] + N.x[e.b]);
    const double my = 0.5 * (N.y[e.a] + N.y[e.b]);
    const double mz = 0.5 * (N.z[e.a] + N.z[e.b]);
    const NodeId nid = static_cast<NodeId>(N.size());
    N.x.push_back(mx);
    N.y.push_back(my);
    N.z.push_back(mz);
    return nid;
}

RecoveryReport BoundaryRecovery::recover(Mesh& mesh)
{
    RecoveryReport rep{};
    rep.edgesMissing = edges_.size();
    rep.facesMissing = faces_.size();

    // Phase 1: edge recovery.
    std::vector<RequiredEdge> pending = edges_;
    while (!pending.empty() && rep.steinerAdded < p_.maxSteinerInsertions) {
        std::vector<RequiredEdge> nextRound;
        for (const auto& e : pending) {
            if (edge_present(mesh, e)) {
                ++rep.edgesRecovered;
                continue;
            }
            // Subdivide: insert Steiner point and replace with two sub-edges.
            const NodeId mid = insert_steiner(mesh, e);
            ++rep.steinerAdded;
            nextRound.push_back({e.a, mid});
            nextRound.push_back({mid, e.b});
        }
        if (nextRound.size() == pending.size())
            break; // no progress
        pending = std::move(nextRound);
    }

    // Phase 2: face recovery (best-effort: subdivide a missing face into
    // 3 sub-faces by centroid Steiner point).
    std::vector<RequiredFace> pendingF = faces_;
    while (!pendingF.empty() && rep.steinerAdded < p_.maxSteinerInsertions) {
        std::vector<RequiredFace> nextRound;
        for (const auto& f : pendingF) {
            if (face_present(mesh, f)) {
                ++rep.facesRecovered;
                continue;
            }
            auto& N = mesh.nodes();
            const double cx = (N.x[f.a] + N.x[f.b] + N.x[f.c]) / 3.0;
            const double cy = (N.y[f.a] + N.y[f.b] + N.y[f.c]) / 3.0;
            const double cz = (N.z[f.a] + N.z[f.b] + N.z[f.c]) / 3.0;
            const NodeId c = static_cast<NodeId>(N.size());
            N.x.push_back(cx);
            N.y.push_back(cy);
            N.z.push_back(cz);
            ++rep.steinerAdded;
            nextRound.push_back({f.a, f.b, c});
            nextRound.push_back({f.b, f.c, c});
            nextRound.push_back({f.c, f.a, c});
        }
        if (nextRound.size() == pendingF.size())
            break;
        pendingF = std::move(nextRound);
    }

    SIMALL_LOG_INFO("Meshing",
                    "BoundaryRecovery: edges ",
                    rep.edgesRecovered,
                    "/",
                    rep.edgesMissing,
                    " faces ",
                    rep.facesRecovered,
                    "/",
                    rep.facesMissing,
                    " Steiner=",
                    rep.steinerAdded);
    return rep;
}

} // namespace simall::meshing
