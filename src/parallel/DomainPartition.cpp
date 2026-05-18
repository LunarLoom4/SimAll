// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/DomainPartition.cpp
// =============================================================================
#include "parallel/DomainPartition.hpp"
#include "core/Logger.hpp"
#include "meshing/PartitionerMetis.hpp"
#include "meshing/PartitionerScotch.hpp"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace simall::parallel {

DomainPartition::DomainPartition(MpiContext& ctx) : ctx_(ctx) {}

void DomainPartition::initialize(DomainPartitionProps props) {
    p_ = props;
    if (p_.nParts <= 0) p_.nParts = ctx_.size();
    SIMALL_LOG_INFO("Parallel",
        "DomainPartition init: kind=",
        (p_.kind == PartitionerKind::Metis ? "metis" : "scotch"),
        " nParts=", p_.nParts);
}

std::size_t DomainPartition::partition(const meshing::Mesh& mesh,
                                       DomainPlan& plan) {
    const auto& Ff = mesh.faces();
    plan.cellRank.clear();
    plan.localCells.clear();
    plan.ghostCells.clear();
    plan.neighbours.clear();

    // 1) Graph partition.
    if (p_.kind == PartitionerKind::Metis) {
        meshing::PartitionerMetis pm;
        meshing::MetisProps mp;
        mp.nParts = p_.nParts;
        mp.imbalanceTol = p_.imbalanceTol;
        pm.initialize(mp);
        plan.edgeCut = pm.partition(mesh, plan.cellRank);
    } else {
        meshing::PartitionerScotch ps;
        meshing::ScotchProps sp;
        sp.nParts = p_.nParts;
        sp.imbalanceTol = p_.imbalanceTol;
        ps.initialize(sp);
        plan.edgeCut = ps.partition(mesh, plan.cellRank);
    }

    const std::int32_t myRank = ctx_.rank();

    // 2) Collect local owned cells and per-neighbour-rank cut-face cells.
    //    Map global cell id → local position in the [owned | ghost] array.
    std::unordered_map<std::int32_t, std::int32_t> g2local;
    for (std::size_t c = 0; c < plan.cellRank.size(); ++c) {
        if (plan.cellRank[c] == myRank) {
            g2local[static_cast<std::int32_t>(c)] =
                static_cast<std::int32_t>(plan.localCells.size());
            plan.localCells.push_back(static_cast<std::int32_t>(c));
        }
    }

    // Per remote rank → (set of local cells whose values are needed there,
    //                    set of remote cells we need from there).
    std::unordered_map<std::int32_t, std::vector<std::int32_t>> sendByRank;
    std::unordered_map<std::int32_t, std::vector<std::int32_t>> recvByRank;
    std::unordered_set<std::int64_t> sentPair, recvPair;
    auto pair_key = [](std::int32_t r, std::int32_t c) {
        return (static_cast<std::int64_t>(r) << 32) ^ static_cast<std::int64_t>(c);
    };

    for (std::size_t f = 0; f < Ff.size(); ++f) {
        const auto o = static_cast<std::int32_t>(Ff.owner[f]);
        const auto n = Ff.neighbor[f];
        if (n == meshing::kBoundaryCell) continue;
        const auto ni = static_cast<std::int32_t>(n);
        const auto ro = plan.cellRank[o];
        const auto rn = plan.cellRank[ni];
        if (ro == rn) continue;
        // Owner is local → its value is needed by rn's rank; we own o.
        if (ro == myRank) {
            const auto k = pair_key(rn, o);
            if (sentPair.insert(k).second) sendByRank[rn].push_back(o);
            const auto kr = pair_key(rn, ni);
            if (recvPair.insert(kr).second) recvByRank[rn].push_back(ni);
        } else if (rn == myRank) {
            const auto k = pair_key(ro, ni);
            if (sentPair.insert(k).second) sendByRank[ro].push_back(ni);
            const auto kr = pair_key(ro, o);
            if (recvPair.insert(kr).second) recvByRank[ro].push_back(o);
        }
    }

    // 3) Append ghost cells to the local layout and translate global→local.
    for (auto& [remote, recvGlobals] : recvByRank) {
        for (auto gc : recvGlobals) {
            if (g2local.find(gc) != g2local.end()) continue;
            g2local[gc] = static_cast<std::int32_t>(plan.localCells.size()
                                                  + plan.ghostCells.size());
            plan.ghostCells.push_back(gc);
        }
    }

    // 4) Materialise NeighbourComm specs.
    std::unordered_set<std::int32_t> ranks;
    for (auto& [r, _] : sendByRank) ranks.insert(r);
    for (auto& [r, _] : recvByRank) ranks.insert(r);
    for (auto r : ranks) {
        NeighbourComm nc{};
        nc.rank = r;
        for (auto gc : sendByRank[r]) nc.sendCells.push_back(g2local[gc]);
        for (auto gc : recvByRank[r]) nc.recvCells.push_back(g2local[gc]);
        std::sort(nc.sendCells.begin(), nc.sendCells.end());
        std::sort(nc.recvCells.begin(), nc.recvCells.end());
        plan.neighbours.push_back(std::move(nc));
    }

    SIMALL_LOG_INFO("Parallel",
        "DomainPartition rank=", myRank,
        " owned=", plan.localCells.size(),
        " ghost=", plan.ghostCells.size(),
        " neighbours=", plan.neighbours.size(),
        " edgeCut=", plan.edgeCut);
    return plan.localCells.size();
}

void DomainPartition::install_into(GhostExchange& gx,
                                   const DomainPlan& plan) const {
    gx.clear_neighbours();
    gx.set_layout(plan.localCells.size(), plan.ghostCells.size());
    for (const auto& nc : plan.neighbours) gx.add_neighbour(nc);
}

}  // namespace simall::parallel
