// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/PartitionerMetis.cpp
// =============================================================================
#include "meshing/PartitionerMetis.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <queue>
#include <random>
#include <vector>

namespace simall::meshing {

void PartitionerMetis::initialize(MetisProps props) {
    p_ = props;
    SIMALL_LOG_INFO("Meshing",
        "PartitionerMetis init: nParts=", p_.nParts,
        " imbalanceTol=", p_.imbalanceTol);
}

std::size_t PartitionerMetis::partition(const Mesh& mesh,
                                        std::vector<std::int32_t>& partOut) {
    const auto& C  = mesh.cells();
    const auto& Ff = mesh.faces();
    const std::int32_t nV = static_cast<std::int32_t>(C.size());

    std::vector<std::int32_t> xadj(nV + 1, 0);
    for (std::size_t f = 0; f < Ff.size(); ++f) {
        if (Ff.neighbor[f] == kBoundaryCell) continue;
        ++xadj[Ff.owner[f] + 1];
        ++xadj[Ff.neighbor[f] + 1];
    }
    for (std::int32_t i = 0; i < nV; ++i) xadj[i+1] += xadj[i];
    std::vector<std::int32_t> adjncy(xadj.back());
    std::vector<std::int32_t> cursor(xadj);
    for (std::size_t f = 0; f < Ff.size(); ++f) {
        if (Ff.neighbor[f] == kBoundaryCell) continue;
        const std::int32_t o = static_cast<std::int32_t>(Ff.owner[f]);
        const std::int32_t n = static_cast<std::int32_t>(Ff.neighbor[f]);
        adjncy[cursor[o]++] = n;
        adjncy[cursor[n]++] = o;
    }
    return partition(xadj, adjncy, partOut);
}

std::size_t PartitionerMetis::partition(const std::vector<std::int32_t>& xadj,
                                        const std::vector<std::int32_t>& adjncy,
                                        std::vector<std::int32_t>& partOut) {
    const std::int32_t nV = static_cast<std::int32_t>(xadj.size()) - 1;
    if (nV <= 0 || p_.nParts <= 1) {
        partOut.assign(std::max(0, nV), 0);
        return 0;
    }
    partOut.assign(nV, 0);

    // -------- Phase 1 — initial growth-region partitioning. --------
    // Pick nParts seed vertices spaced ~ uniformly via BFS-distance heuristic.
    std::mt19937_64 rng(p_.rngSeed);
    std::vector<std::int32_t> seeds;
    seeds.push_back(static_cast<std::int32_t>(rng() % nV));
    while (static_cast<int>(seeds.size()) < p_.nParts) {
        // farthest-point seed.
        std::vector<std::int32_t> d(nV, -1);
        std::queue<std::int32_t> q;
        for (auto s : seeds) { d[s] = 0; q.push(s); }
        std::int32_t farV = seeds.front(), farD = -1;
        while (!q.empty()) {
            const auto u = q.front(); q.pop();
            if (d[u] > farD) { farD = d[u]; farV = u; }
            for (std::int32_t k = xadj[u]; k < xadj[u+1]; ++k) {
                const auto v = adjncy[k];
                if (d[v] < 0) { d[v] = d[u] + 1; q.push(v); }
            }
        }
        seeds.push_back(farV);
    }
    // Multi-source BFS labelling.
    {
        std::vector<std::int32_t> d(nV, std::numeric_limits<std::int32_t>::max());
        std::queue<std::int32_t> q;
        for (int p = 0; p < p_.nParts; ++p) {
            d[seeds[p]] = 0; partOut[seeds[p]] = p; q.push(seeds[p]);
        }
        while (!q.empty()) {
            const auto u = q.front(); q.pop();
            for (std::int32_t k = xadj[u]; k < xadj[u+1]; ++k) {
                const auto v = adjncy[k];
                if (d[v] > d[u] + 1) { d[v] = d[u] + 1;
                                       partOut[v] = partOut[u]; q.push(v); }
            }
        }
    }

    // -------- Phase 2 — Kernighan-Lin / FM refinement passes. --------
    std::vector<std::int32_t> partSize(p_.nParts, 0);
    for (auto p : partOut) ++partSize[p];
    const std::int32_t target = nV / p_.nParts;
    const std::int32_t maxBalanced =
        static_cast<std::int32_t>(target * p_.imbalanceTol + 1.0);

    for (int sweep = 0; sweep < p_.refineSweeps; ++sweep) {
        std::int32_t moves = 0;
        for (std::int32_t u = 0; u < nV; ++u) {
            // External-internal gain per other-part.
            std::vector<std::int32_t> ext(p_.nParts, 0);
            for (std::int32_t k = xadj[u]; k < xadj[u+1]; ++k)
                ++ext[partOut[adjncy[k]]];
            const auto cur = partOut[u];
            std::int32_t bestP = cur;
            std::int32_t bestGain = 0;
            for (int p = 0; p < p_.nParts; ++p) {
                if (p == cur) continue;
                if (partSize[p] + 1 > maxBalanced) continue;
                const std::int32_t gain = ext[p] - ext[cur];
                if (gain > bestGain) { bestGain = gain; bestP = p; }
            }
            if (bestP != cur && bestGain > 0) {
                --partSize[cur]; ++partSize[bestP];
                partOut[u] = bestP; ++moves;
            }
        }
        if (moves == 0) break;
    }

    // Edge cut.
    std::size_t edgeCut = 0;
    for (std::int32_t u = 0; u < nV; ++u)
        for (std::int32_t k = xadj[u]; k < xadj[u+1]; ++k)
            if (adjncy[k] > u && partOut[adjncy[k]] != partOut[u]) ++edgeCut;
    SIMALL_LOG_INFO("Meshing",
        "Metis partition done: nParts=", p_.nParts, " edgeCut=", edgeCut);
    return edgeCut;
}

}  // namespace simall::meshing
