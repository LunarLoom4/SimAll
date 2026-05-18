// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/PartitionerScotch.cpp
// =============================================================================
#include "meshing/PartitionerScotch.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <queue>
#include <random>

namespace simall::meshing {

void PartitionerScotch::initialize(ScotchProps props) {
    p_ = props;
    SIMALL_LOG_INFO("Meshing",
        "PartitionerScotch init: nParts=", p_.nParts,
        " fmPasses=", p_.fmPasses);
}

std::size_t PartitionerScotch::partition(const Mesh& mesh,
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

void PartitionerScotch::bisect(const std::vector<std::int32_t>& xadj,
                               const std::vector<std::int32_t>& adjncy,
                               const std::vector<std::int32_t>& vmap,
                               std::vector<std::int32_t>& partOut,
                               std::int32_t partA, std::int32_t partB,
                               std::uint64_t seed) {
    if (partA == partB) {
        for (auto v : vmap) partOut[v] = partA;
        return;
    }
    const int nLocal = static_cast<int>(vmap.size());
    if (nLocal == 0) return;
    if (nLocal == 1) { partOut[vmap[0]] = partA; return; }

    // BFS-based initial bisection from two farthest seeds.
    std::mt19937_64 rng(seed);
    // Build local-vertex remap: original index → local index in this subgraph.
    // We bisect on the original global graph but only consider vertices in `vmap`.
    std::vector<char> inSet(xadj.size() - 1, 0);
    for (auto v : vmap) inSet[v] = 1;

    auto bfs_far = [&](std::int32_t start) {
        std::vector<std::int32_t> d(xadj.size()-1, -1);
        std::queue<std::int32_t> q;
        d[start] = 0; q.push(start);
        std::int32_t farV = start, farD = 0;
        while (!q.empty()) {
            const auto u = q.front(); q.pop();
            if (d[u] > farD) { farD = d[u]; farV = u; }
            for (std::int32_t k = xadj[u]; k < xadj[u+1]; ++k) {
                const auto v = adjncy[k];
                if (!inSet[v]) continue;
                if (d[v] < 0) { d[v] = d[u] + 1; q.push(v); }
            }
        }
        return farV;
    };
    const auto s0 = vmap[rng() % nLocal];
    const auto s1 = bfs_far(s0);
    const auto s2 = bfs_far(s1);

    // Two-source BFS labelling.
    std::vector<std::int32_t> lab(xadj.size()-1, -1);
    {
        std::queue<std::int32_t> q;
        lab[s1] = 0; q.push(s1);
        lab[s2] = 1; q.push(s2);
        while (!q.empty()) {
            const auto u = q.front(); q.pop();
            for (std::int32_t k = xadj[u]; k < xadj[u+1]; ++k) {
                const auto v = adjncy[k];
                if (!inSet[v] || lab[v] >= 0) continue;
                lab[v] = lab[u]; q.push(v);
            }
        }
    }
    // FM bisection refinement.
    const int target = nLocal / 2;
    const int maxImb = static_cast<int>(target * p_.imbalanceTol + 1);
    int sizeA = 0, sizeB = 0;
    for (auto v : vmap) (lab[v] == 0 ? ++sizeA : ++sizeB);
    for (int pass = 0; pass < p_.fmPasses; ++pass) {
        int moves = 0;
        for (auto u : vmap) {
            int extA = 0, extB = 0;
            for (std::int32_t k = xadj[u]; k < xadj[u+1]; ++k) {
                const auto v = adjncy[k];
                if (!inSet[v]) continue;
                if (lab[v] == 0) ++extA; else ++extB;
            }
            if (lab[u] == 0 && extB - extA > 0 && sizeB + 1 <= maxImb) {
                lab[u] = 1; --sizeA; ++sizeB; ++moves;
            } else if (lab[u] == 1 && extA - extB > 0 && sizeA + 1 <= maxImb) {
                lab[u] = 0; ++sizeA; --sizeB; ++moves;
            }
        }
        if (moves == 0) break;
    }

    // Recurse.
    const int nA = (partB - partA + 1) / 2;
    std::vector<std::int32_t> mapA, mapB;
    for (auto v : vmap) (lab[v] == 0 ? mapA : mapB).push_back(v);
    bisect(xadj, adjncy, mapA, partOut, partA,          partA + nA - 1, seed * 1103515245ULL + 12345);
    bisect(xadj, adjncy, mapB, partOut, partA + nA,     partB,          seed * 6364136223846793005ULL + 1442695040888963407ULL);
}

std::size_t PartitionerScotch::partition(const std::vector<std::int32_t>& xadj,
                                         const std::vector<std::int32_t>& adjncy,
                                         std::vector<std::int32_t>& partOut) {
    const std::int32_t nV = static_cast<std::int32_t>(xadj.size()) - 1;
    if (nV <= 0 || p_.nParts <= 1) {
        partOut.assign(std::max(0, nV), 0);
        return 0;
    }
    partOut.assign(nV, -1);
    std::vector<std::int32_t> vmap(nV);
    for (std::int32_t i = 0; i < nV; ++i) vmap[i] = i;
    bisect(xadj, adjncy, vmap, partOut, 0,
           static_cast<std::int32_t>(p_.nParts - 1), p_.rngSeed);

    // Edge cut.
    std::size_t edgeCut = 0;
    for (std::int32_t u = 0; u < nV; ++u)
        for (std::int32_t k = xadj[u]; k < xadj[u+1]; ++k)
            if (adjncy[k] > u && partOut[adjncy[k]] != partOut[u]) ++edgeCut;
    SIMALL_LOG_INFO("Meshing",
        "Scotch partition done: nParts=", p_.nParts, " edgeCut=", edgeCut);
    return edgeCut;
}

}  // namespace simall::meshing
