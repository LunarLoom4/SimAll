// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MultiblockHex.cpp
// =============================================================================
#include "meshing/MultiblockHex.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace simall::meshing {

namespace {
// Trilinear hex interpolation: corners 0..7 in standard ordering,
//   0: (0,0,0)  1: (1,0,0)  2: (1,1,0)  3: (0,1,0)
//   4: (0,0,1)  5: (1,0,1)  6: (1,1,1)  7: (0,1,1)
util::Vec3d trilinear(const std::array<util::Vec3d,8>& c,
                      double u, double v, double w) {
    auto lerp = [](const util::Vec3d& a, const util::Vec3d& b, double t)
    { return util::Vec3d{a.x+(b.x-a.x)*t, a.y+(b.y-a.y)*t, a.z+(b.z-a.z)*t}; };
    const auto e0 = lerp(c[0], c[1], u);
    const auto e1 = lerp(c[3], c[2], u);
    const auto e2 = lerp(c[4], c[5], u);
    const auto e3 = lerp(c[7], c[6], u);
    const auto f0 = lerp(e0, e1, v);
    const auto f1 = lerp(e2, e3, v);
    return lerp(f0, f1, w);
}

struct Key { long long ix, iy, iz; };
struct KeyHash { size_t operator()(const Key& k) const noexcept
    { return std::hash<long long>{}(k.ix) ^ (std::hash<long long>{}(k.iy)<<1)
           ^ (std::hash<long long>{}(k.iz)<<2); } };
struct KeyEq { bool operator()(const Key& a, const Key& b) const noexcept
    { return a.ix == b.ix && a.iy == b.iy && a.iz == b.iz; } };
}  // namespace

void MultiblockHex::initialize(MultiblockProps props) {
    p_ = props;
    SIMALL_LOG_INFO("Meshing",
        "MultiblockHex init: tol=", p_.matchTol, " reorder=", p_.reorder);
}

std::size_t MultiblockHex::assemble(Mesh& outMesh) {
    if (blocks_.empty()) return 0;
    auto& N = outMesh.nodes();
    auto& F = outMesh.faces();
    auto& C = outMesh.cells();
    if (F.nodeOffsets.empty()) F.nodeOffsets.push_back(0);
    if (C.faceOffsets.empty()) C.faceOffsets.push_back(0);

    // Global node registry keyed by quantised coordinate.
    std::unordered_map<Key, NodeId, KeyHash, KeyEq> registry;
    const double inv = 1.0 / std::max(p_.matchTol, 1e-12);
    auto quantise = [&](const util::Vec3d& p) -> Key {
        return { std::llround(p.x * inv),
                 std::llround(p.y * inv),
                 std::llround(p.z * inv) };
    };
    auto get_or_add = [&](const util::Vec3d& p) -> NodeId {
        const Key k = quantise(p);
        auto it = registry.find(k);
        if (it != registry.end()) return it->second;
        const NodeId id = static_cast<NodeId>(N.size());
        N.x.push_back(p.x); N.y.push_back(p.y); N.z.push_back(p.z);
        registry.emplace(k, id);
        return id;
    };

    // Face registry (for matching cross-block faces): keyed by sorted node id quartet.
    struct FaceKey { std::array<NodeId,4> n; };
    struct FaceHash { size_t operator()(const FaceKey& f) const noexcept
        { size_t h=0; for (auto x : f.n) h ^= std::hash<NodeId>{}(x) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
          return h; } };
    struct FaceEq { bool operator()(const FaceKey& a, const FaceKey& b) const noexcept
        { return a.n == b.n; } };
    std::unordered_map<FaceKey, std::size_t, FaceHash, FaceEq> faceMap;

    // CellOff tracker.
    std::int32_t cellOff = C.faceOffsets.empty() ? 0 : C.faceOffsets.back();
    std::size_t hexCount = 0;

    auto add_or_join_face = [&](NodeId a, NodeId b, NodeId c, NodeId d,
                                CellId owner, ZoneId zone) {
        FaceKey key{{a,b,c,d}};
        std::sort(key.n.begin(), key.n.end());
        auto it = faceMap.find(key);
        if (it != faceMap.end()) {
            // Join: set neighbour, clear boundary tag.
            const std::size_t fi = it->second;
            F.neighbor[fi] = owner;
            F.boundaryZone[fi] = 0;
            C.faceIndices.push_back(static_cast<FaceId>(fi));
            return;
        }
        const std::size_t fi = F.size();
        F.owner.push_back(owner);
        F.neighbor.push_back(kBoundaryCell);
        F.areaX.push_back(0); F.areaY.push_back(0); F.areaZ.push_back(0);
        F.centroidX.push_back(0); F.centroidY.push_back(0); F.centroidZ.push_back(0);
        F.boundaryZone.push_back(zone);
        F.nodeIndices.push_back(a); F.nodeIndices.push_back(b);
        F.nodeIndices.push_back(c); F.nodeIndices.push_back(d);
        F.nodeOffsets.push_back(static_cast<std::int32_t>(F.nodeIndices.size()));
        faceMap.emplace(key, fi);
        C.faceIndices.push_back(static_cast<FaceId>(fi));
    };

    for (const auto& blk : blocks_) {
        const auto Ni = blk.divisions[0];
        const auto Nj = blk.divisions[1];
        const auto Nk = blk.divisions[2];
        if (Ni == 0 || Nj == 0 || Nk == 0) continue;

        // Parametric coordinates in [0,1] for each direction.  When the
        // simpleGrading expansion ratio G = (last cell length)/(first cell
        // length) is 1, the distribution is uniform.  Otherwise cells follow
        // a geometric progression with per-cell ratio r = G^(1/(N-1)), and
        //   u[i] = (r^i - 1) / (r^N - 1)
        // recovers the cumulative parametric position.
        auto make_param = [](std::uint32_t N, double G) {
            std::vector<double> u(N + 1);
            if (N == 0) { u[0] = 0.0; return u; }
            if (N == 1 || !(G > 0.0) || std::abs(G - 1.0) < 1e-12) {
                for (std::uint32_t i = 0; i <= N; ++i)
                    u[i] = double(i) / double(N);
                return u;
            }
            const double r     = std::pow(G, 1.0 / double(N - 1));
            const double denom = std::pow(r, double(N)) - 1.0;
            u[0] = 0.0;
            for (std::uint32_t i = 1; i < N; ++i)
                u[i] = (std::pow(r, double(i)) - 1.0) / denom;
            u[N] = 1.0;
            return u;
        };
        const std::vector<double> uPar = make_param(Ni, blk.grading[0]);
        const std::vector<double> vPar = make_param(Nj, blk.grading[1]);
        const std::vector<double> wPar = make_param(Nk, blk.grading[2]);

        // Build node-id table indexed by (i,j,k) ∈ [0..Ni]·[0..Nj]·[0..Nk].
        std::vector<NodeId> nid((Ni+1)*(Nj+1)*(Nk+1));
        auto NIDX = [&](std::uint32_t i, std::uint32_t j, std::uint32_t k){
            return ((k * (Nj+1)) + j) * (Ni+1) + i;
        };
        for (std::uint32_t k = 0; k <= Nk; ++k)
        for (std::uint32_t j = 0; j <= Nj; ++j)
        for (std::uint32_t i = 0; i <= Ni; ++i) {
            nid[NIDX(i,j,k)] = get_or_add(
                trilinear(blk.corners, uPar[i], vPar[j], wPar[k]));
        }
        // Build hex cells.
        for (std::uint32_t k = 0; k < Nk; ++k)
        for (std::uint32_t j = 0; j < Nj; ++j)
        for (std::uint32_t i = 0; i < Ni; ++i) {
            const CellId cellId = static_cast<CellId>(C.size());
            const NodeId n0 = nid[NIDX(i,   j,   k  )];
            const NodeId n1 = nid[NIDX(i+1, j,   k  )];
            const NodeId n2 = nid[NIDX(i+1, j+1, k  )];
            const NodeId n3 = nid[NIDX(i,   j+1, k  )];
            const NodeId n4 = nid[NIDX(i,   j,   k+1)];
            const NodeId n5 = nid[NIDX(i+1, j,   k+1)];
            const NodeId n6 = nid[NIDX(i+1, j+1, k+1)];
            const NodeId n7 = nid[NIDX(i,   j+1, k+1)];
            C.volume.push_back(0); C.centroidX.push_back(0);
            C.centroidY.push_back(0); C.centroidZ.push_back(0);
            const std::int32_t cellStart = cellOff;
            // Six faces with zone tags from block-face spec on boundary i/j/k extremes.
            const ZoneId zImin = (i==0)    ? blk.faceZones[0] : 0;
            const ZoneId zImax = (i==Ni-1) ? blk.faceZones[1] : 0;
            const ZoneId zJmin = (j==0)    ? blk.faceZones[2] : 0;
            const ZoneId zJmax = (j==Nj-1) ? blk.faceZones[3] : 0;
            const ZoneId zKmin = (k==0)    ? blk.faceZones[4] : 0;
            const ZoneId zKmax = (k==Nk-1) ? blk.faceZones[5] : 0;
            add_or_join_face(n0, n3, n7, n4, cellId, zImin);
            add_or_join_face(n1, n5, n6, n2, cellId, zImax);
            add_or_join_face(n0, n4, n5, n1, cellId, zJmin);
            add_or_join_face(n3, n2, n6, n7, cellId, zJmax);
            add_or_join_face(n0, n1, n2, n3, cellId, zKmin);
            add_or_join_face(n4, n7, n6, n5, cellId, zKmax);
            cellOff = cellStart + 6;
            C.faceOffsets.push_back(cellOff);
            ++hexCount;
        }
    }
    SIMALL_LOG_INFO("Meshing",
        "Multiblock assembled: blocks=", blocks_.size(),
        " hexes=", hexCount, " nodes=", N.size(), " faces=", F.size());
    return hexCount;
}

}  // namespace simall::meshing
