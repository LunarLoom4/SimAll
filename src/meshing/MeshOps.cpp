// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MeshOps.cpp
// Phase  : 23 Pass 11
// =============================================================================

#include "meshing/MeshOps.hpp"

#include "core/Logger.hpp"
#include "meshing/Connectivity.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::meshing::ops
{

// Sentinel for "no remap" / "uninitialised" slots inside compaction passes.
// NOTE: util::kInvalidId is 0, but 0 is a valid mesh element index, so we
// cannot use it here. Use the max value of the id type instead.
inline constexpr NodeId kInvalidId = static_cast<NodeId>(-1);

// =============================================================================
// Affine helpers
// =============================================================================
double Affine::det() const noexcept
{
    const auto& r = R;
    return r[0] * (r[4] * r[8] - r[5] * r[7]) - r[1] * (r[3] * r[8] - r[5] * r[6])
           + r[2] * (r[3] * r[7] - r[4] * r[6]);
}

namespace
{

inline util::Vec3d apply(const Affine& a, double x, double y, double z) noexcept
{
    return util::Vec3d{a.R[0] * x + a.R[1] * y + a.R[2] * z + a.t.x,
                       a.R[3] * x + a.R[4] * y + a.R[5] * z + a.t.y,
                       a.R[6] * x + a.R[7] * y + a.R[8] * z + a.t.z};
}

void reverse_face_orientations(Mesh& m)
{
    auto& f = m.faces();
    const std::size_t nF = f.size();
    for (std::size_t i = 0; i < nF; ++i) {
        const std::int32_t b = f.nodeOffsets[i];
        const std::int32_t e = f.nodeOffsets[i + 1];
        if (e - b < 3)
            continue;
        std::reverse(f.nodeIndices.begin() + b, f.nodeIndices.begin() + e);
    }
}

} // namespace

// =============================================================================
// Geometry transforms
// =============================================================================
void transform_points(Mesh& m, const Affine& a)
{
    auto& n = m.nodes();
    const std::size_t N = n.size();
    for (std::size_t i = 0; i < N; ++i) {
        const auto p = apply(a, n.x[i], n.y[i], n.z[i]);
        n.x[i] = p.x;
        n.y[i] = p.y;
        n.z[i] = p.z;
    }
    if (a.det() < 0.0)
        reverse_face_orientations(m);
    m.compute_geometry();
}

void translate(Mesh& m, double dx, double dy, double dz)
{
    Affine a;
    a.t = {dx, dy, dz};
    transform_points(m, a);
}

void rotate(Mesh& m, util::Vec3d axis, double angleRad, util::Vec3d pivot)
{
    // Rodrigues rotation about a unit axis through pivot.
    const double n = axis.norm();
    if (n <= 0.0)
        return;
    const util::Vec3d k{axis.x / n, axis.y / n, axis.z / n};
    const double c = std::cos(angleRad);
    const double s = std::sin(angleRad);
    const double C = 1.0 - c;

    Affine A;
    A.R = {c + k.x * k.x * C,
           k.x * k.y * C - k.z * s,
           k.x * k.z * C + k.y * s,
           k.y * k.x * C + k.z * s,
           c + k.y * k.y * C,
           k.y * k.z * C - k.x * s,
           k.z * k.x * C - k.y * s,
           k.z * k.y * C + k.x * s,
           c + k.z * k.z * C};
    // Pre/post translate to rotate about pivot: x' = R*(x - p) + p
    A.t = {pivot.x - (A.R[0] * pivot.x + A.R[1] * pivot.y + A.R[2] * pivot.z),
           pivot.y - (A.R[3] * pivot.x + A.R[4] * pivot.y + A.R[5] * pivot.z),
           pivot.z - (A.R[6] * pivot.x + A.R[7] * pivot.y + A.R[8] * pivot.z)};
    transform_points(m, A);
}

void mirror(Mesh& m, util::Vec3d normal, util::Vec3d p0)
{
    const double n2 = normal.dot(normal);
    if (n2 <= 0.0)
        return;
    const util::Vec3d nu{
        normal.x / std::sqrt(n2), normal.y / std::sqrt(n2), normal.z / std::sqrt(n2)};
    // Householder reflection: R = I - 2 n n^T ;  t = 2 (n . p0) n
    Affine A;
    A.R = {1.0 - 2.0 * nu.x * nu.x,
           -2.0 * nu.x * nu.y,
           -2.0 * nu.x * nu.z,
           -2.0 * nu.y * nu.x,
           1.0 - 2.0 * nu.y * nu.y,
           -2.0 * nu.y * nu.z,
           -2.0 * nu.z * nu.x,
           -2.0 * nu.z * nu.y,
           1.0 - 2.0 * nu.z * nu.z};
    const double d = 2.0 * (nu.x * p0.x + nu.y * p0.y + nu.z * p0.z);
    A.t = {d * nu.x, d * nu.y, d * nu.z};
    transform_points(m, A);
}

// =============================================================================
// Topology composition: merge_meshes
// =============================================================================
namespace
{

template <class V, class W> void append_vec(V& dst, const W& src)
{
    dst.insert(dst.end(), src.begin(), src.end());
}

void concat_csr_offsets(std::vector<std::int32_t>& outBuf,
                        const util::aligned_vector<std::int32_t>& aOff,
                        const util::aligned_vector<std::int32_t>& bOff,
                        std::int32_t shift)
{
    outBuf.clear();
    outBuf.reserve(aOff.size() + bOff.size() - 1);
    outBuf.insert(outBuf.end(), aOff.begin(), aOff.end());
    for (std::size_t i = 1; i < bOff.size(); ++i) {
        outBuf.push_back(bOff[i] + shift);
    }
}

} // namespace

void merge_meshes(const Mesh& a, const Mesh& b, Mesh& out)
{
    // Snapshot of a's sizes BEFORE we touch out (handles the &out == &a alias).
    const auto& aN = a.nodes();
    const auto& aF = a.faces();
    const auto& aC = a.cells();
    const auto& bN = b.nodes();
    const auto& bF = b.faces();
    const auto& bC = b.cells();

    NodeStorage nOut;
    FaceStorage fOut;
    CellStorage cOut;

    const std::size_t nA_nodes = aN.size();
    const std::size_t nA_faces = aF.size();
    const std::size_t nA_cells = aC.size();
    const std::size_t nA_faceNodes = aF.nodeIndices.size();
    const std::size_t nA_cellFaces = aC.faceIndices.size();

    // ---- Nodes -----------------------------------------------------------
    nOut.reserve(nA_nodes + bN.size());
    append_vec(nOut.x, aN.x);
    append_vec(nOut.x, bN.x);
    append_vec(nOut.y, aN.y);
    append_vec(nOut.y, bN.y);
    append_vec(nOut.z, aN.z);
    append_vec(nOut.z, bN.z);

    // ---- Faces -----------------------------------------------------------
    fOut.owner.reserve(nA_faces + bF.size());
    fOut.neighbor.reserve(nA_faces + bF.size());
    fOut.boundaryZone.reserve(nA_faces + bF.size());
    append_vec(fOut.owner, aF.owner);
    append_vec(fOut.neighbor, aF.neighbor);
    append_vec(fOut.boundaryZone, aF.boundaryZone);
    for (auto v : bF.owner) {
        fOut.owner.push_back(v + nA_cells);
    }
    for (auto v : bF.neighbor) {
        fOut.neighbor.push_back(v == kBoundaryCell ? kBoundaryCell : v + nA_cells);
    }
    append_vec(fOut.boundaryZone, bF.boundaryZone);

    // face geometry will be recomputed but reserve to keep arrays sized
    fOut.areaX.resize(fOut.owner.size(), 0.0);
    fOut.areaY.resize(fOut.owner.size(), 0.0);
    fOut.areaZ.resize(fOut.owner.size(), 0.0);
    fOut.centroidX.resize(fOut.owner.size(), 0.0);
    fOut.centroidY.resize(fOut.owner.size(), 0.0);
    fOut.centroidZ.resize(fOut.owner.size(), 0.0);

    // CSR: face -> nodes
    fOut.nodeIndices.reserve(nA_faceNodes + bF.nodeIndices.size());
    append_vec(fOut.nodeIndices, aF.nodeIndices);
    for (auto v : bF.nodeIndices)
        fOut.nodeIndices.push_back(v + nA_nodes);

    {
        std::vector<std::int32_t> tmp;
        concat_csr_offsets(
            tmp, aF.nodeOffsets, bF.nodeOffsets, static_cast<std::int32_t>(nA_faceNodes));
        fOut.nodeOffsets.assign(tmp.begin(), tmp.end());
    }

    // ---- Cells -----------------------------------------------------------
    cOut.volume.resize(nA_cells + bC.size(), 0.0);
    cOut.centroidX.resize(nA_cells + bC.size(), 0.0);
    cOut.centroidY.resize(nA_cells + bC.size(), 0.0);
    cOut.centroidZ.resize(nA_cells + bC.size(), 0.0);

    cOut.faceIndices.reserve(nA_cellFaces + bC.faceIndices.size());
    append_vec(cOut.faceIndices, aC.faceIndices);
    for (auto v : bC.faceIndices)
        cOut.faceIndices.push_back(v + nA_faces);

    {
        std::vector<std::int32_t> tmp;
        concat_csr_offsets(
            tmp, aC.faceOffsets, bC.faceOffsets, static_cast<std::int32_t>(nA_cellFaces));
        cOut.faceOffsets.assign(tmp.begin(), tmp.end());
    }

    // ---- Commit to out ---------------------------------------------------
    out.nodes() = std::move(nOut);
    out.faces() = std::move(fOut);
    out.cells() = std::move(cOut);
    out.compute_geometry();
}

// =============================================================================
// stitch_meshes
// =============================================================================
namespace
{

struct GridKey
{
    std::int64_t i, j, k;
    bool operator==(const GridKey& o) const noexcept { return i == o.i && j == o.j && k == o.k; }
};
struct GridKeyHash
{
    std::size_t operator()(const GridKey& g) const noexcept
    {
        // Cantor-style mix
        auto h1 = std::hash<std::int64_t>{}(g.i);
        auto h2 = std::hash<std::int64_t>{}(g.j);
        auto h3 = std::hash<std::int64_t>{}(g.k);
        std::size_t r = h1;
        r ^= h2 + 0x9e3779b97f4a7c15ULL + (r << 6) + (r >> 2);
        r ^= h3 + 0x9e3779b97f4a7c15ULL + (r << 6) + (r >> 2);
        return r;
    }
};

struct NodeKeyHash
{
    std::size_t operator()(const std::vector<NodeId>& v) const noexcept
    {
        std::size_t r = 1469598103934665603ULL;
        for (auto n : v) {
            r ^= std::hash<NodeId>{}(n);
            r *= 1099511628211ULL;
        }
        return r;
    }
};

} // namespace

StitchStats stitch_meshes(const Mesh& a, const Mesh& b, Mesh& out, StitchOptions opt)
{
    StitchStats stats{};

    merge_meshes(a, b, out);

    auto& nodes = out.nodes();
    auto& faces = out.faces();
    auto& cells = out.cells();

    const std::size_t N = nodes.size();
    const double tol = opt.weldTolerance;
    const double tol2 = tol * tol;
    const double inv = (tol > 0.0) ? 1.0 / tol : 1.0;

    // -------- Build spatial hash over all nodes ---------------------------
    std::unordered_map<GridKey, std::vector<NodeId>, GridKeyHash> grid;
    grid.reserve(N * 2);
    for (NodeId i = 0; i < N; ++i) {
        GridKey g{static_cast<std::int64_t>(std::floor(nodes.x[i] * inv)),
                  static_cast<std::int64_t>(std::floor(nodes.y[i] * inv)),
                  static_cast<std::int64_t>(std::floor(nodes.z[i] * inv))};
        grid[g].push_back(i);
    }

    // -------- Find canonical (lowest index) representative per node ------
    std::vector<NodeId> rep(N);
    std::iota(rep.begin(), rep.end(), NodeId{0});

    auto dist2 = [&](NodeId a_, NodeId b_) {
        const double dx = nodes.x[a_] - nodes.x[b_];
        const double dy = nodes.y[a_] - nodes.y[b_];
        const double dz = nodes.z[a_] - nodes.z[b_];
        return dx * dx + dy * dy + dz * dz;
    };

    for (NodeId i = 0; i < N; ++i) {
        if (rep[i] != i)
            continue;
        const GridKey g{static_cast<std::int64_t>(std::floor(nodes.x[i] * inv)),
                        static_cast<std::int64_t>(std::floor(nodes.y[i] * inv)),
                        static_cast<std::int64_t>(std::floor(nodes.z[i] * inv))};
        for (std::int64_t dk = -1; dk <= 1; ++dk)
            for (std::int64_t dj = -1; dj <= 1; ++dj)
                for (std::int64_t di = -1; di <= 1; ++di) {
                    const GridKey gn{g.i + di, g.j + dj, g.k + dk};
                    auto it = grid.find(gn);
                    if (it == grid.end())
                        continue;
                    for (NodeId j : it->second) {
                        if (j <= i)
                            continue;
                        if (rep[j] != j)
                            continue;
                        if (dist2(i, j) <= tol2)
                            rep[j] = i;
                    }
                }
    }

    // Compact node representatives -> new indices
    std::vector<NodeId> newIdx(N, kInvalidId);
    {
        // path compress (rep[rep[..]])
        for (NodeId i = 0; i < N; ++i) {
            NodeId r = i;
            while (rep[r] != r)
                r = rep[r];
            rep[i] = r;
        }
        NodeId cnt = 0;
        for (NodeId i = 0; i < N; ++i) {
            if (rep[i] == i)
                newIdx[i] = cnt++;
        }
        for (NodeId i = 0; i < N; ++i) {
            if (rep[i] != i)
                newIdx[i] = newIdx[rep[i]];
        }
        stats.weldedNodes = N - cnt;

        // Compact node arrays
        NodeStorage nc;
        nc.reserve(cnt);
        nc.x.resize(cnt);
        nc.y.resize(cnt);
        nc.z.resize(cnt);
        for (NodeId i = 0; i < N; ++i) {
            if (rep[i] == i) {
                const NodeId q = newIdx[i];
                nc.x[q] = nodes.x[i];
                nc.y[q] = nodes.y[i];
                nc.z[q] = nodes.z[i];
            }
        }
        nodes = std::move(nc);
    }

    // Rewrite face node indices
    for (auto& v : faces.nodeIndices)
        v = newIdx[v];

    // -------- Detect duplicate boundary faces -----------------------------
    const std::size_t F = faces.size();
    std::vector<std::uint8_t> aliveFace(F, 1);
    std::vector<FaceId> faceRemap(F);
    std::iota(faceRemap.begin(), faceRemap.end(), FaceId{0});

    std::unordered_map<std::vector<NodeId>, FaceId, NodeKeyHash> bndKey;
    bndKey.reserve(F);

    for (FaceId fi = 0; fi < F; ++fi) {
        if (faces.neighbor[fi] != kBoundaryCell)
            continue;
        const std::int32_t s = faces.nodeOffsets[fi];
        const std::int32_t e = faces.nodeOffsets[fi + 1];
        if (e - s < 3)
            continue;
        std::vector<NodeId> key(faces.nodeIndices.begin() + s, faces.nodeIndices.begin() + e);
        std::sort(key.begin(), key.end());
        auto it = bndKey.find(key);
        if (it == bndKey.end()) {
            bndKey.emplace(std::move(key), fi);
        } else {
            const FaceId keep = it->second;
            // Promote `keep` to internal: its neighbor becomes `fi`'s owner.
            faces.neighbor[keep] = faces.owner[fi];
            // Mark `fi` as dead and redirect its remap to `keep`.
            aliveFace[fi] = 0;
            faceRemap[fi] = keep;
            ++stats.deduplicatedFaces;
        }
    }

    if (stats.deduplicatedFaces == 0) {
        out.compute_geometry();
        return stats;
    }

    // -------- Compact face arrays -----------------------------------------
    {
        // First pass: assign new face indices to alive faces.
        std::vector<FaceId> newFi(F, kInvalidId);
        FaceId cnt = 0;
        for (FaceId fi = 0; fi < F; ++fi) {
            if (aliveFace[fi])
                newFi[fi] = cnt++;
        }
        // Resolve faceRemap through aliveFace (a dead face redirects to its
        // live partner; we then translate that live partner to its new id).
        std::vector<FaceId> finalRemap(F);
        for (FaceId fi = 0; fi < F; ++fi) {
            FaceId target = aliveFace[fi] ? fi : faceRemap[fi];
            finalRemap[fi] = newFi[target];
        }

        FaceStorage fc;
        fc.owner.resize(cnt);
        fc.neighbor.resize(cnt);
        fc.areaX.resize(cnt, 0.0);
        fc.areaY.resize(cnt, 0.0);
        fc.areaZ.resize(cnt, 0.0);
        fc.centroidX.resize(cnt, 0.0);
        fc.centroidY.resize(cnt, 0.0);
        fc.centroidZ.resize(cnt, 0.0);
        fc.boundaryZone.resize(cnt, 0);
        // CSR rebuild
        std::vector<std::int32_t> off;
        off.reserve(cnt + 1);
        off.push_back(0);
        fc.nodeIndices.reserve(faces.nodeIndices.size());
        for (FaceId fi = 0; fi < F; ++fi) {
            if (!aliveFace[fi])
                continue;
            const FaceId q = newFi[fi];
            fc.owner[q] = faces.owner[fi];
            fc.neighbor[q] = faces.neighbor[fi];
            fc.boundaryZone[q] = faces.boundaryZone[fi];
            const std::int32_t s = faces.nodeOffsets[fi];
            const std::int32_t e = faces.nodeOffsets[fi + 1];
            for (std::int32_t k = s; k < e; ++k) {
                fc.nodeIndices.push_back(faces.nodeIndices[k]);
            }
            off.push_back(static_cast<std::int32_t>(fc.nodeIndices.size()));
        }
        fc.nodeOffsets.assign(off.begin(), off.end());
        faces = std::move(fc);

        // -------- Rewrite cell face indices --------------------------------
        std::vector<std::int32_t> coff;
        coff.reserve(cells.faceOffsets.size());
        std::vector<FaceId> cidx;
        cidx.reserve(cells.faceIndices.size());
        coff.push_back(0);
        for (std::size_t c = 0; c + 1 < cells.faceOffsets.size(); ++c) {
            const std::int32_t s = cells.faceOffsets[c];
            const std::int32_t e = cells.faceOffsets[c + 1];
            // Use a set-like dedupe so that a cell that now shares the
            // surviving face does not list it twice (shouldn't happen, but
            // protects against weird inputs).
            std::vector<FaceId> seen;
            for (std::int32_t k = s; k < e; ++k) {
                FaceId nf = finalRemap[cells.faceIndices[k]];
                if (nf == kInvalidId)
                    continue;
                if (std::find(seen.begin(), seen.end(), nf) == seen.end()) {
                    seen.push_back(nf);
                }
            }
            for (FaceId nf : seen)
                cidx.push_back(nf);
            coff.push_back(static_cast<std::int32_t>(cidx.size()));
        }
        cells.faceIndices.assign(cidx.begin(), cidx.end());
        cells.faceOffsets.assign(coff.begin(), coff.end());
    }

    out.compute_geometry();
    return stats;
}

// =============================================================================
// Reverse Cuthill-McKee
// =============================================================================
namespace
{

std::size_t cell_bandwidth(const Mesh& m) noexcept
{
    std::size_t bw = 0;
    const auto& f = m.faces();
    for (std::size_t i = 0; i < f.size(); ++i) {
        if (f.neighbor[i] == kBoundaryCell)
            continue;
        const auto a = f.owner[i];
        const auto b = f.neighbor[i];
        const auto d = (a > b) ? (a - b) : (b - a);
        if (d > bw)
            bw = static_cast<std::size_t>(d);
    }
    return bw;
}

} // namespace

RenumberStats renumber_cells_cuthill_mckee(Mesh& m)
{
    RenumberStats stats{};
    stats.bandwidthBefore = cell_bandwidth(m);

    const std::size_t C = m.cells().size();
    if (C == 0) {
        return stats;
    }

    // ---- Build adjacency CSR from interior faces -------------------------
    auto& f = m.faces();
    auto& cells = m.cells();

    std::vector<std::vector<CellId>> adj(C);
    for (std::size_t i = 0; i < f.size(); ++i) {
        if (f.neighbor[i] == kBoundaryCell)
            continue;
        const auto a = f.owner[i];
        const auto b = f.neighbor[i];
        if (a >= C || b >= C)
            continue;
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
    for (auto& v : adj) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    // ---- BFS / Cuthill-McKee ---------------------------------------------
    std::vector<std::uint8_t> seen(C, 0);
    std::vector<CellId> order;
    order.reserve(C);

    // Sort cells by ascending degree for deterministic start choice.
    auto degree = [&](CellId c) { return adj[c].size(); };

    while (order.size() < C) {
        // Find unvisited cell with minimum degree as seed for next component.
        CellId seed = kInvalidId;
        for (CellId c = 0; c < C; ++c) {
            if (seen[c])
                continue;
            if (seed == kInvalidId || degree(c) < degree(seed))
                seed = c;
        }
        if (seed == kInvalidId)
            break;

        std::deque<CellId> q;
        q.push_back(seed);
        seen[seed] = 1;
        while (!q.empty()) {
            CellId c = q.front();
            q.pop_front();
            order.push_back(c);
            // Enqueue unseen neighbours in ascending degree order.
            std::vector<CellId> nbrs;
            nbrs.reserve(adj[c].size());
            for (CellId n : adj[c]) {
                if (!seen[n]) {
                    nbrs.push_back(n);
                    seen[n] = 1;
                }
            }
            std::sort(nbrs.begin(), nbrs.end(), [&](CellId x, CellId y) {
                return degree(x) < degree(y);
            });
            for (CellId n : nbrs)
                q.push_back(n);
        }
    }

    // Reverse for Reverse Cuthill-McKee.
    std::reverse(order.begin(), order.end());

    // ---- Apply permutation ----------------------------------------------
    // newPos[oldCell] = newCellIdx
    std::vector<CellId> newPos(C, kInvalidId);
    for (CellId k = 0; k < order.size(); ++k)
        newPos[order[k]] = k;
    for (CellId c = 0; c < C; ++c) {
        if (newPos[c] == kInvalidId)
            newPos[c] = c; // safety
    }

    // Permute cells.volume + centroids + CSR.
    util::aligned_vector<double> volNew(C, 0.0);
    util::aligned_vector<double> cxNew(C, 0.0), cyNew(C, 0.0), czNew(C, 0.0);
    for (CellId c = 0; c < C; ++c) {
        const CellId q = newPos[c];
        volNew[q] = cells.volume[c];
        cxNew[q] = cells.centroidX[c];
        cyNew[q] = cells.centroidY[c];
        czNew[q] = cells.centroidZ[c];
    }
    cells.volume = std::move(volNew);
    cells.centroidX = std::move(cxNew);
    cells.centroidY = std::move(cyNew);
    cells.centroidZ = std::move(czNew);

    // CSR face list per cell, permuted.
    std::vector<std::int32_t> offNew;
    offNew.reserve(C + 1);
    std::vector<FaceId> idxNew;
    idxNew.reserve(cells.faceIndices.size());
    offNew.push_back(0);
    // We need cell `q`'s face slice, which comes from old cell `order[q]`.
    for (CellId q = 0; q < C; ++q) {
        const CellId c = order[q];
        const std::int32_t s = cells.faceOffsets[c];
        const std::int32_t e = cells.faceOffsets[c + 1];
        for (std::int32_t k = s; k < e; ++k)
            idxNew.push_back(cells.faceIndices[k]);
        offNew.push_back(static_cast<std::int32_t>(idxNew.size()));
    }
    cells.faceIndices.assign(idxNew.begin(), idxNew.end());
    cells.faceOffsets.assign(offNew.begin(), offNew.end());

    // Remap face owner/neighbor through newPos.
    for (auto& v : f.owner)
        v = newPos[v];
    for (auto& v : f.neighbor) {
        if (v != kBoundaryCell)
            v = newPos[v];
    }

    stats.cellsPermuted = order.size();
    stats.bandwidthAfter = cell_bandwidth(m);
    return stats;
}

// =============================================================================
// check_mesh
// =============================================================================
std::string CheckReport::format() const
{
    std::ostringstream os;
    os << "check_mesh: " << (ok ? "OK" : "FAIL") << "\n";
    os << "  nodes            : " << nNodes << "\n";
    os << "  faces            : " << nFaces << " (interior=" << nInteriorFaces
       << ", boundary=" << nBoundaryFaces << ")\n";
    os << "  cells            : " << nCells << "\n";
    os << "  face area        : [" << minFaceArea << ", " << maxFaceArea << "]\n";
    os << "  cell volume      : [" << minCellVolume << ", " << maxCellVolume << "]\n";
    os << "  degenerate faces : " << nDegenerateFaces << "\n";
    os << "  zero-vol cells   : " << nZeroVolumeCells << "\n";
    os << "  negative vol     : " << nNegativeVolume << "\n";
    os << "  orphaned faces   : " << nOrphanedFaces << "\n";
    os << "  max non-ortho deg: " << maxNonOrthoDeg << "\n";
    os << "  max skewness     : " << maxSkewness << "\n";
    os << "  max aspect ratio : " << maxAspectRatio << "\n";
    for (const auto& e : errors)
        os << "  - " << e << "\n";
    return os.str();
}

CheckReport check_mesh(const Mesh& m)
{
    CheckReport r;
    const auto& n = m.nodes();
    const auto& f = m.faces();
    const auto& c = m.cells();

    r.nNodes = n.size();
    r.nFaces = f.size();
    r.nCells = c.size();

    constexpr double kAreaEps = 1.0e-30;
    constexpr double kVolEps = 1.0e-30;

    r.minFaceArea = std::numeric_limits<double>::infinity();
    r.maxFaceArea = -std::numeric_limits<double>::infinity();
    r.minCellVolume = std::numeric_limits<double>::infinity();
    r.maxCellVolume = -std::numeric_limits<double>::infinity();

    for (std::size_t i = 0; i < r.nFaces; ++i) {
        const double a =
            std::sqrt(f.areaX[i] * f.areaX[i] + f.areaY[i] * f.areaY[i] + f.areaZ[i] * f.areaZ[i]);
        if (a < r.minFaceArea)
            r.minFaceArea = a;
        if (a > r.maxFaceArea)
            r.maxFaceArea = a;
        if (a < kAreaEps)
            ++r.nDegenerateFaces;

        if (f.neighbor[i] == kBoundaryCell) {
            ++r.nBoundaryFaces;
        } else {
            ++r.nInteriorFaces;
        }
        if (f.owner[i] >= r.nCells
            || (f.neighbor[i] != kBoundaryCell && f.neighbor[i] >= r.nCells)) {
            ++r.nOrphanedFaces;
        }
    }

    for (std::size_t i = 0; i < r.nCells; ++i) {
        const double v = c.volume[i];
        if (v < r.minCellVolume)
            r.minCellVolume = v;
        if (v > r.maxCellVolume)
            r.maxCellVolume = v;
        if (std::abs(v) < kVolEps)
            ++r.nZeroVolumeCells;
        if (v < 0.0)
            ++r.nNegativeVolume;
    }

    if (r.nFaces == 0)
        r.minFaceArea = r.maxFaceArea = 0.0;
    if (r.nCells == 0)
        r.minCellVolume = r.maxCellVolume = 0.0;

    // Layer the MeshQuality summary on top.
    const auto q = MeshQuality::evaluate(m);
    r.maxNonOrthoDeg = q.maxNonOrtho;
    r.maxSkewness = q.maxSkewness;
    r.maxAspectRatio = q.maxAspect;

    if (r.nDegenerateFaces > 0)
        r.errors.push_back("degenerate faces detected");
    if (r.nZeroVolumeCells > 0)
        r.errors.push_back("zero-volume cells detected");
    if (r.nNegativeVolume > 0)
        r.errors.push_back("negative-volume cells detected");
    if (r.nOrphanedFaces > 0)
        r.errors.push_back("faces reference out-of-range cells");
    r.ok = r.errors.empty();
    return r;
}

// =============================================================================
// Sub-mesh extraction & refinement  (Pass 11b)
// =============================================================================
namespace
{

// Re-emit the (kept) cells through ConnectivityBuilder so the result has
// fully consistent CSR + geometry; this is also how `refine_hex` commits.
void rebuild_from_cell_descriptors(const NodeStorage& ns,
                                   const std::vector<CellDescriptor>& cells,
                                   Mesh& out)
{
    Mesh fresh;
    ConnectivityBuilder::build(fresh, ns, cells);
    out.nodes() = std::move(fresh.nodes());
    out.faces() = std::move(fresh.faces());
    out.cells() = std::move(fresh.cells());
    out.compute_geometry();
}

} // namespace

SplitStats split_mesh(const Mesh& src, const std::vector<std::uint8_t>& cellMask, Mesh& out)
{
    SplitStats stats{};
    const auto& sn = src.nodes();
    const auto& sf = src.faces();
    const auto& sc = src.cells();
    const std::size_t nCells = sc.size();

    if (cellMask.size() != nCells) {
        // Caller error - treat as "select nothing".
        out = Mesh{};
        return stats;
    }

    // Pass 1: count and remap kept nodes.
    std::vector<std::uint8_t> nodeKept(sn.size(), 0);
    for (std::size_t c = 0; c < nCells; ++c) {
        if (!cellMask[c])
            continue;
        ++stats.selectedCells;
        const std::int32_t fb = sc.faceOffsets[c];
        const std::int32_t fe = sc.faceOffsets[c + 1];
        for (std::int32_t fi = fb; fi < fe; ++fi) {
            const FaceId f = sc.faceIndices[fi];
            const std::int32_t nb = sf.nodeOffsets[f];
            const std::int32_t ne = sf.nodeOffsets[f + 1];
            for (std::int32_t ni = nb; ni < ne; ++ni)
                nodeKept[sf.nodeIndices[ni]] = 1;
        }
    }

    NodeStorage outNodes;
    std::vector<NodeId> nodeRemap(sn.size(), 0);
    for (std::size_t i = 0; i < sn.size(); ++i) {
        if (!nodeKept[i])
            continue;
        nodeRemap[i] = static_cast<NodeId>(outNodes.x.size());
        outNodes.x.push_back(sn.x[i]);
        outNodes.y.push_back(sn.y[i]);
        outNodes.z.push_back(sn.z[i]);
    }

    // Pass 2: emit selected cells, remapping node IDs into the new numbering.
    std::vector<CellDescriptor> outCells;
    outCells.reserve(stats.selectedCells);
    for (std::size_t c = 0; c < nCells; ++c) {
        if (!cellMask[c])
            continue;
        CellDescriptor cd;
        const std::int32_t fb = sc.faceOffsets[c];
        const std::int32_t fe = sc.faceOffsets[c + 1];
        cd.faces.reserve(static_cast<std::size_t>(fe - fb));
        for (std::int32_t fi = fb; fi < fe; ++fi) {
            const FaceId f = sc.faceIndices[fi];
            const std::int32_t nb = sf.nodeOffsets[f];
            const std::int32_t ne = sf.nodeOffsets[f + 1];
            std::vector<NodeId> face;
            face.reserve(static_cast<std::size_t>(ne - nb));
            for (std::int32_t ni = nb; ni < ne; ++ni)
                face.push_back(nodeRemap[sf.nodeIndices[ni]]);
            // Preserve original face orientation as authored on this cell:
            // if the source cell was the face's neighbour (not owner), the
            // canonical node ordering points the other way - reverse it so
            // the new face is owned by this cell with the correct outward
            // normal.
            if (sf.neighbor[f] != kBoundaryCell && sf.neighbor[f] == c)
                std::reverse(face.begin(), face.end());
            cd.faces.push_back(std::move(face));
        }
        outCells.push_back(std::move(cd));
    }

    rebuild_from_cell_descriptors(outNodes, outCells, out);
    stats.outputCells = out.cells().size();
    stats.outputFaces = out.faces().size();
    stats.outputNodes = out.nodes().size();
    return stats;
}

// -----------------------------------------------------------------------------
// refine_hex: topological hex detection + 1-to-8 octant split
// -----------------------------------------------------------------------------
namespace
{

struct HexCorners
{
    // Indexed [k][j][i], i/j/k in {0,1}.
    std::array<NodeId, 8> n{};
    NodeId& at(int i, int j, int k) { return n[k * 4 + j * 2 + i]; }
    NodeId at(int i, int j, int k) const { return n[k * 4 + j * 2 + i]; }
};

inline std::pair<NodeId, NodeId> edge_key(NodeId a, NodeId b)
{
    return (a < b) ? std::pair{a, b} : std::pair{b, a};
}

struct PairHash
{
    std::size_t operator()(const std::pair<NodeId, NodeId>& p) const noexcept
    {
        return std::hash<NodeId>{}(p.first) * 1469598103934665603ull
               ^ std::hash<NodeId>{}(p.second);
    }
};

struct VecHash
{
    std::size_t operator()(const std::vector<NodeId>& v) const noexcept
    {
        std::size_t h = 1469598103934665603ull;
        for (NodeId n : v)
            h ^= n + 0x9e3779b97f4a7c15ull + (h << 12) + (h >> 4);
        return h;
    }
};

/// Detect a topological hex.  Reads the 6 face vertex lists from `cellFaces`
/// (each must be a CCW 4-cycle from `mesh.faces()`); on success populates
/// `H` with the 8 corners in (i,j,k) order.  Returns false if any face is
/// not a quad, the unique node count is not 8, or any node does not have
/// degree-3 in the cell.
bool detect_hex(const std::vector<std::vector<NodeId>>& cellFaces, HexCorners& H)
{
    if (cellFaces.size() != 6)
        return false;
    for (const auto& f : cellFaces)
        if (f.size() != 4)
            return false;

    // Unique node set + per-node face membership (face indices 0..5).
    std::unordered_map<NodeId, std::array<int, 3>> nodeFaces; // value = face indices
    std::unordered_map<NodeId, int> nodeCount;
    for (int fi = 0; fi < 6; ++fi) {
        for (NodeId n : cellFaces[fi]) {
            int& c = nodeCount[n];
            if (c < 3)
                nodeFaces[n][c] = fi;
            ++c;
        }
    }
    if (nodeCount.size() != 8)
        return false;
    for (auto& [n, c] : nodeCount)
        if (c != 3)
            return false;

    // Pick N0 = lowest-id node.
    NodeId N0 = std::numeric_limits<NodeId>::max();
    for (auto& [n, _] : nodeCount)
        N0 = std::min(N0, n);

    const auto& F = nodeFaces[N0];
    const int f0 = F[0], fA_cand1 = F[1], fA_cand2 = F[2];

    // N0's two cycle-neighbours in f0.
    const auto& v0 = cellFaces[f0];
    int idx0 = -1;
    for (int i = 0; i < 4; ++i)
        if (v0[i] == N0) {
            idx0 = i;
            break;
        }
    if (idx0 < 0)
        return false;
    const NodeId Ea = v0[(idx0 + 3) & 3];     // predecessor in cycle
    const NodeId Eb = v0[(idx0 + 1) & 3];     // successor in cycle
    const NodeId opp_f0 = v0[(idx0 + 2) & 3]; // diagonal of f0

    // Determine fA = the OTHER face (among N0's 3) that also contains Ea,
    // and fB = the OTHER face that also contains Eb.
    auto contains = [&](int fi, NodeId n) {
        for (NodeId v : cellFaces[fi])
            if (v == n)
                return true;
        return false;
    };
    const int fA = contains(fA_cand1, Ea) ? fA_cand1 : fA_cand2;
    const int fB = (fA == fA_cand1) ? fA_cand2 : fA_cand1;
    if (!contains(fA, Ea) || !contains(fB, Eb))
        return false;

    // E_z = node in fA ∩ fB, !=N0. (The corner along the third axis.)
    NodeId Ec = std::numeric_limits<NodeId>::max();
    for (NodeId v : cellFaces[fA])
        if (v != N0 && contains(fB, v)) {
            Ec = v;
            break;
        }
    if (Ec == std::numeric_limits<NodeId>::max())
        return false;

    // Diagonal of fA = corner (1,0,1).  Found by index-of-N0 in fA + 2.
    const auto& vA = cellFaces[fA];
    int idxA = -1;
    for (int i = 0; i < 4; ++i)
        if (vA[i] == N0) {
            idxA = i;
            break;
        }
    if (idxA < 0)
        return false;
    const NodeId opp_fA = vA[(idxA + 2) & 3];

    const auto& vB = cellFaces[fB];
    int idxB = -1;
    for (int i = 0; i < 4; ++i)
        if (vB[i] == N0) {
            idxB = i;
            break;
        }
    if (idxB < 0)
        return false;
    const NodeId opp_fB = vB[(idxB + 2) & 3];

    // (1,1,1) = unique remaining node.
    std::array<NodeId, 7> known{N0, Ea, Eb, Ec, opp_f0, opp_fA, opp_fB};
    NodeId far = std::numeric_limits<NodeId>::max();
    for (auto& [n, _] : nodeCount) {
        bool seen = false;
        for (NodeId k : known)
            if (k == n) {
                seen = true;
                break;
            }
        if (!seen) {
            far = n;
            break;
        }
    }
    if (far == std::numeric_limits<NodeId>::max())
        return false;

    // Label corners. The mapping is:
    //   (0,0,0) = N0       (1,0,0) = Ea      (0,1,0) = Eb      (0,0,1) = Ec
    //   (1,1,0) = opp_f0   (1,0,1) = opp_fA  (0,1,1) = opp_fB  (1,1,1) = far
    H.at(0, 0, 0) = N0;
    H.at(1, 0, 0) = Ea;
    H.at(0, 1, 0) = Eb;
    H.at(0, 0, 1) = Ec;
    H.at(1, 1, 0) = opp_f0;
    H.at(1, 0, 1) = opp_fA;
    H.at(0, 1, 1) = opp_fB;
    H.at(1, 1, 1) = far;
    return true;
}

/// Pull cell `c`'s face vertex lists from the mesh, oriented as authored
/// (owner cells see the face CCW; for cells that are the face's neighbour
/// we reverse the loop so all 6 faces have outward normals from this cell).
std::vector<std::vector<NodeId>> gather_cell_faces(const Mesh& m, CellId c)
{
    const auto& sf = m.faces();
    const auto& sc = m.cells();
    const std::int32_t fb = sc.faceOffsets[c];
    const std::int32_t fe = sc.faceOffsets[c + 1];
    std::vector<std::vector<NodeId>> out;
    out.reserve(static_cast<std::size_t>(fe - fb));
    for (std::int32_t fi = fb; fi < fe; ++fi) {
        const FaceId f = sc.faceIndices[fi];
        const std::int32_t nb = sf.nodeOffsets[f];
        const std::int32_t ne = sf.nodeOffsets[f + 1];
        std::vector<NodeId> face;
        face.reserve(static_cast<std::size_t>(ne - nb));
        for (std::int32_t ni = nb; ni < ne; ++ni)
            face.push_back(sf.nodeIndices[ni]);
        if (sf.neighbor[f] != kBoundaryCell && sf.neighbor[f] == c)
            std::reverse(face.begin(), face.end());
        out.push_back(std::move(face));
    }
    return out;
}

} // namespace

RefineStats refine_hex(Mesh& m)
{
    RefineStats stats{};
    const std::size_t nCells = m.cells().size();
    stats.cellsBefore = nCells;
    stats.nodesBefore = m.nodes().size();

    // Detect every hex up-front so we can reject the whole pass if any cell
    // fails - but per the header contract, we instead PASS THROUGH rejected
    // cells (so partial-hex meshes get partially refined).
    std::vector<HexCorners> hex(nCells);
    std::vector<std::uint8_t> isHex(nCells, 0);
    for (std::size_t c = 0; c < nCells; ++c) {
        auto cellFaces = gather_cell_faces(m, c);
        if (detect_hex(cellFaces, hex[c])) {
            isHex[c] = 1;
            ++stats.cellsRefined;
        } else {
            ++stats.rejectedCells;
        }
    }

    // ---- Build the augmented node set --------------------------------------
    NodeStorage ns = m.nodes(); // start with all original nodes

    std::unordered_map<std::pair<NodeId, NodeId>, NodeId, PairHash> edgeMid;
    std::unordered_map<std::vector<NodeId>, NodeId, VecHash> faceMid;
    std::vector<NodeId> cellMid(nCells, std::numeric_limits<NodeId>::max());

    auto get_edge_mid = [&](NodeId a, NodeId b) -> NodeId {
        auto key = edge_key(a, b);
        auto it = edgeMid.find(key);
        if (it != edgeMid.end())
            return it->second;
        const NodeId id = static_cast<NodeId>(ns.x.size());
        ns.x.push_back(0.5 * (m.nodes().x[a] + m.nodes().x[b]));
        ns.y.push_back(0.5 * (m.nodes().y[a] + m.nodes().y[b]));
        ns.z.push_back(0.5 * (m.nodes().z[a] + m.nodes().z[b]));
        edgeMid.emplace(key, id);
        return id;
    };
    auto get_face_mid = [&](std::array<NodeId, 4> f) -> NodeId {
        std::vector<NodeId> key(f.begin(), f.end());
        std::sort(key.begin(), key.end());
        auto it = faceMid.find(key);
        if (it != faceMid.end())
            return it->second;
        const NodeId id = static_cast<NodeId>(ns.x.size());
        double cx = 0, cy = 0, cz = 0;
        for (NodeId v : f) {
            cx += m.nodes().x[v];
            cy += m.nodes().y[v];
            cz += m.nodes().z[v];
        }
        ns.x.push_back(0.25 * cx);
        ns.y.push_back(0.25 * cy);
        ns.z.push_back(0.25 * cz);
        faceMid.emplace(std::move(key), id);
        return id;
    };
    auto get_cell_mid = [&](CellId c) -> NodeId {
        if (cellMid[c] != std::numeric_limits<NodeId>::max())
            return cellMid[c];
        const HexCorners& H = hex[c];
        double cx = 0, cy = 0, cz = 0;
        for (NodeId v : H.n) {
            cx += m.nodes().x[v];
            cy += m.nodes().y[v];
            cz += m.nodes().z[v];
        }
        const NodeId id = static_cast<NodeId>(ns.x.size());
        ns.x.push_back(0.125 * cx);
        ns.y.push_back(0.125 * cy);
        ns.z.push_back(0.125 * cz);
        cellMid[c] = id;
        return id;
    };

    // ---- Emit refined cells -------------------------------------------------
    std::vector<CellDescriptor> outCells;
    outCells.reserve(stats.cellsRefined * 8 + stats.rejectedCells);

    // Helper: 27-point table for cell c.  Index by (i,j,k), i/j/k in {0,1,2}
    // where 0/2 are corners, 1 is the midpoint along that axis.
    auto build_p27 = [&](CellId c, std::array<NodeId, 27>& P) {
        const HexCorners& H = hex[c];
        auto C = [&](int i, int j, int k) { return H.at(i, j, k); };
        auto put = [&](int i, int j, int k, NodeId id) { P[k * 9 + j * 3 + i] = id; };
        // 8 corners
        put(0, 0, 0, C(0, 0, 0));
        put(2, 0, 0, C(1, 0, 0));
        put(0, 2, 0, C(0, 1, 0));
        put(2, 2, 0, C(1, 1, 0));
        put(0, 0, 2, C(0, 0, 1));
        put(2, 0, 2, C(1, 0, 1));
        put(0, 2, 2, C(0, 1, 1));
        put(2, 2, 2, C(1, 1, 1));
        // 12 edge midpoints
        put(1, 0, 0, get_edge_mid(C(0, 0, 0), C(1, 0, 0)));
        put(1, 2, 0, get_edge_mid(C(0, 1, 0), C(1, 1, 0)));
        put(1, 0, 2, get_edge_mid(C(0, 0, 1), C(1, 0, 1)));
        put(1, 2, 2, get_edge_mid(C(0, 1, 1), C(1, 1, 1)));
        put(0, 1, 0, get_edge_mid(C(0, 0, 0), C(0, 1, 0)));
        put(2, 1, 0, get_edge_mid(C(1, 0, 0), C(1, 1, 0)));
        put(0, 1, 2, get_edge_mid(C(0, 0, 1), C(0, 1, 1)));
        put(2, 1, 2, get_edge_mid(C(1, 0, 1), C(1, 1, 1)));
        put(0, 0, 1, get_edge_mid(C(0, 0, 0), C(0, 0, 1)));
        put(2, 0, 1, get_edge_mid(C(1, 0, 0), C(1, 0, 1)));
        put(0, 2, 1, get_edge_mid(C(0, 1, 0), C(0, 1, 1)));
        put(2, 2, 1, get_edge_mid(C(1, 1, 0), C(1, 1, 1)));
        // 6 face centroids
        put(1, 1, 0, get_face_mid({C(0, 0, 0), C(1, 0, 0), C(1, 1, 0), C(0, 1, 0)}));
        put(1, 1, 2, get_face_mid({C(0, 0, 1), C(1, 0, 1), C(1, 1, 1), C(0, 1, 1)}));
        put(1, 0, 1, get_face_mid({C(0, 0, 0), C(1, 0, 0), C(1, 0, 1), C(0, 0, 1)}));
        put(1, 2, 1, get_face_mid({C(0, 1, 0), C(1, 1, 0), C(1, 1, 1), C(0, 1, 1)}));
        put(0, 1, 1, get_face_mid({C(0, 0, 0), C(0, 1, 0), C(0, 1, 1), C(0, 0, 1)}));
        put(2, 1, 1, get_face_mid({C(1, 0, 0), C(1, 1, 0), C(1, 1, 1), C(1, 0, 1)}));
        // 1 cell centroid
        put(1, 1, 1, get_cell_mid(c));
    };

    for (std::size_t c = 0; c < nCells; ++c) {
        if (!isHex[c]) {
            outCells.push_back(CellDescriptor{gather_cell_faces(m, c), 0});
            continue;
        }
        std::array<NodeId, 27> P{};
        build_p27(static_cast<CellId>(c), P);
        auto N = [&](int i, int j, int k) { return P[k * 9 + j * 3 + i]; };

        // 8 sub-hexes at octants (oi, oj, ok) in {0,1}.
        for (int ok = 0; ok < 2; ++ok)
            for (int oj = 0; oj < 2; ++oj)
                for (int oi = 0; oi < 2; ++oi) {
                    const NodeId c000 = N(oi, oj, ok);
                    const NodeId c100 = N(oi + 1, oj, ok);
                    const NodeId c010 = N(oi, oj + 1, ok);
                    const NodeId c110 = N(oi + 1, oj + 1, ok);
                    const NodeId c001 = N(oi, oj, ok + 1);
                    const NodeId c101 = N(oi + 1, oj, ok + 1);
                    const NodeId c011 = N(oi, oj + 1, ok + 1);
                    const NodeId c111 = N(oi + 1, oj + 1, ok + 1);
                    CellDescriptor cd;
                    // VTK-hex convention used elsewhere (CartesianMesher).
                    cd.faces = {
                        {c000, c010, c110, c100}, // -Z
                        {c001, c101, c111, c011}, // +Z
                        {c000, c100, c101, c001}, // -Y
                        {c010, c011, c111, c110}, // +Y
                        {c000, c001, c011, c010}, // -X
                        {c100, c110, c111, c101}, // +X
                    };
                    outCells.push_back(std::move(cd));
                }
    }

    rebuild_from_cell_descriptors(ns, outCells, m);
    stats.cellsAfter = m.cells().size();
    stats.nodesAfter = m.nodes().size();
    return stats;
}

// =============================================================================
// extract_subdomain (Pass 20): owned + optional one-deep halo via split_mesh
// =============================================================================
SubdomainStats extract_subdomain(const Mesh& src,
                                 const std::vector<std::int32_t>& cellRank,
                                 std::int32_t rank,
                                 bool includeGhostLayer,
                                 Mesh& out)
{
    const auto& C = src.cells();
    const auto& F = src.faces();
    const std::size_t nC = C.size();
    if (cellRank.size() != nC) {
        throw std::invalid_argument(
            "extract_subdomain: cellRank.size() must equal mesh cell count");
    }

    SubdomainStats stats;
    std::vector<std::uint8_t> mask(nC, 0);

    for (std::size_t c = 0; c < nC; ++c) {
        if (cellRank[c] == rank) {
            mask[c] = 1;
            stats.ownedCells++;
        }
    }

    // Pass 1: count interface faces (owned <-> non-owned interior faces).
    const std::size_t nF = F.size();
    for (std::size_t f = 0; f < nF; ++f) {
        const auto o = F.owner[f];
        const auto nb = F.neighbor[f];
        if (nb == kBoundaryCell)
            continue;
        const bool oOwned = (cellRank[o] == rank);
        const bool nbOwned = (cellRank[nb] == rank);
        if (oOwned != nbOwned)
            stats.interfaceFaces++;
    }

    // Pass 2: extend mask with one-deep halo of non-owned cells.
    if (includeGhostLayer) {
        for (std::size_t f = 0; f < nF; ++f) {
            const auto o = F.owner[f];
            const auto nb = F.neighbor[f];
            if (nb == kBoundaryCell)
                continue;
            const bool oOwned = (cellRank[o] == rank);
            const bool nbOwned = (cellRank[nb] == rank);
            if (oOwned && !nbOwned && !mask[nb]) {
                mask[nb] = 1;
                stats.ghostCells++;
            } else if (nbOwned && !oOwned && !mask[o]) {
                mask[o] = 1;
                stats.ghostCells++;
            }
        }
    }

    stats.split = split_mesh(src, mask, out);
    return stats;
}

} // namespace simall::meshing::ops
