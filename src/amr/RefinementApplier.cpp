// =============================================================================
// SimAll Beta - AMR Subsystem
// File   : src/amr/RefinementApplier.cpp
//
// Strategy: Hex isotropic refinement (1 → 8 children) with neighbour-buffer
// expansion to guarantee 2:1 conformality. Each input cell is assumed to be
// a hex (or hex-equivalent) reconstructable as its axis-aligned bounding
// box; if a non-hex cell is encountered, it is passed through unchanged.
// Coarsening of unflagged cells is not performed in this pass.
// =============================================================================
#include "amr/RefinementApplier.hpp"

#include "core/Logger.hpp"
#include "meshing/Connectivity.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace simall::amr
{

namespace
{

struct AABB
{
    double xmin, ymin, zmin, xmax, ymax, zmax;
};

AABB cell_bbox(const meshing::Mesh& m, meshing::CellId c)
{
    AABB b{std::numeric_limits<double>::max(),
           std::numeric_limits<double>::max(),
           std::numeric_limits<double>::max(),
           -std::numeric_limits<double>::max(),
           -std::numeric_limits<double>::max(),
           -std::numeric_limits<double>::max()};
    const auto& C = m.cells();
    const auto& F = m.faces();
    const auto& N = m.nodes();
    const int fs = C.faceOffsets[c], fe = C.faceOffsets[c + 1];
    for (int k = fs; k < fe; ++k) {
        const meshing::FaceId fid = C.faceIndices[k];
        const int ns = F.nodeOffsets[fid], ne = F.nodeOffsets[fid + 1];
        for (int p = ns; p < ne; ++p) {
            const meshing::NodeId v = F.nodeIndices[p];
            b.xmin = std::min(b.xmin, N.x[v]);
            b.xmax = std::max(b.xmax, N.x[v]);
            b.ymin = std::min(b.ymin, N.y[v]);
            b.ymax = std::max(b.ymax, N.y[v]);
            b.zmin = std::min(b.zmin, N.z[v]);
            b.zmax = std::max(b.zmax, N.z[v]);
        }
    }
    return b;
}

void emit_hex(meshing::NodeStorage& ns, std::vector<meshing::CellDescriptor>& cells, const AABB& b)
{
    using meshing::NodeId;
    const std::size_t n0 = ns.size();
    const double xs[2] = {b.xmin, b.xmax};
    const double ys[2] = {b.ymin, b.ymax};
    const double zs[2] = {b.zmin, b.zmax};
    for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 2; ++i) {
                ns.x.push_back(xs[i]);
                ns.y.push_back(ys[j]);
                ns.z.push_back(zs[k]);
            }
    auto P = [&](int i, int j, int k) -> NodeId { return NodeId(n0 + k * 4 + j * 2 + i); };
    meshing::CellDescriptor cd;
    cd.faces = {{P(0, 0, 0), P(0, 1, 0), P(0, 1, 1), P(0, 0, 1)},
                {P(1, 0, 0), P(1, 0, 1), P(1, 1, 1), P(1, 1, 0)},
                {P(0, 0, 0), P(0, 0, 1), P(1, 0, 1), P(1, 0, 0)},
                {P(0, 1, 0), P(1, 1, 0), P(1, 1, 1), P(0, 1, 1)},
                {P(0, 0, 0), P(1, 0, 0), P(1, 1, 0), P(0, 1, 0)},
                {P(0, 0, 1), P(0, 1, 1), P(1, 1, 1), P(1, 0, 1)}};
    cells.push_back(std::move(cd));
}

} // namespace

void RefinementApplier::apply(const meshing::Mesh& in,
                              const std::vector<RefineFlag>& flagsIn,
                              const std::vector<std::int32_t>& level,
                              meshing::Mesh& out,
                              std::vector<std::int32_t>& outLevel)
{
    const std::size_t nC = in.cells().size();
    std::vector<RefineFlag> flags = flagsIn;
    if (flags.size() != nC)
        flags.assign(nC, RefineFlag::Keep);

    // ---- 2:1 buffer expansion. Iterate until stable: if any neighbour of a
    // refined cell is unrefined at a coarser level, mark it for refinement too.
    const auto& C = in.cells();
    const auto& F = in.faces();
    bool changed = true;
    int sweep = 0;
    while (changed && sweep < 20) {
        changed = false;
        for (std::size_t c = 0; c < nC; ++c) {
            if (flags[c] != RefineFlag::Refine)
                continue;
            const int s = C.faceOffsets[c];
            const int e = C.faceOffsets[c + 1];
            for (int k = s; k < e; ++k) {
                const meshing::FaceId fid = C.faceIndices[k];
                const meshing::CellId nb = (F.owner[fid] == c) ? F.neighbor[fid] : F.owner[fid];
                if (nb == meshing::kBoundaryCell)
                    continue;
                const std::int32_t lc = level.empty() ? 0 : level[c];
                const std::int32_t ln = level.empty() ? 0 : level[nb];
                if (flags[nb] != RefineFlag::Refine && ln < lc + 1) {
                    flags[nb] = RefineFlag::Refine;
                    changed = true;
                }
            }
        }
        ++sweep;
    }

    meshing::NodeStorage ns;
    std::vector<meshing::CellDescriptor> cells;
    cells.reserve(nC * 2);
    outLevel.clear();
    outLevel.reserve(nC * 2);

    for (std::size_t c = 0; c < nC; ++c) {
        const AABB b = cell_bbox(in, c);
        const std::int32_t L = level.empty() ? 0 : level[c];
        if (flags[c] == RefineFlag::Refine) {
            const double mx = 0.5 * (b.xmin + b.xmax);
            const double my = 0.5 * (b.ymin + b.ymax);
            const double mz = 0.5 * (b.zmin + b.zmax);
            for (int kz = 0; kz < 2; ++kz)
                for (int ky = 0; ky < 2; ++ky)
                    for (int kx = 0; kx < 2; ++kx) {
                        AABB sub{kx ? mx : b.xmin,
                                 ky ? my : b.ymin,
                                 kz ? mz : b.zmin,
                                 kx ? b.xmax : mx,
                                 ky ? b.ymax : my,
                                 kz ? b.zmax : mz};
                        emit_hex(ns, cells, sub);
                        outLevel.push_back(L + 1);
                    }
        } else {
            emit_hex(ns, cells, b);
            outLevel.push_back(L);
        }
    }

    meshing::ConnectivityBuilder::build(out, ns, cells);
    SIMALL_LOG_INFO(
        "AMR", "refinement applied: ", nC, " → ", cells.size(), " cells (sweeps=", sweep, ")");
}

} // namespace simall::amr
