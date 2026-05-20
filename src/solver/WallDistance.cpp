// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/WallDistance.cpp
//
// References:
//   - Spalding (1994): "Calculation of turbulent heat transfer in cluttered
//     spaces", Proc. 10th Int. Heat Transfer Conf.
//   - Tucker (2003): "Differential equation-based wall distance computation
//     for DES and RANS", J. Comput. Phys. 190, 229-248.
//   - Wukie & Orkwis (2018): "Robust wall-distance via Poisson smoothing".
// =============================================================================
#include "solver/WallDistance.hpp"

#include "core/Logger.hpp"
#include "solver/CSRMatrix.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

namespace simall::solver
{

namespace
{

// ----- closest point on a triangle (Ericson, RTCD §5.1.5) -----------------
util::Vec3d closest_on_triangle(const util::Vec3d& p,
                                const util::Vec3d& a,
                                const util::Vec3d& b,
                                const util::Vec3d& c)
{
    const double abx = b.x - a.x, aby = b.y - a.y, abz = b.z - a.z;
    const double acx = c.x - a.x, acy = c.y - a.y, acz = c.z - a.z;
    const double apx = p.x - a.x, apy = p.y - a.y, apz = p.z - a.z;
    const double d1 = abx * apx + aby * apy + abz * apz;
    const double d2 = acx * apx + acy * apy + acz * apz;
    if (d1 <= 0 && d2 <= 0)
        return a;
    const double bpx = p.x - b.x, bpy = p.y - b.y, bpz = p.z - b.z;
    const double d3 = abx * bpx + aby * bpy + abz * bpz;
    const double d4 = acx * bpx + acy * bpy + acz * bpz;
    if (d3 >= 0 && d4 <= d3)
        return b;
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        const double v = d1 / (d1 - d3);
        return {a.x + v * abx, a.y + v * aby, a.z + v * abz};
    }
    const double cpx = p.x - c.x, cpy = p.y - c.y, cpz = p.z - c.z;
    const double d5 = abx * cpx + aby * cpy + abz * cpz;
    const double d6 = acx * cpx + acy * cpy + acz * cpz;
    if (d6 >= 0 && d5 <= d6)
        return c;
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        const double w = d2 / (d2 - d6);
        return {a.x + w * acx, a.y + w * acy, a.z + w * acz};
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return {b.x + w * (c.x - b.x), b.y + w * (c.y - b.y), b.z + w * (c.z - b.z)};
    }
    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom, w = vc * denom;
    return {a.x + abx * v + acx * w, a.y + aby * v + acy * w, a.z + abz * v + acz * w};
}

bool is_wall_bc(BCType t)
{
    return t == BCType::Wall || t == BCType::NoSlipWall || t == BCType::MovingWall
           || t == BCType::SlipWall;
}

// ============================================================ Exact
class ExactNearest final : public IWallDistance
{
public:
    void compute(const meshing::Mesh& mesh,
                 const std::vector<BoundarySpec>& bcs,
                 ScalarField& out) override
    {
        const auto& C = mesh.cells();
        const auto& F = mesh.faces();
        const auto& N = mesh.nodes();
        const std::size_t nC = C.size();
        out.assign(nC, std::numeric_limits<double>::max());

        // Collect wall faces.
        std::unordered_set<meshing::ZoneId> wallZones;
        for (const auto& b : bcs)
            if (is_wall_bc(b.type))
                wallZones.insert(b.zone);

        struct Tri
        {
            util::Vec3d a, b, c;
        };
        std::vector<Tri> wallTris;
        for (std::size_t f = 0; f < F.size(); ++f) {
            if (!wallZones.count(F.boundaryZone[f]))
                continue;
            const int beg = F.nodeOffsets[f], end = F.nodeOffsets[f + 1];
            // Fan-triangulate the polygonal face from node 0.
            const meshing::NodeId n0 = F.nodeIndices[beg];
            for (int k = beg + 1; k + 1 < end; ++k) {
                wallTris.push_back(
                    {{N.x[n0], N.y[n0], N.z[n0]},
                     {N.x[F.nodeIndices[k]], N.y[F.nodeIndices[k]], N.z[F.nodeIndices[k]]},
                     {N.x[F.nodeIndices[k + 1]],
                      N.y[F.nodeIndices[k + 1]],
                      N.z[F.nodeIndices[k + 1]]}});
            }
        }
        if (wallTris.empty()) {
            std::fill(out.begin(), out.end(), 0.0);
            return;
        }

        for (std::size_t c = 0; c < nC; ++c) {
            const util::Vec3d p{C.centroidX[c], C.centroidY[c], C.centroidZ[c]};
            double best = std::numeric_limits<double>::max();
            for (const Tri& t : wallTris) {
                const util::Vec3d q = closest_on_triangle(p, t.a, t.b, t.c);
                const double dx = p.x - q.x, dy = p.y - q.y, dz = p.z - q.z;
                const double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < best)
                    best = d2;
            }
            out[c] = std::sqrt(best);
        }
        SIMALL_LOG_INFO("Solver",
                        "Wall distance (exact): ",
                        wallTris.size(),
                        " wall triangles vs ",
                        nC,
                        " cells");
    }
};

// ============================================================ Poisson
class PoissonWallDistance final : public IWallDistance
{
public:
    explicit PoissonWallDistance(ILinearSolver& s) : sol_(s) {}

    void compute(const meshing::Mesh& mesh,
                 const std::vector<BoundarySpec>& bcs,
                 ScalarField& out) override
    {
        const auto& C = mesh.cells();
        const auto& F = mesh.faces();
        const std::size_t nC = C.size();

        std::unordered_set<meshing::ZoneId> wallZones;
        for (const auto& b : bcs)
            if (is_wall_bc(b.type))
                wallZones.insert(b.zone);

        // Build CSR sparsity (self + face-neighbours).
        CSRMatrix A;
        A.rowPtr.assign(nC + 1, 0);
        std::vector<std::vector<int>> nbr(nC);
        for (std::size_t c = 0; c < nC; ++c) {
            std::vector<int> ns{static_cast<int>(c)};
            const int beg = C.faceOffsets[c], end = C.faceOffsets[c + 1];
            for (int k = beg; k < end; ++k) {
                const meshing::FaceId fid = C.faceIndices[k];
                const meshing::CellId oth = (F.owner[fid] == c) ? F.neighbor[fid] : F.owner[fid];
                if (oth != meshing::kBoundaryCell)
                    ns.push_back(static_cast<int>(oth));
            }
            std::sort(ns.begin(), ns.end());
            ns.erase(std::unique(ns.begin(), ns.end()), ns.end());
            nbr[c] = std::move(ns);
            A.rowPtr[c + 1] = A.rowPtr[c] + static_cast<int>(nbr[c].size());
        }
        A.colIdx.reserve(A.rowPtr[nC]);
        for (const auto& v : nbr)
            for (int n : v)
                A.colIdx.push_back(n);
        A.values.assign(A.colIdx.size(), 0.0);

        auto find = [&](int row, int col) -> int {
            for (int k = A.rowPtr[row]; k < A.rowPtr[row + 1]; ++k)
                if (A.colIdx[k] == col)
                    return k;
            return -1;
        };
        auto add = [&](int row, int col, double v) {
            const int k = find(row, col);
            if (k >= 0)
                A.values[k] += v;
        };

        util::aligned_vector<double> b(nC, 1.0); // RHS = +1·V_c
        for (std::size_t c = 0; c < nC; ++c)
            b[c] *= C.volume[c];

        for (std::size_t f = 0; f < F.size(); ++f) {
            const meshing::CellId o = F.owner[f];
            const meshing::CellId n = F.neighbor[f];
            const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
            const double Amag2 = Ax * Ax + Ay * Ay + Az * Az;
            if (n != meshing::kBoundaryCell) {
                const double dx = C.centroidX[n] - C.centroidX[o];
                const double dy = C.centroidY[n] - C.centroidY[o];
                const double dz = C.centroidZ[n] - C.centroidZ[o];
                const double dn = std::abs(dx * Ax + dy * Ay + dz * Az);
                const double D = Amag2 / std::max(dn, 1e-30);
                add(o, o, D);
                add(o, n, -D);
                add(n, n, D);
                add(n, o, -D);
            } else if (wallZones.count(F.boundaryZone[f])) {
                // Dirichlet φ = 0 at wall.
                const double dx = F.centroidX[f] - C.centroidX[o];
                const double dy = F.centroidY[f] - C.centroidY[o];
                const double dz = F.centroidZ[f] - C.centroidZ[o];
                const double dn = std::abs(dx * Ax + dy * Ay + dz * Az);
                const double D = Amag2 / std::max(dn, 1e-30);
                add(o, o, D); // RHS contribution is 0 (φ_b = 0)
            }
            // Other boundaries: ∂φ/∂n = 0 (no contribution).
        }

        util::aligned_vector<double> phi(nC, 0.0);
        sol_.solve(A, b, phi);

        // d = -|∇φ| + √(|∇φ|² + 2φ)
        LeastSquaresGradient G(mesh);
        VectorField gPhi;
        G.evaluate(phi, gPhi);
        out.resize(nC);
        for (std::size_t c = 0; c < nC; ++c) {
            const double gmag =
                std::sqrt(gPhi.x[c] * gPhi.x[c] + gPhi.y[c] * gPhi.y[c] + gPhi.z[c] * gPhi.z[c]);
            const double inside = gmag * gmag + 2.0 * std::max(phi[c], 0.0);
            out[c] = -gmag + std::sqrt(inside);
            if (out[c] < 0)
                out[c] = 0;
        }
        SIMALL_LOG_INFO("Solver", "Wall distance (Poisson-Eikonal): solved for ", nC, " cells");
    }

private:
    ILinearSolver& sol_;
};

} // namespace

std::unique_ptr<IWallDistance> make_wall_distance_exact()
{
    return std::make_unique<ExactNearest>();
}
std::unique_ptr<IWallDistance> make_wall_distance_poisson(ILinearSolver& s)
{
    return std::make_unique<PoissonWallDistance>(s);
}

} // namespace simall::solver
