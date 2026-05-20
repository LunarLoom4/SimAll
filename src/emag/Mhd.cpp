// =============================================================================
// SimAll Beta - Electromagnetics Subsystem
// File   : src/emag/Mhd.cpp
//
// Implementation strategy
//   φ Poisson:   ∇·(σ ∇φ) = ∇·(σ u × B)
//   discretised face-by-face on the polyhedral mesh with central diffusion
//   plus the (u × B) flux divergence as a right-hand side.
//
// After the φ solve the cell-centred current density follows from
//   J = σ (-∇φ + u × B)
// using a least-squares gradient of φ. Finally the Lorentz body force
//   f = J × B
// is written to the FieldRegistry vector channel "S_Lorentz" for the
// momentum solver to consume on its next outer iteration.
// =============================================================================
#include "emag/Mhd.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::emag
{

namespace
{
inline int find_col(const solver::CSRMatrix& A, int r, int c)
{
    for (int k = A.rowPtr[r]; k < A.rowPtr[r + 1]; ++k)
        if (A.colIdx[k] == c)
            return k;
    return -1;
}
inline void add_at(solver::CSRMatrix& A, int r, int c, double v)
{
    const int k = find_col(A, r, c);
    if (k >= 0)
        A.values[k] += v;
}
const MhdProps::PhiBC* find_bc(const std::vector<MhdProps::PhiBC>& v, meshing::ZoneId z)
{
    for (const auto& b : v)
        if (b.zone == z)
            return &b;
    return nullptr;
}
} // namespace

Mhd::Mhd(meshing::Mesh& m, solver::FieldRegistry& f, solver::ILinearSolver& l)
    : mesh_(m), F_(f), lin_(l)
{
}

void Mhd::initialize(MhdProps p)
{
    props_ = std::move(p);
    const std::size_t nC = mesh_.cells().size();
    F_.scalar("phi_e", nC);
    F_.vector("J_e", nC);
    F_.vector("S_Lorentz", nC);
}

void Mhd::build_sparsity()
{
    const auto& C = mesh_.cells();
    const auto& F = mesh_.faces();
    const std::size_t nC = C.size();
    A_.rowPtr.assign(nC + 1, 0);
    A_.colIdx.clear();
    for (std::size_t c = 0; c < nC; ++c) {
        std::vector<int> ns{static_cast<int>(c)};
        const int s = C.faceOffsets[c], e = C.faceOffsets[c + 1];
        for (int k = s; k < e; ++k) {
            const meshing::FaceId fid = C.faceIndices[k];
            const meshing::CellId oth = (F.owner[fid] == c) ? F.neighbor[fid] : F.owner[fid];
            if (oth != meshing::kBoundaryCell)
                ns.push_back(static_cast<int>(oth));
        }
        std::sort(ns.begin(), ns.end());
        ns.erase(std::unique(ns.begin(), ns.end()), ns.end());
        for (int n : ns)
            A_.colIdx.push_back(n);
        A_.rowPtr[c + 1] = static_cast<int>(A_.colIdx.size());
    }
    A_.values.assign(A_.colIdx.size(), 0.0);
    rhs_.assign(nC, 0.0);
    sparsity_built_ = true;
}

double Mhd::solve_iteration()
{
    if (!sparsity_built_)
        build_sparsity();
    const auto& C = mesh_.cells();
    const auto& Ff = mesh_.faces();
    const std::size_t nC = C.size();

    const auto* U = F_.find_vector("U");
    auto& phi = *F_.find_scalar("phi_e");
    auto& J = *F_.find_vector("J_e");
    auto& Sl = *F_.find_vector("S_Lorentz");
    const double sigma = props_.sigma;

    // Precompute B per cell.
    util::aligned_vector<double> Bx(nC), By(nC), Bz(nC);
    for (std::size_t c = 0; c < nC; ++c) {
        util::Vec3d B = props_.Bfield
                            ? props_.Bfield({C.centroidX[c], C.centroidY[c], C.centroidZ[c]})
                            : util::Vec3d{0, 0, 0};
        Bx[c] = B.x;
        By[c] = B.y;
        Bz[c] = B.z;
    }

    std::fill(A_.values.begin(), A_.values.end(), 0.0);
    std::fill(rhs_.begin(), rhs_.end(), 0.0);

    // Assemble Poisson with RHS from ∇·(σ u × B).
    for (std::size_t f = 0; f < Ff.size(); ++f) {
        const meshing::CellId o = Ff.owner[f];
        const meshing::CellId n = Ff.neighbor[f];
        const double Ax = Ff.areaX[f], Ay = Ff.areaY[f], Az = Ff.areaZ[f];
        const double Amag2 = Ax * Ax + Ay * Ay + Az * Az;
        // Face-averaged u × B (face centre via owner / neighbour blend).
        double uXBx, uXBy, uXBz;
        if (U) {
            double ux, uy, uz, bx, by, bz;
            if (n != meshing::kBoundaryCell) {
                ux = 0.5 * (U->x[o] + U->x[n]);
                uy = 0.5 * (U->y[o] + U->y[n]);
                uz = 0.5 * (U->z[o] + U->z[n]);
                bx = 0.5 * (Bx[o] + Bx[n]);
                by = 0.5 * (By[o] + By[n]);
                bz = 0.5 * (Bz[o] + Bz[n]);
            } else {
                ux = U->x[o];
                uy = U->y[o];
                uz = U->z[o];
                bx = Bx[o];
                by = By[o];
                bz = Bz[o];
            }
            uXBx = uy * bz - uz * by;
            uXBy = uz * bx - ux * bz;
            uXBz = ux * by - uy * bx;
        } else {
            uXBx = uXBy = uXBz = 0.0;
        }
        const double srcFlux = sigma * (uXBx * Ax + uXBy * Ay + uXBz * Az);

        if (n != meshing::kBoundaryCell) {
            const double dx = C.centroidX[n] - C.centroidX[o];
            const double dy = C.centroidY[n] - C.centroidY[o];
            const double dz = C.centroidZ[n] - C.centroidZ[o];
            const double dn = std::abs(dx * Ax + dy * Ay + dz * Az);
            const double D = sigma * Amag2 / std::max(dn, 1e-30);
            add_at(A_, o, o, +D);
            add_at(A_, o, n, -D);
            add_at(A_, n, n, +D);
            add_at(A_, n, o, -D);
            rhs_[o] += srcFlux;
            rhs_[n] -= srcFlux;
        } else {
            const auto* bc = find_bc(props_.bcs, Ff.boundaryZone[f]);
            const double dx = Ff.centroidX[f] - C.centroidX[o];
            const double dy = Ff.centroidY[f] - C.centroidY[o];
            const double dz = Ff.centroidZ[f] - C.centroidZ[o];
            const double dn = std::abs(dx * Ax + dy * Ay + dz * Az);
            const double D = sigma * Amag2 / std::max(dn, 1e-30);
            if (bc && bc->isDirichlet) {
                add_at(A_, o, o, +D);
                rhs_[o] += D * bc->value + srcFlux;
            } else {
                // Insulator (Neumann): J·n = 0 ⇒ -σ ∂φ/∂n + σ (u×B)·n = 0
                // Move RHS contribution only.
                rhs_[o] += srcFlux;
            }
        }
    }

    // Pin solution if no Dirichlet anywhere (singular Poisson).
    bool anyDirichlet = false;
    for (const auto& b : props_.bcs)
        if (b.isDirichlet) {
            anyDirichlet = true;
            break;
        }
    if (!anyDirichlet) {
        const int k = find_col(A_, 0, 0);
        if (k >= 0)
            A_.values[k] += 1e30;
    }

    util::aligned_vector<double> sol = phi;
    lin_.solve(A_, rhs_, sol);
    phi = std::move(sol);

    // J = σ (-∇φ + u × B)
    solver::LeastSquaresGradient G(mesh_);
    solver::VectorField gphi;
    G.evaluate(phi, gphi);
    for (std::size_t c = 0; c < nC; ++c) {
        const double ux = U ? U->x[c] : 0;
        const double uy = U ? U->y[c] : 0;
        const double uz = U ? U->z[c] : 0;
        const double uxBx = uy * Bz[c] - uz * By[c];
        const double uxBy = uz * Bx[c] - ux * Bz[c];
        const double uxBz = ux * By[c] - uy * Bx[c];
        J.x[c] = sigma * (-gphi.x[c] + uxBx);
        J.y[c] = sigma * (-gphi.y[c] + uxBy);
        J.z[c] = sigma * (-gphi.z[c] + uxBz);
        // Lorentz body force f = J × B
        Sl.x[c] = J.y[c] * Bz[c] - J.z[c] * By[c];
        Sl.y[c] = J.z[c] * Bx[c] - J.x[c] * Bz[c];
        Sl.z[c] = J.x[c] * By[c] - J.y[c] * Bx[c];
    }

    // Diagnostic residual.
    util::aligned_vector<double> Ax(nC, 0);
    A_.spmv(phi, Ax);
    double r = 0;
    for (std::size_t i = 0; i < nC; ++i) {
        const double d = Ax[i] - rhs_[i];
        r += d * d;
    }
    return std::sqrt(r);
}

} // namespace simall::emag
