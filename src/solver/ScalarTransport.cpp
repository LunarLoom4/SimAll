// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/ScalarTransport.cpp
// =============================================================================
#include "solver/ScalarTransport.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::solver
{

namespace
{
inline int find_col(const CSRMatrix& A, int row, int col)
{
    for (int k = A.rowPtr[row]; k < A.rowPtr[row + 1]; ++k)
        if (A.colIdx[k] == col)
            return k;
    return -1;
}
inline void add_at(CSRMatrix& A, int r, int c, double v)
{
    const int k = find_col(A, r, c);
    if (k >= 0)
        A.values[k] += v;
}
const ScalarBC* find_bc(const std::vector<ScalarBC>& v, meshing::ZoneId z)
{
    for (const auto& b : v)
        if (b.zone == z)
            return &b;
    return nullptr;
}

// ---------------------------- TVD limiters --------------------------------
inline double psi_vanLeer(double r)
{
    return (r + std::abs(r)) / (1.0 + std::abs(r) + 1e-30);
}
inline double psi_superBee(double r)
{
    return std::max({0.0, std::min(2.0 * r, 1.0), std::min(r, 2.0)});
}
inline double psi_muscl(double r)
{
    // Symmetric MUSCL (van Leer-MUSCL family) k = 1/3 -> kappa scheme,
    // but a generic minmod-2 form is robust on unstructured meshes:
    return std::max(0.0, std::min({2.0 * r, 0.5 * (1.0 + r), 2.0}));
}
inline double limiter(SpatialScheme s, double r)
{
    switch (s) {
    case SpatialScheme::MUSCL:
        return psi_muscl(r);
    case SpatialScheme::SuperBeeTVD:
        return psi_superBee(r);
    case SpatialScheme::VanLeerTVD:
        return psi_vanLeer(r);
    default:
        return 0.0; // pure upwind
    }
}
} // namespace

ScalarTransport::ScalarTransport(meshing::Mesh& m, FieldRegistry& f, ILinearSolver& l)
    : mesh_(m), F_(f), lin_(l)
{
}

void ScalarTransport::set_zone_bc(ScalarBC bc)
{
    for (auto& b : bcs_)
        if (b.zone == bc.zone) {
            b = bc;
            return;
        }
    bcs_.push_back(bc);
}

void ScalarTransport::build_sparsity()
{
    const auto& C = mesh_.cells();
    const auto& Ff = mesh_.faces();
    const std::size_t nC = C.size();
    A_.rowPtr.assign(nC + 1, 0);
    A_.colIdx.clear();
    for (std::size_t c = 0; c < nC; ++c) {
        std::vector<int> ns{static_cast<int>(c)};
        const int beg = C.faceOffsets[c], end = C.faceOffsets[c + 1];
        for (int k = beg; k < end; ++k) {
            const meshing::FaceId fid = C.faceIndices[k];
            const meshing::CellId oth = (Ff.owner[fid] == c) ? Ff.neighbor[fid] : Ff.owner[fid];
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

double ScalarTransport::solve_iteration()
{
    if (!sparsity_built_)
        build_sparsity();
    const auto& C = mesh_.cells();
    const auto& Ff = mesh_.faces();
    const auto& U = *F_.find_vector("U");
    auto& phi = *F_.find_scalar(name_);
    const std::size_t nC = C.size();

    std::fill(A_.values.begin(), A_.values.end(), 0.0);
    rhs_.assign(nC, 0.0);

    // High-order deferred-correction: pre-compute ∇φ once per iteration
    // when a TVD/MUSCL scheme is requested. Matrix stays 1st-order upwind
    // (→ diagonally dominant); high-order minus upwind face value goes to
    // rhs so the converged solution carries 2nd-order spatial accuracy.
    const bool useTVD = (scheme_ == SpatialScheme::MUSCL || scheme_ == SpatialScheme::SuperBeeTVD
                         || scheme_ == SpatialScheme::VanLeerTVD);
    VectorField gPhi;
    if (useTVD) {
        LeastSquaresGradient G(mesh_);
        G.evaluate(phi, gPhi);
    }

    for (std::size_t f = 0; f < Ff.size(); ++f) {
        const meshing::CellId o = Ff.owner[f];
        const meshing::CellId n = Ff.neighbor[f];
        const double Ax = Ff.areaX[f], Ay = Ff.areaY[f], Az = Ff.areaZ[f];
        const double Amag2 = Ax * Ax + Ay * Ay + Az * Az;
        double Ufx, Ufy, Ufz;
        if (n != meshing::kBoundaryCell) {
            Ufx = 0.5 * (U.x[o] + U.x[n]);
            Ufy = 0.5 * (U.y[o] + U.y[n]);
            Ufz = 0.5 * (U.z[o] + U.z[n]);
        } else {
            Ufx = U.x[o];
            Ufy = U.y[o];
            Ufz = U.z[o];
        }
        const double Fm = rho_ * (Ufx * Ax + Ufy * Ay + Ufz * Az);
        const double ConvO = std::max(Fm, 0.0);
        const double ConvN = std::max(-Fm, 0.0);

        if (n != meshing::kBoundaryCell) {
            const double dx = C.centroidX[n] - C.centroidX[o];
            const double dy = C.centroidY[n] - C.centroidY[o];
            const double dz = C.centroidZ[n] - C.centroidZ[o];
            const double dn = std::abs(dx * Ax + dy * Ay + dz * Az);
            const double D = gamma_ * Amag2 / std::max(dn, 1e-30);
            add_at(A_, o, o, ConvO + D);
            add_at(A_, o, n, -(ConvN + D));
            add_at(A_, n, n, ConvN + D);
            add_at(A_, n, o, -(ConvO + D));

            // -------- TVD / MUSCL deferred correction (internal face) --------
            if (useTVD) {
                // Upwind / downwind cells determined by mass flux direction.
                const meshing::CellId U = (Fm >= 0) ? o : n;
                const meshing::CellId D2 = (Fm >= 0) ? n : o;
                const double phiU = phi[U], phiD = phi[D2];
                const double dPhi = phiD - phiU;
                // r = (∇φ_U · d) / (φ_D - φ_U) where d = x_D - x_U
                const double rdx = C.centroidX[D2] - C.centroidX[U];
                const double rdy = C.centroidY[D2] - C.centroidY[U];
                const double rdz = C.centroidZ[D2] - C.centroidZ[U];
                const double gdotd = gPhi.x[U] * rdx + gPhi.y[U] * rdy + gPhi.z[U] * rdz;
                // 2 ∇φ_U·d / Δφ - 1   (Sweby r definition for unstructured cells)
                const double r = (std::abs(dPhi) > 1e-30) ? (2.0 * gdotd / dPhi - 1.0) : 0.0;
                const double psi = limiter(scheme_, r);
                // Face value (cell-to-face vector ≈ ½ d for centred placement):
                const double phiHO = phiU + 0.5 * psi * dPhi;
                const double dPhi_face = phiHO - phiU; // high-order - upwind
                // Move flux correction Fm * dPhi_face from matrix to RHS:
                // for owner: +Fm * φ_face contributes; we already used +Fm * φ_U.
                rhs_[o] -= Fm * dPhi_face;
                rhs_[n] += Fm * dPhi_face;
            }
        } else {
            const auto* bc = find_bc(bcs_, Ff.boundaryZone[f]);
            const double dx = Ff.centroidX[f] - C.centroidX[o];
            const double dy = Ff.centroidY[f] - C.centroidY[o];
            const double dz = Ff.centroidZ[f] - C.centroidZ[o];
            const double dn = std::abs(dx * Ax + dy * Ay + dz * Az);
            const double D = gamma_ * Amag2 / std::max(dn, 1e-30);
            if (!bc || bc->kind == ScalarBC::Kind::Neumann) {
                rhs_[o] += (bc ? bc->value : 0.0) * std::sqrt(Amag2);
                add_at(A_, o, o, ConvO);
            } else if (bc->kind == ScalarBC::Kind::Dirichlet) {
                add_at(A_, o, o, ConvO + D);
                rhs_[o] += (ConvN + D) * bc->value;
            } else { // Robin: h(φ - φ_∞)
                const double h = bc->value;
                const double Afmag = std::sqrt(Amag2);
                add_at(A_, o, o, ConvO + h * Afmag);
                rhs_[o] += h * Afmag * bc->valueB;
            }
        }
    }

    // Optional explicit source term (e.g., viscous heating).
    if (!sourceField_.empty()) {
        if (auto* S = F_.find_scalar(sourceField_)) {
            for (std::size_t c = 0; c < nC; ++c)
                rhs_[c] += (*S)[c] * C.volume[c];
        }
    }

    // Patankar under-relaxation.
    for (std::size_t c = 0; c < nC; ++c) {
        const int kd = find_col(A_, c, c);
        if (kd < 0)
            continue;
        const double a = A_.values[kd];
        const double a2 = a / urf_;
        A_.values[kd] = a2;
        rhs_[c] += (a2 - a) * phi[c];
    }

    util::aligned_vector<double> sol = phi;
    lin_.solve(A_, rhs_, sol);
    phi = std::move(sol);

    util::aligned_vector<double> Ax(nC, 0);
    A_.spmv(phi, Ax);
    double r = 0;
    for (std::size_t i = 0; i < nC; ++i) {
        const double d = Ax[i] - rhs_[i];
        r += d * d;
    }
    return std::sqrt(r);
}

} // namespace simall::solver
