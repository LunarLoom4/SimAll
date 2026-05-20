// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/FluxCorrectedTransport.cpp
// =============================================================================
#include "solver/FluxCorrectedTransport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace simall::solver
{

namespace
{
const BoundarySpec* find_bc(const std::vector<BoundarySpec>& bcs, meshing::ZoneId z)
{
    for (const auto& b : bcs)
        if (b.zone == z)
            return &b;
    return nullptr;
}

// Boundary face value of phi: Dirichlet from BC, else owner extrapolation.
inline double boundary_phi(const BoundarySpec* bc, double phiOwner)
{
    if (!bc)
        return phiOwner;
    if (bc->type == BCType::VelocityInlet || bc->type == BCType::PressureInlet
        || bc->type == BCType::MassFlowInlet) {
        return bc->scalarValue; // BC-supplied scalar value at inlet
    }
    return phiOwner; // outflow / wall / symmetry → zero-grad
}
} // namespace

void FluxCorrectedTransport::compute_limiter(const util::aligned_vector<double>& phi,
                                             const util::aligned_vector<double>& massFlux,
                                             double dt,
                                             util::aligned_vector<double>& alpha)
{
    const auto& F = mesh_.faces();
    const auto& C = mesh_.cells();
    const std::size_t nC = C.size();
    const std::size_t nF = F.size();
    alpha.assign(nF, 1.0);

    // ---- 1) Low-order upwind flux per face: F_L = max(F,0) φ_O + min(F,0) φ_N
    util::aligned_vector<double> FL(nF, 0.0), FH(nF, 0.0), A(nF, 0.0);
    for (std::size_t f = 0; f < nF; ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        const double Ff = massFlux[f];
        double phiN;
        if (n != meshing::kBoundaryCell)
            phiN = phi[n];
        else
            phiN = boundary_phi(find_bc(bcs_, F.boundaryZone[f]), phi[o]);
        // Upwind
        FL[f] = (Ff >= 0.0) ? Ff * phi[o] : Ff * phiN;
        // Central
        FH[f] = Ff * 0.5 * (phi[o] + phiN);
        A[f] = FH[f] - FL[f];
    }

    // ---- 2) Low-order advanced field φ^low ----
    util::aligned_vector<double> phiLow = phi;
    for (std::size_t f = 0; f < nF; ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        const double dV_o = dt / std::max(C.volume[o], 1.0e-30);
        phiLow[o] -= dV_o * FL[f];
        if (n != meshing::kBoundaryCell) {
            const double dV_n = dt / std::max(C.volume[n], 1.0e-30);
            phiLow[n] += dV_n * FL[f];
        }
    }

    // ---- 3) Per-cell antidiffusive accumulators P+, P-  +  bounds φ_max, φ_min ----
    util::aligned_vector<double> Pplus(nC, 0.0), Pminus(nC, 0.0);
    util::aligned_vector<double> phiMax(nC), phiMin(nC);
    for (std::size_t c = 0; c < nC; ++c) {
        phiMax[c] = phiMin[c] = phi[c];
    }
    // Initial-bound widening by neighbour values
    for (std::size_t f = 0; f < nF; ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        double phiN;
        if (n != meshing::kBoundaryCell)
            phiN = phi[n];
        else
            phiN = boundary_phi(find_bc(bcs_, F.boundaryZone[f]), phi[o]);
        if (phiN > phiMax[o])
            phiMax[o] = phiN;
        if (phiN < phiMin[o])
            phiMin[o] = phiN;
        if (n != meshing::kBoundaryCell) {
            if (phi[o] > phiMax[n])
                phiMax[n] = phi[o];
            if (phi[o] < phiMin[n])
                phiMin[n] = phi[o];
        }
    }
    // Accumulate P+/-: positive A_f *received* by cell increases P+; outgoing increases P-.
    // Convention: A_f acts owner→neighbour. So owner LOSES A_f, neighbour GAINS A_f.
    for (std::size_t f = 0; f < nF; ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        const double Af = A[f];
        if (Af >= 0.0) {
            // Owner sees outgoing antidiffusion of magnitude |Af|
            Pminus[o] += Af;
            if (n != meshing::kBoundaryCell)
                Pplus[n] += Af;
        } else {
            Pplus[o] += -Af;
            if (n != meshing::kBoundaryCell)
                Pminus[n] += -Af;
        }
    }

    // ---- 4) Per-cell Q+/-, R+/- ----
    util::aligned_vector<double> Rplus(nC, 1.0), Rminus(nC, 1.0);
    for (std::size_t c = 0; c < nC; ++c) {
        const double Vdt = C.volume[c] / std::max(dt, 1.0e-30);
        const double Qp = (phiMax[c] - phiLow[c]) * Vdt;
        const double Qm = (phiLow[c] - phiMin[c]) * Vdt;
        Rplus[c] = (Pplus[c] > 1.0e-30) ? std::min(1.0, Qp / Pplus[c]) : 1.0;
        Rminus[c] = (Pminus[c] > 1.0e-30) ? std::min(1.0, Qm / Pminus[c]) : 1.0;
        Rplus[c] = std::max(0.0, Rplus[c]);
        Rminus[c] = std::max(0.0, Rminus[c]);
    }

    // ---- 5) Face limiter α_f ----
    for (std::size_t f = 0; f < nF; ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        const double Af = A[f];
        double a;
        if (Af >= 0.0) {
            // Antidiffusion sends from o to n.
            const double Rn = (n != meshing::kBoundaryCell) ? Rplus[n] : 1.0;
            a = std::min(Rn, Rminus[o]);
        } else {
            const double Rn = (n != meshing::kBoundaryCell) ? Rminus[n] : 1.0;
            a = std::min(Rplus[o], Rn);
        }
        alpha[f] = std::clamp(a, 0.0, 1.0);
    }
}

void FluxCorrectedTransport::advance(util::aligned_vector<double>& phi,
                                     const util::aligned_vector<double>& massFlux,
                                     double dt)
{
    const auto& F = mesh_.faces();
    const auto& C = mesh_.cells();
    const std::size_t nF = F.size();

    util::aligned_vector<double> alpha;
    compute_limiter(phi, massFlux, dt, alpha);

    // Recompute fluxes and apply F = F_L + α (F_H - F_L)
    for (std::size_t f = 0; f < nF; ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        const double Ff = massFlux[f];
        double phiN;
        if (n != meshing::kBoundaryCell)
            phiN = phi[n];
        else
            phiN = boundary_phi(find_bc(bcs_, F.boundaryZone[f]), phi[o]);
        const double FL = (Ff >= 0.0) ? Ff * phi[o] : Ff * phiN;
        const double FH = Ff * 0.5 * (phi[o] + phiN);
        const double Fcorr = FL + alpha[f] * (FH - FL);
        const double dV_o = dt / std::max(C.volume[o], 1.0e-30);
        phi[o] -= dV_o * Fcorr;
        if (n != meshing::kBoundaryCell) {
            const double dV_n = dt / std::max(C.volume[n], 1.0e-30);
            phi[n] += dV_n * Fcorr;
        }
    }
}

} // namespace simall::solver
