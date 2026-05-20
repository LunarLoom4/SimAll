// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/SimpleAlgorithm.cpp
//
// Reference: Patankar (1980); Ferziger & Perić (2002) ch 7-8; Jasak PhD (1996).
// =============================================================================
#include "solver/SimpleAlgorithm.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

namespace simall::solver
{

namespace
{
// Linear-search lookup of column index inside a CSR row (rows are short).
inline int find_col(const CSRMatrix& A, int row, int col)
{
    for (int k = A.rowPtr[row]; k < A.rowPtr[row + 1]; ++k)
        if (A.colIdx[k] == col)
            return k;
    return -1;
}
inline void add_at(CSRMatrix& A, int row, int col, double v)
{
    const int k = find_col(A, row, col);
    if (k >= 0)
        A.values[k] += v;
}
} // namespace

SimpleAlgorithm::SimpleAlgorithm(meshing::Mesh& m,
                                 FieldRegistry& f,
                                 const std::vector<BoundarySpec>& b,
                                 ILinearSolver& sm,
                                 ILinearSolver& sp,
                                 SimpleOptions o)
    : mesh_(m), F_(f), bcs_(b), linMom_(sm), linP_(sp), opt_(o)
{
}

// ============================================================ sparsity
void SimpleAlgorithm::build_sparsity()
{
    const auto& C = mesh_.cells();
    const auto& Ff = mesh_.faces();
    const std::size_t nC = C.size();
    // Periodic pairs must be known before sparsity so twin cells are
    // included in each owner's column list.
    build_periodic_pairs();

    A_mom_.rowPtr.assign(nC + 1, 0);
    A_p_.rowPtr.assign(nC + 1, 0);
    A_mom_.colIdx.clear();
    A_p_.colIdx.clear();

    for (std::size_t c = 0; c < nC; ++c) {
        std::vector<int> nbrs{static_cast<int>(c)};
        const int beg = C.faceOffsets[c], end = C.faceOffsets[c + 1];
        for (int k = beg; k < end; ++k) {
            const meshing::FaceId fid = C.faceIndices[k];
            const meshing::CellId oth = (Ff.owner[fid] == c) ? Ff.neighbor[fid] : Ff.owner[fid];
            if (oth != meshing::kBoundaryCell) {
                nbrs.push_back(static_cast<int>(oth));
            } else if (!periodicTwin_.empty() && periodicTwin_[fid] >= 0) {
                nbrs.push_back(static_cast<int>(Ff.owner[periodicTwin_[fid]]));
            }
        }
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
        for (int n : nbrs) {
            A_mom_.colIdx.push_back(n);
            A_p_.colIdx.push_back(n);
        }
        A_mom_.rowPtr[c + 1] = static_cast<int>(A_mom_.colIdx.size());
        A_p_.rowPtr[c + 1] = A_mom_.rowPtr[c + 1];
    }
    A_mom_.values.assign(A_mom_.colIdx.size(), 0.0);
    A_p_.values.assign(A_p_.colIdx.size(), 0.0);
    aP_.assign(nC, 0.0);
    rhs_.assign(nC, 0.0);
    sol_.assign(nC, 0.0);
    pPrime_.assign(nC, 0.0);
    // BDF2 history (zero == cold start, falls back to implicit Euler on step 1)
    Ux_n_.assign(nC, 0.0);
    Uy_n_.assign(nC, 0.0);
    Uz_n_.assign(nC, 0.0);
    Ux_nm1_.assign(nC, 0.0);
    Uy_nm1_.assign(nC, 0.0);
    Uz_nm1_.assign(nC, 0.0);
    sparsity_built_ = true;
}

// ============================================================ periodic pairs
void SimpleAlgorithm::build_periodic_pairs()
{
    const auto& F = mesh_.faces();
    if (externalTwinSet_ && externalTwin_.size() == F.size()) {
        periodicTwin_ = externalTwin_;
        SIMALL_LOG_INFO("Solver",
                        "Periodic pairing: using externally supplied twin table (",
                        F.size(),
                        " faces)");
        return;
    }
    periodicTwin_.assign(F.size(), -1);

    // Collect periodic faces by zone.
    std::vector<int> periodicZones;
    for (const auto& b : bcs_)
        if (b.type == BCType::Periodic)
            periodicZones.push_back(static_cast<int>(b.zone));
    if (periodicZones.size() < 2)
        return;

    // Group all periodic boundary faces by zone.
    std::vector<std::vector<int>> facesByZone(periodicZones.size());
    auto zoneIdx = [&](meshing::ZoneId z) -> int {
        for (std::size_t i = 0; i < periodicZones.size(); ++i)
            if (periodicZones[i] == static_cast<int>(z))
                return static_cast<int>(i);
        return -1;
    };
    for (std::size_t f = 0; f < F.size(); ++f) {
        if (F.neighbor[f] != meshing::kBoundaryCell)
            continue;
        const int zi = zoneIdx(F.boundaryZone[f]);
        if (zi >= 0)
            facesByZone[zi].push_back(static_cast<int>(f));
    }

    // Pair zones two-by-two: (0,1), (2,3), (4,5). Match faces by centroid
    // after projecting out the periodic translation vector. Heuristic: pair
    // the closest tangential-coordinate face.
    for (std::size_t z = 0; z + 1 < periodicZones.size(); z += 2) {
        const auto& A = facesByZone[z];
        const auto& B = facesByZone[z + 1];
        if (A.empty() || B.empty())
            continue;
        // Compute translation as mean(B) - mean(A).
        double Ax = 0, Ay = 0, Az = 0, Bx = 0, By = 0, Bz = 0;
        for (int f : A) {
            Ax += F.centroidX[f];
            Ay += F.centroidY[f];
            Az += F.centroidZ[f];
        }
        for (int f : B) {
            Bx += F.centroidX[f];
            By += F.centroidY[f];
            Bz += F.centroidZ[f];
        }
        const double invA = 1.0 / A.size();
        const double invB = 1.0 / B.size();
        const double tx = Bx * invB - Ax * invA, ty = By * invB - Ay * invA,
                     tz = Bz * invB - Az * invA;
        // For each face in A, find closest in B after subtracting t.
        for (int fa : A) {
            int best = -1;
            double bestd2 = std::numeric_limits<double>::max();
            const double ax = F.centroidX[fa] + tx;
            const double ay = F.centroidY[fa] + ty;
            const double az = F.centroidZ[fa] + tz;
            for (int fb : B) {
                const double d2 = (F.centroidX[fb] - ax) * (F.centroidX[fb] - ax)
                                  + (F.centroidY[fb] - ay) * (F.centroidY[fb] - ay)
                                  + (F.centroidZ[fb] - az) * (F.centroidZ[fb] - az);
                if (d2 < bestd2) {
                    bestd2 = d2;
                    best = fb;
                }
            }
            if (best >= 0) {
                periodicTwin_[fa] = best;
                periodicTwin_[best] = fa;
            }
        }
    }
    SIMALL_LOG_INFO(
        "Solver", "Periodic face pairs established for ", periodicZones.size(), " zones");
}

// ============================================================ boundary helpers
namespace
{
const BoundarySpec* find_bc(const std::vector<BoundarySpec>& bcs, meshing::ZoneId z)
{
    for (const auto& b : bcs)
        if (b.zone == z)
            return &b;
    return nullptr;
}
} // namespace

// ============================================================ mass fluxes
void SimpleAlgorithm::compute_face_fluxes()
{
    const auto& F = mesh_.faces();
    const auto& U = *F_.find_vector("U");
    phiFlux_.assign(F.size(), 0.0);
    for (std::size_t f = 0; f < F.size(); ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        double Ux, Uy, Uz;
        if (n != meshing::kBoundaryCell) {
            Ux = 0.5 * (U.x[o] + U.x[n]);
            Uy = 0.5 * (U.y[o] + U.y[n]);
            Uz = 0.5 * (U.z[o] + U.z[n]);
        } else if (!periodicTwin_.empty() && periodicTwin_[f] >= 0) {
            const meshing::CellId twin = F.owner[periodicTwin_[f]];
            Ux = 0.5 * (U.x[o] + U.x[twin]);
            Uy = 0.5 * (U.y[o] + U.y[twin]);
            Uz = 0.5 * (U.z[o] + U.z[twin]);
        } else {
            const auto* bc = find_bc(bcs_, F.boundaryZone[f]);
            if (bc && (bc->type == BCType::VelocityInlet || bc->type == BCType::MovingWall)) {
                Ux = bc->vectorValue[0];
                Uy = bc->vectorValue[1];
                Uz = bc->vectorValue[2];
            } else if (bc
                       && (bc->type == BCType::NoSlipWall || bc->type == BCType::Wall
                           || bc->type == BCType::Symmetry)) {
                Ux = Uy = Uz = 0.0;
            } else {
                Ux = U.x[o];
                Uy = U.y[o];
                Uz = U.z[o]; // outflow extrapolation
            }
        }
        phiFlux_[f] = opt_.rho * (Ux * F.areaX[f] + Uy * F.areaY[f] + Uz * F.areaZ[f]);
    }
}

// ============================================================ momentum
void SimpleAlgorithm::assemble_momentum(int comp)
{
    const auto& C = mesh_.cells();
    const auto& F = mesh_.faces();
    const auto& U = *F_.find_vector("U");
    const auto& p = *F_.find_scalar("p");
    const std::size_t nC = C.size();

    std::fill(A_mom_.values.begin(), A_mom_.values.end(), 0.0);
    rhs_.assign(nC, 0.0);

    const auto& Ucomp = (comp == 0) ? U.x : (comp == 1) ? U.y : U.z;

    for (std::size_t f = 0; f < F.size(); ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
        const double Amag2 = Ax * Ax + Ay * Ay + Az * Az;
        const double Ff = phiFlux_[f];
        const double ConvO = std::max(Ff, 0.0);
        const double ConvN = std::max(-Ff, 0.0);

        if (n != meshing::kBoundaryCell) {
            const double dx = C.centroidX[n] - C.centroidX[o];
            const double dy = C.centroidY[n] - C.centroidY[o];
            const double dz = C.centroidZ[n] - C.centroidZ[o];
            const double ddotA = std::abs(dx * Ax + dy * Ay + dz * Az);
            const double D = opt_.mu * Amag2 / std::max(ddotA, 1e-30);
            add_at(A_mom_, o, o, ConvO + D);
            add_at(A_mom_, o, n, -(ConvN + D));
            add_at(A_mom_, n, n, ConvN + D);
            add_at(A_mom_, n, o, -(ConvO + D));
        } else {
            const auto* bc = find_bc(bcs_, F.boundaryZone[f]);
            const double dx = F.centroidX[f] - C.centroidX[o];
            const double dy = F.centroidY[f] - C.centroidY[o];
            const double dz = F.centroidZ[f] - C.centroidZ[o];
            const double ddotA = std::abs(dx * Ax + dy * Ay + dz * Az);
            const double D = opt_.mu * Amag2 / std::max(ddotA, 1e-30);
            if (!periodicTwin_.empty() && periodicTwin_[f] >= 0) {
                // Periodic: paired with cell on opposite side. Treat as
                // half-of-internal-face coupling (the twin face supplies
                // the other half on its own row).
                const meshing::CellId tw = F.owner[periodicTwin_[f]];
                add_at(A_mom_, static_cast<int>(o), static_cast<int>(o), ConvO + D);
                add_at(A_mom_, static_cast<int>(o), static_cast<int>(tw), -(ConvN + D));
            } else if (bc
                       && (bc->type == BCType::NoSlipWall || bc->type == BCType::Wall
                           || bc->type == BCType::VelocityInlet
                           || bc->type == BCType::MovingWall)) {
                const double Ub = (bc->type == BCType::NoSlipWall || bc->type == BCType::Wall)
                                      ? 0.0
                                      : bc->vectorValue[comp];
                add_at(A_mom_, o, o, ConvO + D);
                rhs_[o] += (ConvN + D) * Ub;
            } else if (bc && bc->type == BCType::Symmetry) {
                // Zero normal gradient → no contribution to diagonal or rhs
                // for the parallel components; normal component is zeroed in
                // post-correction.
            } else {
                // Outflow / pressure-outlet: zero-gradient → no extra coeff.
                add_at(A_mom_, o, o, ConvO); // upwind only
            }
        }
    }

    // ---- pressure gradient source: -∂p/∂x_k * V_c -------------------------
    {
        LeastSquaresGradient G(mesh_);
        VectorField gP;
        G.evaluate(p, gP);
        for (std::size_t c = 0; c < nC; ++c) {
            const double gk = (comp == 0) ? gP.x[c] : (comp == 1) ? gP.y[c] : gP.z[c];
            rhs_[c] -= gk * C.volume[c];
        }
    }

    // ---- transient term (BDF2 or implicit Euler) --------------------------
    // BDF2:  (3/2 U^{n+1} - 2 U^n + 1/2 U^{n-1}) / Δt
    // BE  :  (U^{n+1} - U^n) / Δt
    if (opt_.dt > 0.0) {
        const auto& Un = (comp == 0) ? Ux_n_ : (comp == 1) ? Uy_n_ : Uz_n_;
        const auto& Unm1 = (comp == 0) ? Ux_nm1_ : (comp == 1) ? Uy_nm1_ : Uz_nm1_;
        const bool useBDF2 = (opt_.timeScheme == TemporalScheme::BDF2) && (timeStep_ >= 1);
        for (std::size_t c = 0; c < nC; ++c) {
            const double Vrho_dt = opt_.rho * C.volume[c] / opt_.dt;
            const int kd = find_col(A_mom_, static_cast<int>(c), static_cast<int>(c));
            if (kd < 0)
                continue;
            if (useBDF2) {
                A_mom_.values[kd] += 1.5 * Vrho_dt;
                rhs_[c] += (2.0 * Un[c] - 0.5 * Unm1[c]) * Vrho_dt;
            } else {
                A_mom_.values[kd] += Vrho_dt;
                rhs_[c] += Un[c] * Vrho_dt;
            }
        }
    }

    // ---- under-relaxation (Patankar) --------------------------------------
    for (std::size_t c = 0; c < nC; ++c) {
        const int kd = find_col(A_mom_, c, c);
        if (kd < 0)
            continue;
        const double aDiag = A_mom_.values[kd];
        const double scaled = aDiag / opt_.urfU;
        A_mom_.values[kd] = scaled;
        rhs_[c] += (scaled - aDiag) * Ucomp[c];
    }

    // capture diagonal for Rhie-Chow / pressure correction
    if (comp == 0) {
        for (std::size_t c = 0; c < nC; ++c)
            aP_[c] = A_mom_.values[find_col(A_mom_, c, c)];
    }
}

// ============================================================ pressure
void SimpleAlgorithm::assemble_pressure_correction()
{
    const auto& C = mesh_.cells();
    const auto& F = mesh_.faces();
    const std::size_t nC = C.size();
    std::fill(A_p_.values.begin(), A_p_.values.end(), 0.0);
    rhs_.assign(nC, 0.0);
    bool anyDirichlet = false;

    for (std::size_t f = 0; f < F.size(); ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
        const double Amag2 = Ax * Ax + Ay * Ay + Az * Az;

        if (n != meshing::kBoundaryCell) {
            const double dx = C.centroidX[n] - C.centroidX[o];
            const double dy = C.centroidY[n] - C.centroidY[o];
            const double dz = C.centroidZ[n] - C.centroidZ[o];
            const double ddotA = std::abs(dx * Ax + dy * Ay + dz * Az);
            // SIMPLE: 1/aP. SIMPLEC: 1/(aP - Σ a_NB). aP_off-diag-sum is
            // recomputed from the momentum CSR row (sum |off-diagonals|).
            auto denom_for = [&](meshing::CellId c) -> double {
                if (opt_.algorithm == PvCouplingVariant::SIMPLEC) {
                    double sumOff = 0.0;
                    for (int k = A_mom_.rowPtr[c]; k < A_mom_.rowPtr[c + 1]; ++k) {
                        if (A_mom_.colIdx[k] != static_cast<int>(c))
                            sumOff += std::abs(A_mom_.values[k]);
                    }
                    return std::max(aP_[c] - sumOff, 1e-30);
                }
                return std::max(aP_[c], 1e-30);
            };
            const double dOwn = denom_for(o), dNbr = denom_for(n);
            const double aPf = 0.5 * (dOwn + dNbr);
            const double coef = opt_.rho * Amag2 / std::max(aPf * ddotA, 1e-30);
            add_at(A_p_, o, o, coef);
            add_at(A_p_, o, n, -coef);
            add_at(A_p_, n, n, coef);
            add_at(A_p_, n, o, -coef);
        } else {
            const auto* bc = find_bc(bcs_, F.boundaryZone[f]);
            // Same SIMPLEC-aware denominator for boundary contributions.
            auto denom_for_b = [&](meshing::CellId c) -> double {
                if (opt_.algorithm == PvCouplingVariant::SIMPLEC) {
                    double sumOff = 0.0;
                    for (int k = A_mom_.rowPtr[c]; k < A_mom_.rowPtr[c + 1]; ++k)
                        if (A_mom_.colIdx[k] != static_cast<int>(c))
                            sumOff += std::abs(A_mom_.values[k]);
                    return std::max(aP_[c] - sumOff, 1e-30);
                }
                return std::max(aP_[c], 1e-30);
            };
            if (!periodicTwin_.empty() && periodicTwin_[f] >= 0) {
                const meshing::CellId tw = F.owner[periodicTwin_[f]];
                const double dx = F.centroidX[f] - C.centroidX[o];
                const double dy = F.centroidY[f] - C.centroidY[o];
                const double dz = F.centroidZ[f] - C.centroidZ[o];
                const double ddotA = std::abs(dx * Ax + dy * Ay + dz * Az);
                const double coef = opt_.rho * Amag2 / std::max(denom_for_b(o) * ddotA, 1e-30);
                add_at(A_p_, static_cast<int>(o), static_cast<int>(o), coef);
                add_at(A_p_, static_cast<int>(o), static_cast<int>(tw), -coef);
            } else if (bc
                       && (bc->type == BCType::PressureOutlet
                           || bc->type == BCType::PressureInlet)) {
                const double dx = F.centroidX[f] - C.centroidX[o];
                const double dy = F.centroidY[f] - C.centroidY[o];
                const double dz = F.centroidZ[f] - C.centroidZ[o];
                const double ddotA = std::abs(dx * Ax + dy * Ay + dz * Az);
                const double coef = opt_.rho * Amag2 / std::max(denom_for_b(o) * ddotA, 1e-30);
                add_at(A_p_, o, o, coef); // Dirichlet p'=0 anchors the system
                anyDirichlet = true;
            }
            // Walls / symmetry / velocity-inlet: zero-gradient on p' (no coef).
        }
        // RHS: -divergence of the (provisional) face flux
        rhs_[o] -= phiFlux_[f];
        if (n != meshing::kBoundaryCell)
            rhs_[n] += phiFlux_[f];
    }

    // If no Dirichlet anchor (all-wall domain), pin cell 0 to remove the
    // compatible-but-singular null space.
    if (!anyDirichlet && nC > 0) {
        const int kd = find_col(A_p_, 0, 0);
        if (kd >= 0)
            A_p_.values[kd] *= 1.0e6;
        rhs_[0] = 0.0;
    }
}

// ============================================================ correction
void SimpleAlgorithm::correct_fields()
{
    auto& U = *F_.find_vector("U");
    auto& p = *F_.find_scalar("p");
    const auto& C = mesh_.cells();
    const auto& F = mesh_.faces();
    const std::size_t nC = C.size();

    // p ← p + urfP * p'
    for (std::size_t c = 0; c < nC; ++c)
        p[c] += opt_.urfP * pPrime_[c];

    // U ← U - (V/aP) * ∇p'
    LeastSquaresGradient G(mesh_);
    VectorField gPp;
    G.evaluate(pPrime_, gPp);
    for (std::size_t c = 0; c < nC; ++c) {
        const double w = C.volume[c] / std::max(aP_[c], 1e-30);
        U.x[c] -= w * gPp.x[c];
        U.y[c] -= w * gPp.y[c];
        U.z[c] -= w * gPp.z[c];
    }

    // Face-flux correction: F ← F + coef * (p'_o - p'_n)
    for (std::size_t f = 0; f < F.size(); ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
        const double Amag2 = Ax * Ax + Ay * Ay + Az * Az;
        if (n != meshing::kBoundaryCell) {
            const double dx = C.centroidX[n] - C.centroidX[o];
            const double dy = C.centroidY[n] - C.centroidY[o];
            const double dz = C.centroidZ[n] - C.centroidZ[o];
            const double ddotA = std::abs(dx * Ax + dy * Ay + dz * Az);
            const double aPf = 0.5 * (aP_[o] + aP_[n]);
            const double coef = opt_.rho * Amag2 / std::max(aPf * ddotA, 1e-30);
            phiFlux_[f] += coef * (pPrime_[o] - pPrime_[n]);
        }
    }
}

// ============================================================ helper: solve one momentum component
double SimpleAlgorithm::solve_momentum_component(int k)
{
    auto& U = *F_.find_vector("U");
    util::aligned_vector<double>* Ucomp[3] = {&U.x, &U.y, &U.z};
    sol_ = *Ucomp[k];
    linMom_.solve(A_mom_, rhs_, sol_);
    util::aligned_vector<double> Ax(sol_.size(), 0);
    A_mom_.spmv(sol_, Ax);
    double s = 0;
    for (std::size_t i = 0; i < sol_.size(); ++i) {
        const double d = Ax[i] - rhs_[i];
        s += d * d;
    }
    *Ucomp[k] = std::move(sol_);
    return std::sqrt(s);
}

// ============================================================ helper: solve p' & correct
double SimpleAlgorithm::solve_pressure_correction_and_correct()
{
    pPrime_.assign(mesh_.cells().size(), 0.0);
    linP_.solve(A_p_, rhs_, pPrime_);
    double cs = 0;
    for (double v : rhs_)
        cs += v * v;
    correct_fields();
    return std::sqrt(cs);
}

// ============================================================ iterate
SimpleResiduals SimpleAlgorithm::iterate()
{
    if (!sparsity_built_)
        build_sparsity();
    SimpleResiduals res{};

    compute_face_fluxes();
    for (int k = 0; k < 3; ++k) {
        assemble_momentum(k);
        res.mom[k] = solve_momentum_component(k);
    }
    // Halo exchange of velocity components so that downstream face-flux
    // and gradient computations see consistent neighbour-rank data.
    if (sync_) {
        if (auto* U = F_.find_vector("U"))
            sync_->sync_vector(U->x, U->y, U->z);
    }
    compute_face_fluxes();
    assemble_pressure_correction();
    res.cont = solve_pressure_correction_and_correct();
    // Halo exchange of corrected pressure / velocity for next outer iter.
    if (sync_) {
        if (auto* p = F_.find_scalar("p"))
            sync_->sync_scalar(*p);
        if (auto* U = F_.find_vector("U"))
            sync_->sync_vector(U->x, U->y, U->z);
    }
    return res;
}

// ============================================================ advance_time_step
SimpleResiduals SimpleAlgorithm::advance_time_step(int nInnerIters)
{
    if (!sparsity_built_)
        build_sparsity();
    // Roll history: U^{n-1} <- U^n; U^n <- current U.
    auto& U = *F_.find_vector("U");
    if (timeStep_ >= 1) {
        Ux_nm1_ = Ux_n_;
        Uy_nm1_ = Uy_n_;
        Uz_nm1_ = Uz_n_;
    }
    Ux_n_ = U.x;
    Uy_n_ = U.y;
    Uz_n_ = U.z;

    SimpleResiduals last{};
    for (int k = 0; k < std::max(1, nInnerIters); ++k)
        last = iterate();
    ++timeStep_;
    return last;
}

} // namespace simall::solver
