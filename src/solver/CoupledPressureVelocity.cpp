// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/CoupledPressureVelocity.cpp
//
// Monolithic 4N × 4N coupled solve.  Sparsity per cell:
//   4 rows × (4 × (1 + |neighbours|)) columns each.
// =============================================================================
#include "solver/CoupledPressureVelocity.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace simall::solver {

namespace {
const BoundarySpec* find_bc(const std::vector<BoundarySpec>& bcs,
                            meshing::ZoneId z) {
    for (const auto& b : bcs) if (b.zone == z) return &b;
    return nullptr;
}
inline int find_col(const CSRMatrix& A, int row, int col) {
    for (int k = A.rowPtr[row]; k < A.rowPtr[row + 1]; ++k)
        if (A.colIdx[k] == col) return k;
    return -1;
}
inline void add_at(CSRMatrix& A, int row, int col, double v) {
    const int k = find_col(A, row, col);
    if (k >= 0) A.values[k] += v;
}
}  // namespace

CoupledPressureVelocity::CoupledPressureVelocity(
        meshing::Mesh& mesh,
        FieldRegistry& fields,
        const std::vector<BoundarySpec>& boundaries,
        ILinearSolver& blockSolver,
        CoupledPvOptions opts)
    : mesh_(mesh), F_(fields), bcs_(boundaries),
      lin_(blockSolver), opt_(opts) {}

// ============================================================ sparsity
void CoupledPressureVelocity::build_sparsity() {
    const auto& C = mesh_.cells();
    const auto& Fa = mesh_.faces();
    const std::size_t nC = C.size();
    const std::size_t N  = 4 * nC;

    A_block_.rowPtr.assign(N + 1, 0);
    A_block_.colIdx.clear();

    // Build the list of cell-neighbours for each cell once.
    std::vector<std::vector<int>> nbrs(nC);
    for (std::size_t c = 0; c < nC; ++c) {
        nbrs[c].push_back(static_cast<int>(c));
        const int beg = C.faceOffsets[c], end = C.faceOffsets[c + 1];
        for (int k = beg; k < end; ++k) {
            const meshing::FaceId fid = C.faceIndices[k];
            const meshing::CellId oth =
                (Fa.owner[fid] == c) ? Fa.neighbor[fid] : Fa.owner[fid];
            if (oth != meshing::kBoundaryCell)
                nbrs[c].push_back(static_cast<int>(oth));
        }
        std::sort(nbrs[c].begin(), nbrs[c].end());
        nbrs[c].erase(std::unique(nbrs[c].begin(), nbrs[c].end()), nbrs[c].end());
    }

    // For each cell c and equation eq, the row 4c+eq has columns
    // { 4*n + e for n in nbrs(c), e in {0,1,2,3} } — full 4-block coupling
    // for every (c, n) pair.
    for (std::size_t c = 0; c < nC; ++c) {
        for (int eq = 0; eq < 4; ++eq) {
            const int row = 4 * static_cast<int>(c) + eq;
            for (int n : nbrs[c])
                for (int e = 0; e < 4; ++e)
                    A_block_.colIdx.push_back(4 * n + e);
            A_block_.rowPtr[row + 1] = static_cast<int>(A_block_.colIdx.size());
        }
    }
    A_block_.values.assign(A_block_.colIdx.size(), 0.0);
    rhs_.assign(N, 0.0);
    x_  .assign(N, 0.0);
    aP_mom_.assign(nC, 0.0);

    Ux_n_  .assign(nC, 0.0); Uy_n_  .assign(nC, 0.0); Uz_n_  .assign(nC, 0.0);
    Ux_nm1_.assign(nC, 0.0); Uy_nm1_.assign(nC, 0.0); Uz_nm1_.assign(nC, 0.0);

    sparsity_built_ = true;
}

// ============================================================ assemble
void CoupledPressureVelocity::assemble() {
    const auto& C = mesh_.cells();
    const auto& Fa = mesh_.faces();
    const std::size_t nC = C.size();
    auto& U = *F_.find_vector("U");
    auto& p = *F_.find_scalar("p");

    std::fill(A_block_.values.begin(), A_block_.values.end(), 0.0);
    std::fill(rhs_.begin(), rhs_.end(), 0.0);
    std::fill(aP_mom_.begin(), aP_mom_.end(), 0.0);

    // ----- helper to write into block (cell row, cell col, eqRow, eqCol)
    auto add_block = [&](int rowCell, int colCell, int eqR, int eqC, double v) {
        const int row = 4 * rowCell + eqR;
        const int col = 4 * colCell + eqC;
        add_at(A_block_, row, col, v);
    };

    bool anyPressureDirichlet = false;

    // ----- face loop: convection + diffusion + pressure-gradient + continuity
    for (std::size_t f = 0; f < Fa.size(); ++f) {
        const meshing::CellId o = Fa.owner[f];
        const meshing::CellId n = Fa.neighbor[f];
        const double Ax = Fa.areaX[f], Ay = Fa.areaY[f], Az = Fa.areaZ[f];
        const double Amag2 = Ax*Ax + Ay*Ay + Az*Az;

        if (n != meshing::kBoundaryCell) {
            // ---- geometry ----
            const double dx = C.centroidX[n] - C.centroidX[o];
            const double dy = C.centroidY[n] - C.centroidY[o];
            const double dz = C.centroidZ[n] - C.centroidZ[o];
            const double ddotA = std::abs(dx*Ax + dy*Ay + dz*Az);
            const double D = opt_.mu * Amag2 / std::max(ddotA, 1.0e-30);

            // ---- mass flux estimate (Picard linearisation: use current U)
            const double Uxf = 0.5 * (U.x[o] + U.x[n]);
            const double Uyf = 0.5 * (U.y[o] + U.y[n]);
            const double Uzf = 0.5 * (U.z[o] + U.z[n]);
            const double Ff  = opt_.rho * (Uxf*Ax + Uyf*Ay + Uzf*Az);
            const double ConvO = std::max(Ff, 0.0);
            const double ConvN = std::max(-Ff, 0.0);

            // ---- momentum (3 equations × upwind convection + diffusion) ----
            for (int eq = 0; eq < 3; ++eq) {
                add_block(o, o, eq, eq,  (ConvO + D));
                add_block(o, n, eq, eq, -(ConvN + D));
                add_block(n, n, eq, eq,  (ConvN + D));
                add_block(n, o, eq, eq, -(ConvO + D));
            }

            // ---- pressure-gradient coupling in momentum rows: ----
            // ∂p/∂x_eq · V = (p_n - p_o) * A_eq  (in coupled form, distributed
            // half-half between cells). We use a finite-volume face integral:
            //   ∫ p n̂ dA ≈ ½(p_o + p_n) · A_eq      (centred — explicit non-symm)
            // which becomes in matrix form:
            //   add  +½ A_eq to (row=mom_eq of o, col=p of o), (col=p of n)
            //   add  -½ A_eq to (row=mom_eq of n, col=p of o), (col=p of n)
            const double Acomp[3] = { Ax, Ay, Az };
            for (int eq = 0; eq < 3; ++eq) {
                add_block(o, o, eq, 3,  0.5 * Acomp[eq]);
                add_block(o, n, eq, 3,  0.5 * Acomp[eq]);
                add_block(n, o, eq, 3, -0.5 * Acomp[eq]);
                add_block(n, n, eq, 3, -0.5 * Acomp[eq]);
            }

            // ---- continuity coupling: divergence of velocity ----
            // ∫ ρ U · n̂ dA ≈ ρ · ½(U_o + U_n) · A
            // → in matrix form (row=p of o):  +½ ρ A_eq · u_eq{o,n}
            //                    (row=p of n):  -½ ρ A_eq · u_eq{o,n}
            for (int eq = 0; eq < 3; ++eq) {
                add_block(o, o, 3, eq,  0.5 * opt_.rho * Acomp[eq]);
                add_block(o, n, 3, eq,  0.5 * opt_.rho * Acomp[eq]);
                add_block(n, o, 3, eq, -0.5 * opt_.rho * Acomp[eq]);
                add_block(n, n, 3, eq, -0.5 * opt_.rho * Acomp[eq]);
            }

            // ---- continuity pressure-Laplacian stabilisation (Rhie-Chow) ----
            // To avoid odd-even decoupling of pressure on collocated grids
            // we add a small Laplacian on the pressure block in the
            // continuity equation:
            //    -(ρ V/aP)_f · (p_n - p_o)/|d_ON|² · |A|² · sign
            // For initial Picard iterations aP_mom_ is zero; we accumulate
            // it below from the momentum diagonal we just added.
            // The actual stabilising coefficient is filled in a 2nd pass
            // after aP_mom_ is finalised (see below).
        } else {
            // ---- boundary contributions ----
            const auto* bc = find_bc(bcs_, Fa.boundaryZone[f]);
            const double dx = Fa.centroidX[f] - C.centroidX[o];
            const double dy = Fa.centroidY[f] - C.centroidY[o];
            const double dz = Fa.centroidZ[f] - C.centroidZ[o];
            const double ddotA = std::abs(dx*Ax + dy*Ay + dz*Az);
            const double D = opt_.mu * Amag2 / std::max(ddotA, 1.0e-30);
            const double Acomp[3] = { Ax, Ay, Az };

            if (bc && (bc->type == BCType::NoSlipWall ||
                       bc->type == BCType::Wall ||
                       bc->type == BCType::VelocityInlet ||
                       bc->type == BCType::MovingWall)) {
                // Dirichlet velocity (zero for walls; specified for inlet)
                for (int eq = 0; eq < 3; ++eq) {
                    add_block(o, o, eq, eq, D);   // diffusion only (upwind conv from Uf=Ub absorbed in RHS)
                    const double Ub = (bc->type == BCType::NoSlipWall ||
                                       bc->type == BCType::Wall) ? 0.0
                                      : bc->vectorValue[eq];
                    rhs_[4*o + eq] += D * Ub;
                    // No-slip walls contribute zero mass flux → no continuity coupling.
                    // Inlets contribute known velocity → known mass flux into RHS of continuity:
                    if (bc->type == BCType::VelocityInlet ||
                        bc->type == BCType::MovingWall) {
                        // continuity RHS: -ρ U_b · A   (because it appears as
                        // sum of outgoing fluxes on LHS).
                        rhs_[4*o + 3] -= opt_.rho * Ub * Acomp[eq];
                    }
                }
                // Boundary pressure-gradient term: ∫ p n̂ dA ≈ p_o · A (zero
                // grad ⇒ extrapolate from owner; the integral contributes
                // p_o · A_eq to momentum row of o).
                for (int eq = 0; eq < 3; ++eq)
                    add_block(o, o, eq, 3, Acomp[eq]);
            } else if (bc && (bc->type == BCType::PressureOutlet ||
                              bc->type == BCType::PressureInlet)) {
                // Dirichlet pressure: p_b = bc->scalarValue.
                // Momentum: zero-grad extrapolation on U → upwind conv only.
                const double Uxb = U.x[o], Uyb = U.y[o], Uzb = U.z[o];
                const double Fb = opt_.rho * (Uxb*Ax + Uyb*Ay + Uzb*Az);
                const double conv = std::max(Fb, 0.0);
                for (int eq = 0; eq < 3; ++eq)
                    add_block(o, o, eq, eq, conv);   // upwind diagonal
                // Pressure-gradient: p_b · A → moved fully to RHS of momentum.
                for (int eq = 0; eq < 3; ++eq)
                    rhs_[4*o + eq] -= bc->scalarValue * Acomp[eq];
                // Continuity: outflow conv velocity contributes to LHS as
                // ρ · U_o · A_eq fully on owner.
                for (int eq = 0; eq < 3; ++eq)
                    add_block(o, o, 3, eq, opt_.rho * Acomp[eq]);
                // Pin pressure DOF for this cell? We anchor through the
                // strict diagonal-boost path; here just flag that an outlet exists.
                anyPressureDirichlet = true;
            } else if (bc && bc->type == BCType::Symmetry) {
                // Zero normal velocity + zero tangential gradient → no LHS
                // contribution; tangential momentum is naturally satisfied.
            } else {
                // Generic outflow: zero-grad on everything (extrapolate).
                const double Uxb = U.x[o], Uyb = U.y[o], Uzb = U.z[o];
                const double Fb = opt_.rho * (Uxb*Ax + Uyb*Ay + Uzb*Az);
                const double conv = std::max(Fb, 0.0);
                for (int eq = 0; eq < 3; ++eq)
                    add_block(o, o, eq, eq, conv);
                for (int eq = 0; eq < 3; ++eq)
                    add_block(o, o, eq, 3, Acomp[eq]);
                for (int eq = 0; eq < 3; ++eq)
                    add_block(o, o, 3, eq, opt_.rho * Acomp[eq]);
            }
        }
    }

    // ---- transient terms (BDF2/Euler) on momentum rows ----
    if (opt_.dt > 0.0) {
        const util::aligned_vector<double>* Un[3]   = { &Ux_n_,   &Uy_n_,   &Uz_n_ };
        const util::aligned_vector<double>* Unm1[3] = { &Ux_nm1_, &Uy_nm1_, &Uz_nm1_ };
        const bool useBDF2 =
            (opt_.timeScheme == TemporalScheme::BDF2) && (timeStep_ >= 1);
        for (std::size_t c = 0; c < nC; ++c) {
            const double Vrho_dt = opt_.rho * C.volume[c] / opt_.dt;
            for (int eq = 0; eq < 3; ++eq) {
                const int row = 4 * c + eq;
                const int kd  = find_col(A_block_, row, row);
                if (kd < 0) continue;
                if (useBDF2) {
                    A_block_.values[kd] += 1.5 * Vrho_dt;
                    rhs_[row] += (2.0 * (*Un[eq])[c] - 0.5 * (*Unm1[eq])[c]) * Vrho_dt;
                } else {
                    A_block_.values[kd] += Vrho_dt;
                    rhs_[row] += (*Un[eq])[c] * Vrho_dt;
                }
            }
        }
    }

    // ---- Capture momentum diagonals for Rhie-Chow stabilisation ----
    for (std::size_t c = 0; c < nC; ++c) {
        // Average of u,v,w diagonals
        double sum = 0.0; int cnt = 0;
        for (int eq = 0; eq < 3; ++eq) {
            const int row = 4*c + eq;
            const int kd  = find_col(A_block_, row, row);
            if (kd >= 0) { sum += A_block_.values[kd]; ++cnt; }
        }
        aP_mom_[c] = (cnt > 0) ? (sum / cnt) : 1.0;
    }

    // ---- Rhie-Chow pressure-Laplacian stabilisation in continuity rows ----
    // Δp term:  add  +c · (p_o - p_n)  to continuity row of o
    //           with c = ρ · |A|² / (aPf · |d_ON ·  Â|)
    for (std::size_t f = 0; f < Fa.size(); ++f) {
        const meshing::CellId o = Fa.owner[f];
        const meshing::CellId n = Fa.neighbor[f];
        if (n == meshing::kBoundaryCell) continue;
        const double Ax = Fa.areaX[f], Ay = Fa.areaY[f], Az = Fa.areaZ[f];
        const double Amag2 = Ax*Ax + Ay*Ay + Az*Az;
        const double dx = C.centroidX[n] - C.centroidX[o];
        const double dy = C.centroidY[n] - C.centroidY[o];
        const double dz = C.centroidZ[n] - C.centroidZ[o];
        const double ddotA = std::abs(dx*Ax + dy*Ay + dz*Az);
        const double aPf = 0.5 * (aP_mom_[o] + aP_mom_[n]);
        const double coef = opt_.rho * Amag2 / std::max(aPf * ddotA, 1.0e-30);
        // continuity row += coef * (p_o - p_n)  with sign as in pressure-correction
        add_block(o, o, 3, 3,  coef);
        add_block(o, n, 3, 3, -coef);
        add_block(n, n, 3, 3,  coef);
        add_block(n, o, 3, 3, -coef);
    }

    // ---- Pin pressure if no Dirichlet anchor ----
    if (opt_.anchorPressureIfNoOutlet && !anyPressureDirichlet && nC > 0) {
        // Replace continuity row of cell 0 with: p_0 = 0
        const int row = 4 * 0 + 3;
        for (int k = A_block_.rowPtr[row]; k < A_block_.rowPtr[row + 1]; ++k)
            A_block_.values[k] = 0.0;
        const int kd = find_col(A_block_, row, row);
        if (kd >= 0) A_block_.values[kd] = 1.0;
        rhs_[row] = 0.0;
    }
}

// ============================================================ apply
void CoupledPressureVelocity::apply_solution(
        const util::aligned_vector<double>& x,
        CoupledPvResiduals& res) {
    auto& U = *F_.find_vector("U");
    auto& p = *F_.find_scalar("p");
    const std::size_t nC = mesh_.cells().size();

    // Residual = ||A x - b||
    util::aligned_vector<double> Ax(x.size(), 0.0);
    A_block_.spmv(x, Ax);
    double sm[3] = {0,0,0}, sc = 0;
    for (std::size_t c = 0; c < nC; ++c) {
        for (int eq = 0; eq < 3; ++eq) {
            const double d = Ax[4*c+eq] - rhs_[4*c+eq];
            sm[eq] += d*d;
        }
        const double d = Ax[4*c+3] - rhs_[4*c+3];
        sc += d*d;
    }
    for (int eq = 0; eq < 3; ++eq) res.mom[eq] = std::sqrt(sm[eq]);
    res.cont = std::sqrt(sc);

    // Copy out solution
    for (std::size_t c = 0; c < nC; ++c) {
        U.x[c] = x[4*c + 0];
        U.y[c] = x[4*c + 1];
        U.z[c] = x[4*c + 2];
        p  [c] = x[4*c + 3];
    }
}

// ============================================================ iterate
CoupledPvResiduals CoupledPressureVelocity::iterate() {
    if (!sparsity_built_) build_sparsity();
    // Initial guess = current fields (warm start)
    auto& U = *F_.find_vector("U");
    auto& p = *F_.find_scalar("p");
    const std::size_t nC = mesh_.cells().size();
    x_.assign(4 * nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c) {
        x_[4*c + 0] = U.x[c];
        x_[4*c + 1] = U.y[c];
        x_[4*c + 2] = U.z[c];
        x_[4*c + 3] = p  [c];
    }
    assemble();
    lin_.solve(A_block_, rhs_, x_);
    CoupledPvResiduals res{};
    apply_solution(x_, res);
    return res;
}

// ============================================================ advance_time_step
CoupledPvResiduals CoupledPressureVelocity::advance_time_step(int nOuter) {
    if (!sparsity_built_) build_sparsity();
    auto& U = *F_.find_vector("U");
    if (timeStep_ >= 1) {
        Ux_nm1_ = Ux_n_; Uy_nm1_ = Uy_n_; Uz_nm1_ = Uz_n_;
    }
    Ux_n_ = U.x; Uy_n_ = U.y; Uz_n_ = U.z;
    CoupledPvResiduals last{};
    for (int k = 0; k < std::max(1, nOuter); ++k) last = iterate();
    ++timeStep_;
    return last;
}

}  // namespace simall::solver
