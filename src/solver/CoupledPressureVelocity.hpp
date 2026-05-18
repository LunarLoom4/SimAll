// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/CoupledPressureVelocity.hpp
// Phase  : 6.5 — Fully-coupled (monolithic) block solve of momentum +
//                continuity in a single 4N × 4N system.
//
// Unknown ordering (interleaved 4-per-cell, cell-major):
//   row index for cell c, equation eq ∈ {0:u, 1:v, 2:w, 3:p}  →  4·c + eq
//
// Equation rows:
//   eq=0:  momentum-x:  Σ_NB a_NB · u_NB + (a_P/urfU) u_P + ∂p/∂x · V = b_u
//   eq=1:  momentum-y:  …                                              = b_v
//   eq=2:  momentum-z:  …                                              = b_w
//   eq=3:  continuity:  Σ_faces ρ · (V/aP)_f · ∇p|_f · A_f
//                       - Σ_faces ρ · ⟨U⟩_f · A_f = 0
//
// The pressure-gradient coupling in momentum rows and the velocity-divergence
// coupling in continuity rows make the system non-symmetric but block-diag-
// dominant. We solve with the existing ILinearSolver (GMRES/BiCGSTAB/TFQMR).
//
// Advantages vs SIMPLE/PISO:
//   - Unconditionally stable; no urf needed; converges in O(10) outer iters
//     for steady flows where SIMPLE needs O(10²-10³).
//   - Robust on highly-skewed meshes and at low Reynolds.
// Costs:
//   - 4×4 block per cell pair → memory ≈ 16× SIMPLE momentum CSR.
//   - Linear solve much more expensive per iteration.
//
// References:
//   Mazhar (2016) "A coupled pressure-velocity solver for OpenFOAM";
//   Darwish & Moukalled, "The Finite Volume Method in CFD" §15.3.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "FieldRegistry.hpp"
#include "solver/Solver.hpp"
#include "solver/CSRMatrix.hpp"
#include "solver/LinearSolvers.hpp"

namespace simall::solver {

struct CoupledPvOptions {
    double rho        = 1.0;
    double mu         = 1.0e-3;
    double dt         = 0.0;       // 0 → steady
    TemporalScheme    timeScheme = TemporalScheme::ImplicitEuler;
    // Pressure-equation Dirichlet anchor: if no pressure-outlet BC exists,
    // cell 0's continuity row is replaced by p_0 = 0 to remove the null space.
    bool   anchorPressureIfNoOutlet = true;
};

struct CoupledPvResiduals {
    double mom[3] = {0,0,0};
    double cont   = 0.0;
};

class CoupledPressureVelocity {
public:
    CoupledPressureVelocity(meshing::Mesh& mesh,
                            FieldRegistry& fields,
                            const std::vector<BoundarySpec>& boundaries,
                            ILinearSolver& blockSolver,
                            CoupledPvOptions opts);

    /// One outer (Picard) iteration: assemble + monolithic solve + apply.
    CoupledPvResiduals iterate();

    /// Advance one physical time step (BDF2 if opts.timeScheme=BDF2 and
    /// step≥1). Internally rolls U^{n-1}, U^n history.
    CoupledPvResiduals advance_time_step(int nOuter = 1);

private:
    void build_sparsity();
    void assemble();          // fills A_block_ and rhs_
    void apply_solution(const util::aligned_vector<double>& x,
                        CoupledPvResiduals& res);

    meshing::Mesh&                   mesh_;
    FieldRegistry&                   F_;
    const std::vector<BoundarySpec>& bcs_;
    ILinearSolver&                   lin_;
    CoupledPvOptions                 opt_;

    CSRMatrix                    A_block_;
    util::aligned_vector<double> rhs_;
    util::aligned_vector<double> x_;
    util::aligned_vector<double> aP_mom_;   // momentum diagonal (for Rhie-Chow weights)

    // BDF2 history
    util::aligned_vector<double> Ux_n_, Uy_n_, Uz_n_;
    util::aligned_vector<double> Ux_nm1_, Uy_nm1_, Uz_nm1_;
    int timeStep_ = 0;

    bool sparsity_built_ = false;
};

}  // namespace simall::solver
