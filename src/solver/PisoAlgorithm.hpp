// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/PisoAlgorithm.hpp
// Phase  : 6.5 — Pressure-Implicit with Splitting of Operators
//
// PISO (Issa 1985) is a non-iterative variant of pressure-velocity coupling
// designed for transient flows. Each physical time step performs:
//
//   1) Momentum predictor:  solve A · U* = b(p^n)   (no under-relaxation)
//   2) Pressure correction #1:  solve continuity for p'_1, correct U*→U** & p
//   3) Pressure correction #2 (and optional #3): solve again with the
//      latest U to capture neighbour-coefficient contributions that the
//      first corrector treated explicitly.
//
// Number of correctors:  typically 2 for first-order, 3 for second-order in
// time. We default to 2 and expose pisoCorrectors.
//
// Implemented as a SimpleAlgorithm subclass that overrides iterate() —
// reuses build_sparsity, assemble_momentum, assemble_pressure_correction,
// correct_fields verbatim.
//
// References:
//   Issa, J. Comp. Phys. 62, 40-65 (1985).
//   Ferziger & Perić, §7.5.4.
// =============================================================================
#pragma once

#include "solver/SimpleAlgorithm.hpp"

namespace simall::solver {

struct PisoOptions {
    int    nCorrectors  = 2;     // number of pressure correctors (≥1)
    int    nNonOrthCorr = 0;     // extra non-orthogonal pressure passes per corrector
    double rho          = 1.0;
    double mu           = 1.0e-3;
    double dt           = 1.0e-3;   // PISO assumes transient (dt > 0)
    TemporalScheme timeScheme = TemporalScheme::BDF2;
};

class PisoAlgorithm : public SimpleAlgorithm {
public:
    PisoAlgorithm(meshing::Mesh& mesh,
                  FieldRegistry& fields,
                  const std::vector<BoundarySpec>& boundaries,
                  ILinearSolver& momentumSolver,
                  ILinearSolver& pressureSolver,
                  PisoOptions opts);

    /// One PISO outer iteration = one predictor + nCorrectors correctors.
    /// For transient analysis call advance_time_step(1) per Δt — PISO is
    /// non-iterative by design.
    SimpleResiduals iterate() override;

    int n_correctors() const noexcept { return piso_.nCorrectors; }

private:
    PisoOptions piso_;

    static SimpleOptions to_simple_opts(const PisoOptions& p) {
        SimpleOptions s;
        s.urfU = 1.0;                 // PISO uses no momentum under-relaxation
        s.urfP = 1.0;                 // … and no pressure under-relaxation
        s.rho  = p.rho;
        s.mu   = p.mu;
        s.dt   = p.dt;
        s.timeScheme = p.timeScheme;
        s.algorithm  = PvCouplingVariant::SIMPLE;  // standard p' denom = 1/aP
        return s;
    }
};

}  // namespace simall::solver
