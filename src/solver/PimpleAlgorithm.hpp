// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/PimpleAlgorithm.hpp
// Phase  : 22 Pass 1 — Foundation / PIMPLE (hybrid PISO-inside-SIMPLE).
//
// PIMPLE combines PISO's non-iterative inner pressure correctors with
// SIMPLE's outer iteration. Each physical time step performs:
//
//   for o = 1 .. nOuterCorrectors:                          // SIMPLE outer
//       if momentumPredictor:
//           assemble & solve momentum     (with urfU under-relaxation)
//       compute_face_fluxes()
//       for c = 1 .. nCorrectors:                            // PISO inner
//           for n = 0 .. nNonOrthCorr:                       // non-orth pass
//               assemble_pressure_correction()
//               solve_pressure_correction_and_correct()      (with urfP)
//       check residuals < tolU / tolP  → early exit
//
// Degenerate cases:
//   nOuterCorrectors == 1                                 → standard PISO
//   nCorrectors      == 1, nOuterCorrectors → many        → standard SIMPLE
//   momentumPredictor == false, nOuterCorrectors == 1     → PISO-no-predictor
//
// References:
//   - Issa, J. Comp. Phys. 62, 40-65 (1985).        [PISO]
//   - Patankar, "Numerical Heat Transfer" (1980).   [SIMPLE]
//   - Jasak PhD §3.8 (1996).                        [PIMPLE in OpenFOAM]
//   - OpenFOAM User Guide §4.3 (pimpleFoam).
// =============================================================================
#pragma once

#include "solver/SimpleAlgorithm.hpp"

namespace simall::solver {

struct PimpleOptions {
    int    nOuterCorrectors = 2;      // SIMPLE outer loop count (≥1)
    int    nCorrectors      = 2;      // PISO inner pressure correctors (≥1)
    int    nNonOrthCorr     = 0;      // extra non-orthogonal pressure passes
    bool   momentumPredictor = true;  // skip predictor when false (PISO-style)

    // Under-relaxation only applied when nOuterCorrectors > 1. When ==1
    // PIMPLE behaves as pure PISO (urf forced to 1.0 internally).
    double urfU             = 0.7;
    double urfP             = 0.3;

    // Outer-loop early-exit tolerances (compared against
    // SimpleResiduals from the most recent outer iteration). 0 disables.
    double tolU             = 0.0;
    double tolP             = 0.0;

    // Physical parameters (forwarded to SimpleOptions).
    double rho              = 1.0;
    double mu               = 1.0e-3;
    double dt               = 1.0e-3;   // transient — PIMPLE wants dt > 0
    TemporalScheme timeScheme = TemporalScheme::BDF2;
};

class PimpleAlgorithm : public SimpleAlgorithm {
public:
    PimpleAlgorithm(meshing::Mesh& mesh,
                    FieldRegistry& fields,
                    const std::vector<BoundarySpec>& boundaries,
                    ILinearSolver& momentumSolver,
                    ILinearSolver& pressureSolver,
                    PimpleOptions opts);

    /// One PIMPLE iteration call performs the full nOuterCorrectors-deep
    /// outer loop. Returns residuals from the FINAL outer corrector.
    /// (For driver compatibility with SimpleAlgorithm::advance_time_step,
    /// pass nInnerIters=1 there — the multi-corrector behaviour is now
    /// encapsulated inside iterate().)
    SimpleResiduals iterate() override;

    int  n_outer_correctors() const noexcept { return pimple_.nOuterCorrectors; }
    int  n_correctors()       const noexcept { return pimple_.nCorrectors;      }
    int  n_non_orth_corr()    const noexcept { return pimple_.nNonOrthCorr;     }
    bool momentum_predictor() const noexcept { return pimple_.momentumPredictor;}

    /// Returns the outer-corrector index at which the most recent iterate()
    /// invocation exited (0-based; equals nOuterCorrectors-1 if no early
    /// exit). Useful for monitoring outer-loop saturation.
    int last_outer_iters() const noexcept { return lastOuterIters_; }

private:
    PimpleOptions pimple_;
    int           lastOuterIters_ = 0;

    static SimpleOptions to_simple_opts(const PimpleOptions& p) {
        SimpleOptions s;
        // When nOuterCorrectors==1 PIMPLE ≡ PISO: turn off under-relaxation.
        const bool pisoMode = (p.nOuterCorrectors <= 1);
        s.urfU = pisoMode ? 1.0 : p.urfU;
        s.urfP = pisoMode ? 1.0 : p.urfP;
        s.rho  = p.rho;
        s.mu   = p.mu;
        s.dt   = p.dt;
        s.timeScheme = p.timeScheme;
        s.algorithm  = PvCouplingVariant::SIMPLE;  // standard p' denom
        return s;
    }
};

}  // namespace simall::solver
