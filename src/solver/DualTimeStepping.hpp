// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/DualTimeStepping.hpp
// Phase  : 6.5 — Pseudo-time (dual time-stepping) wrapper.
//
// For each physical time step Δt, the underlying coupling algorithm
// (SIMPLE/SIMPLEC/PISO/CoupledPV) is iterated in a pseudo-time τ until the
// pseudo-time residual ∂φ/∂τ drops below a tolerance.  This converts an
// inherently iterative steady scheme into a time-accurate one (for SIMPLE/
// SIMPLEC; PISO and Coupled are already time-accurate, but dual-time-stepping
// still benefits them through stiffness damping).
//
// Convergence is tracked through the maximum reduction factor of momentum
// and continuity residuals between successive pseudo-iterations:
//
//      ρ_k = max_i (r_i^{k+1} / r_i^k)
//
// We stop when r_i < absTol or k ≥ maxPseudoIters.
//
// Reference:  Jameson, "Time-Dependent Calculations Using Multigrid…"
//             AIAA 91-1596 (1991);  Pulliam & Steger (1985).
// =============================================================================
#pragma once

#include "solver/SimpleAlgorithm.hpp"

#include <functional>
#include <vector>

namespace simall::solver
{

struct DualTimeOptions
{
    int maxPseudoIters = 50;
    double absTol = 1.0e-6; // absolute convergence on max(mom, cont)
    double relTol = 1.0e-3; // relative drop from first pseudo-iter
    bool verbose = false;
};

struct DualTimeReport
{
    int iters = 0;
    double finalMom = 0.0;
    double finalCont = 0.0;
    bool converged = false;
};

class DualTimeStepping
{
public:
    DualTimeStepping(SimpleAlgorithm& alg, DualTimeOptions opts = {}) : alg_(alg), opt_(opts) {}

    /// Advance ONE physical time step using pseudo-time sub-iterations.
    /// Internally calls alg_.iterate() until convergence; finalises by
    /// calling alg_.advance_time_step(0) — no, that would double the history
    /// roll, so we roll history *here* and then iterate, mirroring what
    /// SimpleAlgorithm::advance_time_step does, but with convergence-driven
    /// inner-iteration count.
    DualTimeReport step();

    /// Run pseudo-time iterations to convergence for STEADY problems
    /// (Δt = 0 in the underlying algorithm). Does not roll the BDF2 history.
    DualTimeReport solve_steady();

private:
    SimpleAlgorithm& alg_;
    DualTimeOptions opt_;
};

} // namespace simall::solver
