// =============================================================================
// SimAll Beta -- Parametric CAD / Sketcher
// File   : src/cad_sketch/ConstraintSolver.hpp
// Phase  : 23 Pass 23.1
//
// Constraint solver for the parametric sketcher.
//
// Algorithm:  damped Gauss-Newton (Levenberg-Marquardt).
//   1. Pack all *non-fixed* parameters of the sketch into x  in R^n.
//   2. For every constraint compute residual vector r(x) in R^m and analytic
//      (or finite-difference fallback) Jacobian J(x) in R^(m x n).
//   3. Solve the damped normal equations  (J^T J + lambda I) dx = -J^T r
//      using Eigen's column-pivoting QR, which simultaneously gives us
//      a numerical rank estimate => DOF and degeneracy diagnostics.
//   4. Adapt lambda Levenberg-style: shrink on accepted steps, grow on
//      rejected ones (residual norm increased).  Terminate when residual
//      norm < tol, or when step length below tol, or max-iters hit.
//
// Output: a `SolveReport` carrying status, iteration count, final residual
// norm, residual history, computed rank/DOF, and a list of parameter ids
// the QR pivoting flagged as rank-deficient (the degeneracy reporter).
// Callers can then highlight over- or under-constrained sketches in the UI.
// =============================================================================
#pragma once

#include "cad_sketch/Sketch.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace simall::cad::sketch {

enum class SolveStatus : std::uint8_t {
    Converged           = 0,   // residual norm below tolerance
    StalledMaxIters     = 1,   // ran out of iterations
    StalledStepTooSmall = 2,   // step length below tolerance but residual still > tol
    Underconstrained    = 3,   // rank(J) < unknowns -- multiple solutions
    Overconstrained     = 4,   // rank(J) == unknowns but residual cannot be driven to 0
    Singular            = 5,   // numerical breakdown (NaN / Inf encountered)
    Empty               = 6    // no constraints (or no unknowns) to solve
};

struct SolveReport {
    SolveStatus              status         {SolveStatus::Empty};
    std::size_t              iterations     {0};
    double                   residual_norm  {0.0};
    std::vector<double>      residual_history{};  // ||r||_2 each iteration
    std::size_t              unknowns       {0};  // n
    std::size_t              residual_count {0};  // m
    std::size_t              rank           {0};  // numerical rank of J at final x
    std::size_t              dof            {0};  // unknowns - rank
    std::vector<ParameterId> degenerate_parameters{};  // pivot-rejected columns
    std::string              message        {};

    [[nodiscard]] bool ok() const noexcept {
        return status == SolveStatus::Converged;
    }
};

struct SolverConfig {
    std::size_t max_iterations    {100};
    double      residual_tol      {1e-10};
    double      step_tol          {1e-12};
    double      initial_damping   {1e-3};   // initial Levenberg lambda
    double      damping_grow      {10.0};   // on rejected step
    double      damping_shrink    {0.1};    // on accepted step
    double      rank_tol          {1e-9};   // threshold for QR rank estimate
    bool        use_finite_diff_jacobian{false};  // for testing / debugging
    double      fd_eps            {1e-7};
};

// -----------------------------------------------------------------------------
// Solves the sketch in-place: mutates the parameter values of `s` so that
// constraint residuals are (approximately) satisfied.  Fixed parameters are
// never modified.  Returns a SolveReport describing the outcome.
// -----------------------------------------------------------------------------
[[nodiscard]] SolveReport solve(Sketch& s, const SolverConfig& cfg = {});

// -----------------------------------------------------------------------------
// Convenience: evaluate the current residual vector of a sketch without
// running any iterations.  Useful for unit tests and dirty/clean diagnostics.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<double> evaluate_residuals(const Sketch& s);

}  // namespace simall::cad::sketch
