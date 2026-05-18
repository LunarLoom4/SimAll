// =============================================================================
// SimAll Beta - Optimization Subsystem
// File   : src/optimization/MmaOptimizer.hpp
// Week   : 18
//
// Method of Moving Asymptotes (Svanberg, 1987 / 2002).  Solves the
// inequality-constrained NLP
//
//     min   f(x)              x ∈ R^n
//     s.t.  g_j(x) ≤ 0        j = 1..m
//           x_low ≤ x ≤ x_up
//
// by sequential convex approximations using moving asymptotes L_i, U_i.
// The asymptotes adapt to the iterate history: when consecutive moves of
// x_i agree in sign the asymptotes expand (s = 1.2); when they oscillate
// they contract (s = 0.7).  This is the canonical "industrial MMA" used
// by TOSCA, OptiStruct, Tosca-Flow.
//
// This implementation:
//   * solves the convex MMA subproblem with a simple primal projected
//     gradient (sufficient for ≤500-DV problems typical of FFD shape opt);
//   * applies the original Svanberg defaults for asymptote update;
//   * never requires user-supplied second derivatives.
// =============================================================================
#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace simall::optimization {

struct MmaProblem {
    std::size_t n = 0;                 // design variables
    std::size_t m = 0;                 // inequality constraints
    std::vector<double> xLow;
    std::vector<double> xUp;
    std::vector<double> x0;            // initial guess

    /// Evaluator: fills f, g (size m), df (size n), dg (size m*n row-major).
    std::function<void(const std::vector<double>& x,
                       double&                    f,
                       std::vector<double>&       g,
                       std::vector<double>&       df,
                       std::vector<double>&       dg)> evaluate;
};

struct MmaOptions {
    std::size_t maxIter   = 100;
    double      moveLimit = 0.2;          // fraction of (xUp - xLow)
    double      tolX      = 1e-6;
    double      asyInit   = 0.5;
    double      asyIncr   = 1.2;
    double      asyDecr   = 0.7;
    std::size_t subIters  = 60;
};

struct MmaResult {
    bool                ok = false;
    std::string         error;
    std::size_t         iterations = 0;
    double              f         = 0.0;
    std::vector<double> x;
    std::vector<double> g;
    std::vector<double> history;        // f at each outer iter
};

[[nodiscard]] MmaResult run_mma(MmaProblem  problem,
                                 MmaOptions  opt = {});

}  // namespace simall::optimization
