// =============================================================================
// SimAll Beta - Optimization Subsystem
// File   : src/optimization/CmaEs.hpp
// Week   : 18
//
// Covariance Matrix Adaptation Evolution Strategy (Hansen 2003).  Black-
// box optimisation of a (possibly noisy) scalar function:
//
//     min  f(x)        x ∈ R^n,    box-bounded
//
// CMA-ES samples a population of λ candidate vectors from a multivariate
// Gaussian N(m, σ²C), evaluates f, picks the best μ, then updates the
// mean m, step size σ, and covariance C using the well-known evolution
// path heuristics.  Default hyperparameters match the original paper.
//
// Industrial-CFD use: gradient-free shape & operating-point optimisation,
// design-of-experiments completion, kriging-surrogate inner loop.
// =============================================================================
#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace simall::optimization {

struct CmaEsOptions {
    std::size_t        dim         = 0;
    std::size_t        maxGen      = 200;
    std::size_t        lambda      = 0;   // 0 → 4 + ⌊3·ln n⌋
    double             sigma0      = 0.3;
    double             tolFitness  = 1e-10;
    double             tolSigma    = 1e-12;
    std::uint64_t      seed        = 0xC0FFEE;
    std::vector<double> mean0;            // size = dim (defaults to 0)
    std::vector<double> xLow;             // size = dim (defaults to -inf)
    std::vector<double> xUp;              // size = dim (defaults to +inf)
    std::function<double(const std::vector<double>&)> evaluate;
};

struct CmaEsResult {
    bool                  ok = false;
    std::string           error;
    std::size_t           generations = 0;
    double                bestF = 0.0;
    std::vector<double>   bestX;
    std::vector<double>   history;   // best f per generation
};

[[nodiscard]] CmaEsResult run_cmaes(CmaEsOptions opt);

}  // namespace simall::optimization
